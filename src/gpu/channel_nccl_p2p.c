#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(HAVE_NCCL) || defined(HAVE_RCCL)
#if defined(HAVE_RCCL)
#include <hip/hip_runtime_api.h>
#if defined(__has_include)
#if __has_include(<rccl/rccl.h>)
#include <rccl/rccl.h>
#elif __has_include(<nccl.h>)
#include <nccl.h>
#else
#include <rccl.h>
#endif
#else
#include <rccl/rccl.h>
#endif
#define channelStream_t hipStream_t
#define channelError_t hipError_t
#define channelStreamCreateWithFlags hipStreamCreateWithFlags
#define channelStreamNonBlocking hipStreamNonBlocking
#define channelStreamSynchronize hipStreamSynchronize
#define channelStreamDestroy hipStreamDestroy
#define channelSuccess hipSuccess
#define channelEvent_t hipEvent_t
#define channelEventCreateWithFlags hipEventCreateWithFlags
#define channelEventDisableTiming hipEventDisableTiming
#define channelEventRecord hipEventRecord
#define channelEventSynchronize hipEventSynchronize
#define channelEventDestroy hipEventDestroy
#define channelEventQuery hipEventQuery
#define channelErrorNotReady hipErrorNotReady
#else
#include <cuda_runtime_api.h>
#include <nccl.h>
#define channelStream_t cudaStream_t
#define channelError_t cudaError_t
#define channelStreamCreateWithFlags cudaStreamCreateWithFlags
#define channelStreamNonBlocking cudaStreamNonBlocking
#define channelStreamSynchronize cudaStreamSynchronize
#define channelStreamDestroy cudaStreamDestroy
#define channelSuccess cudaSuccess
#define channelEvent_t cudaEvent_t
#define channelEventCreateWithFlags cudaEventCreateWithFlags
#define channelEventDisableTiming cudaEventDisableTiming
#define channelEventRecord cudaEventRecord
#define channelEventSynchronize cudaEventSynchronize
#define channelEventDestroy cudaEventDestroy
#define channelEventQuery cudaEventQuery
#define channelErrorNotReady cudaErrorNotReady
#endif

#define CHANNEL_NCCL_SLOTS 4
typedef struct {
  const void *buf;
  size_t bytes;
  int peer;
} channel_nccl_pending_send;

typedef struct {
  ncclComm_t comm;
  channelStream_t stream;
  int nranks;
  /* Completion events of asynchronous alltoalls, one per buffer slot, so a
     wait for slot s does not also wait for a later collective on the stream. */
  channelEvent_t slot_event[CHANNEL_NCCL_SLOTS];
  int slot_pending[CHANNEL_NCCL_SLOTS];
  /* Point-to-point queue (channel_nccl_p2p_*): sends waiting for the next
     group, and one completion event ("ticket") per posted receive. */
  channel_nccl_pending_send *pending;
  int npending, cap_pending;
  channelEvent_t *ticket_event;
  int nticket, cap_ticket;
} channel_nccl_context;

static channel_nccl_context *channel_nccl_default_context = NULL;

static int channel_nccl_check(ncclResult_t result) {
  if (result == ncclSuccess) return 0;
  fprintf(stderr, "NCCL/RCCL error: %s\n", ncclGetErrorString(result));
  return (int)result;
}

int channel_nccl_get_unique_id(void *id_bytes) {
  ncclUniqueId id;
  int status = channel_nccl_check(ncclGetUniqueId(&id));
  if (status != 0) return status;
  memcpy(id_bytes, &id, sizeof(id));
  return 0;
}

int channel_nccl_context_create(int nranks, int rank, const void *id_bytes, void **ctx_out) {
  ncclUniqueId id;
  channel_nccl_context *ctx = NULL;
  if (ctx_out == NULL) return -2;
  if (*ctx_out != NULL) return 0;
  ctx = (channel_nccl_context *)malloc(sizeof(*ctx));
  if (ctx == NULL) return -3;
  memcpy(&id, id_bytes, sizeof(id));
  int status = channel_nccl_check(ncclCommInitRank(&ctx->comm, nranks, id, rank));
  if (status != 0) return status;
  channelError_t stream_status = channelStreamCreateWithFlags(&ctx->stream, channelStreamNonBlocking);
  if (stream_status != channelSuccess) return 100000 + (int)stream_status;
  ctx->nranks = nranks;
  ctx->pending = NULL;
  ctx->npending = ctx->cap_pending = 0;
  ctx->ticket_event = NULL;
  ctx->nticket = ctx->cap_ticket = 0;
  for (int s = 0; s < CHANNEL_NCCL_SLOTS; ++s) {
    ctx->slot_pending[s] = 0;
    stream_status = channelEventCreateWithFlags(&ctx->slot_event[s], channelEventDisableTiming);
    if (stream_status != channelSuccess) return 100000 + (int)stream_status;
  }
  *ctx_out = ctx;
  return 0;
}

int channel_nccl_context_sendrecv(void *ctx_ptr, const void *sendbuf, size_t send_elems, int send_peer,
                                  void *recvbuf, size_t recv_elems, int recv_peer) {
  int status;
  channel_nccl_context *ctx = (channel_nccl_context *)ctx_ptr;
  if (ctx == NULL) return -1;
  status = channel_nccl_check(ncclGroupStart());
  if (status != 0) return status;
  if (recv_peer >= 0 && recv_elems > 0) {
    status = channel_nccl_check(ncclRecv(recvbuf, recv_elems*16, ncclUint8, recv_peer,
                                         ctx->comm, ctx->stream));
    if (status != 0) return status;
  }
  if (send_peer >= 0 && send_elems > 0) {
    status = channel_nccl_check(ncclSend(sendbuf, send_elems*16, ncclUint8, send_peer,
                                         ctx->comm, ctx->stream));
    if (status != 0) return status;
  }
  status = channel_nccl_check(ncclGroupEnd());
  if (status != 0) return status;
  channelError_t stream_status = channelStreamSynchronize(ctx->stream);
  if (stream_status != channelSuccess) return 100000 + (int)stream_status;
  return 0;
}

/* Enqueue an alltoall on the context stream without waiting for it. */
static int channel_nccl_enqueue_alltoall(channel_nccl_context *ctx, const void *sendbuf, void *recvbuf,
                                         size_t count_elems) {
  int status;
  size_t count_bytes = count_elems*16;
#if defined(HAVE_NCCL_ALLTOALL)
  status = channel_nccl_check(ncclAlltoAll(sendbuf, recvbuf, count_bytes,
                                           ncclUint8, ctx->comm, ctx->stream));
  return status;
#else
  int peer;
  const char *send_bytes = (const char *)sendbuf;
  char *recv_bytes = (char *)recvbuf;
  status = channel_nccl_check(ncclGroupStart());
  if (status != 0) return status;
  for (peer = 0; peer < ctx->nranks; ++peer) {
    status = channel_nccl_check(ncclRecv(recv_bytes + (size_t)peer*count_bytes, count_bytes,
                                         ncclUint8, peer, ctx->comm, ctx->stream));
    if (status != 0) return status;
    status = channel_nccl_check(ncclSend(send_bytes + (size_t)peer*count_bytes, count_bytes,
                                         ncclUint8, peer, ctx->comm, ctx->stream));
    if (status != 0) return status;
  }
  return channel_nccl_check(ncclGroupEnd());
#endif
}

int channel_nccl_context_alltoall(void *ctx_ptr, const void *sendbuf, void *recvbuf, size_t count_elems) {
  int status;
  channel_nccl_context *ctx = (channel_nccl_context *)ctx_ptr;
  if (ctx == NULL) return -1;
  status = channel_nccl_enqueue_alltoall(ctx, sendbuf, recvbuf, count_elems);
  if (status != 0) return status;
  channelError_t stream_status = channelStreamSynchronize(ctx->stream);
  if (stream_status != channelSuccess) return 100000 + (int)stream_status;
  return 0;
}

/* Asynchronous form: enqueue, record the slot's event and return at once.
   channel_nccl_context_wait_slot blocks until that collective has finished.
   The send buffer must be complete before the call (OpenMP target regions
   and the FFT wrappers synchronise before returning). */
int channel_nccl_context_alltoall_start(void *ctx_ptr, const void *sendbuf, void *recvbuf, size_t count_elems,
                                        int slot) {
  int status;
  channel_nccl_context *ctx = (channel_nccl_context *)ctx_ptr;
  if (ctx == NULL) return -1;
  if (slot < 0 || slot >= CHANNEL_NCCL_SLOTS) return -4;
  status = channel_nccl_enqueue_alltoall(ctx, sendbuf, recvbuf, count_elems);
  if (status != 0) return status;
  channelError_t ev_status = channelEventRecord(ctx->slot_event[slot], ctx->stream);
  if (ev_status != channelSuccess) return 100000 + (int)ev_status;
  ctx->slot_pending[slot] = 1;
  return 0;
}

int channel_nccl_context_wait_slot(void *ctx_ptr, int slot) {
  channel_nccl_context *ctx = (channel_nccl_context *)ctx_ptr;
  if (ctx == NULL) return 0;
  if (slot < 0 || slot >= CHANNEL_NCCL_SLOTS) return -4;
  if (!ctx->slot_pending[slot]) return 0;
  channelError_t ev_status = channelEventSynchronize(ctx->slot_event[slot]);
  if (ev_status != channelSuccess) return 100000 + (int)ev_status;
  ctx->slot_pending[slot] = 0;
  return 0;
}

/* Unlike the alltoall above -- which NCCL only provides natively behind
   HAVE_NCCL_ALLTOALL and which we otherwise emulate with grouped send/recv --
   ncclAllGather is a native NCCL collective on every version we support.  The
   y-Schur allgather branch was hard-wired to MPI_Allgather because no allgather
   was ever exported here; this is that missing path.
   count_elems counts double complex, hence the 16 bytes each, and it is the
   per-rank CONTRIBUTION: the receive buffer holds nranks times that. */
int channel_nccl_context_allgather(void *ctx_ptr, const void *sendbuf, void *recvbuf, size_t count_elems) {
  int status;
  size_t count_bytes = count_elems*16;
  channel_nccl_context *ctx = (channel_nccl_context *)ctx_ptr;
  if (ctx == NULL) return -1;
  status = channel_nccl_check(ncclAllGather(sendbuf, recvbuf, count_bytes,
                                            ncclUint8, ctx->comm, ctx->stream));
  if (status != 0) return status;
  channelError_t stream_status = channelStreamSynchronize(ctx->stream);
  if (stream_status != channelSuccess) return 100000 + (int)stream_status;
  return 0;
}

/* ---------------------------------------------------------------------------
   Nonblocking point-to-point, the stream-ordered counterpart of MPI_Isend /
   MPI_Irecv / MPI_Test / MPI_Wait used by the pipelined y solve.

   NCCL runs a context's operations strictly in stream order and pairs the
   sends and receives between two ranks in the order they were issued, so:
   - a send is only queued here; it goes out in one group with the next
     receive posted on the same context (or at a flush), so the transfer and
     the receive for the following batch progress together;
   - a receive is posted when the caller asks (it is enqueued behind
     everything already on the stream), and an event recorded after it is
     the ticket that test/wait look at;
   - nothing synchronises the stream until channel_nccl_p2p_complete.
   Send buffers must stay untouched until then; receive buffers are ready
   once their ticket has completed.  Counts are double complex (16 bytes).
   ------------------------------------------------------------------------- */
static int channel_nccl_p2p_group(channel_nccl_context *ctx, void *recvbuf, size_t recv_bytes, int recv_peer,
                                  const void *sendbuf, size_t send_bytes, int send_peer) {
  int status;
  if (ctx->npending == 0 && recv_peer < 0 && send_peer < 0) return 0;
  status = channel_nccl_check(ncclGroupStart());
  if (status != 0) return status;
  for (int i = 0; i < ctx->npending; ++i) {
    status = channel_nccl_check(ncclSend(ctx->pending[i].buf, ctx->pending[i].bytes, ncclUint8,
                                         ctx->pending[i].peer, ctx->comm, ctx->stream));
    if (status != 0) return status;
  }
  ctx->npending = 0;
  if (send_peer >= 0) {
    status = channel_nccl_check(ncclSend(sendbuf, send_bytes, ncclUint8, send_peer, ctx->comm, ctx->stream));
    if (status != 0) return status;
  }
  if (recv_peer >= 0) {
    status = channel_nccl_check(ncclRecv(recvbuf, recv_bytes, ncclUint8, recv_peer, ctx->comm, ctx->stream));
    if (status != 0) return status;
  }
  return channel_nccl_check(ncclGroupEnd());
}

/* Records a new ticket event on the stream and returns its index. */
static int channel_nccl_p2p_ticket(channel_nccl_context *ctx, int *ticket) {
  channelError_t err;
  if (ctx->nticket == ctx->cap_ticket) {
    int cap = ctx->cap_ticket ? 2*ctx->cap_ticket : 32;
    channelEvent_t *grown = (channelEvent_t *)realloc(ctx->ticket_event, cap*sizeof(*grown));
    if (grown == NULL) return -3;
    ctx->ticket_event = grown;
    for (int t = ctx->cap_ticket; t < cap; ++t) {
      err = channelEventCreateWithFlags(&ctx->ticket_event[t], channelEventDisableTiming);
      if (err != channelSuccess) return 100000 + (int)err;
    }
    ctx->cap_ticket = cap;
  }
  err = channelEventRecord(ctx->ticket_event[ctx->nticket], ctx->stream);
  if (err != channelSuccess) return 100000 + (int)err;
  *ticket = ctx->nticket++;
  return 0;
}

int channel_nccl_p2p_isend(void *ctx_ptr, const void *buf, size_t elems, int peer) {
  channel_nccl_context *ctx = (channel_nccl_context *)ctx_ptr;
  if (ctx == NULL) return -1;
  if (ctx->npending == ctx->cap_pending) {
    int cap = ctx->cap_pending ? 2*ctx->cap_pending : 32;
    channel_nccl_pending_send *grown = (channel_nccl_pending_send *)realloc(ctx->pending, cap*sizeof(*grown));
    if (grown == NULL) return -3;
    ctx->pending = grown;
    ctx->cap_pending = cap;
  }
  ctx->pending[ctx->npending].buf = buf;
  ctx->pending[ctx->npending].bytes = elems*16;
  ctx->pending[ctx->npending].peer = peer;
  ctx->npending++;
  return 0;
}

/* CHANNEL_NCCL_P2P_GROUP=0 sends the queued messages in a group of their own
   before posting a receive, instead of together with it (experiment knob). */
static int channel_nccl_p2p_group_sends(void) {
  static int value = -1;
  if (value < 0) {
    const char *env = getenv("CHANNEL_NCCL_P2P_GROUP");
    value = (env != NULL && env[0] == '0') ? 0 : 1;
  }
  return value;
}

int channel_nccl_p2p_irecv(void *ctx_ptr, void *buf, size_t elems, int peer, int *ticket) {
  channel_nccl_context *ctx = (channel_nccl_context *)ctx_ptr;
  int status;
  if (ctx == NULL) return -1;
  if (!channel_nccl_p2p_group_sends()) {
    status = channel_nccl_p2p_group(ctx, NULL, 0, -1, NULL, 0, -1);
    if (status != 0) return status;
  }
  status = channel_nccl_p2p_group(ctx, buf, elems*16, peer, NULL, 0, -1);
  if (status != 0) return status;
  return channel_nccl_p2p_ticket(ctx, ticket);
}

/* Both directions in one group, waited for before returning (the meeting
   exchange of the twisted pipeline). */
int channel_nccl_p2p_sendrecv(void *ctx_ptr, const void *sendbuf, size_t send_elems, int send_peer,
                              void *recvbuf, size_t recv_elems, int recv_peer) {
  channel_nccl_context *ctx = (channel_nccl_context *)ctx_ptr;
  int ticket;
  if (ctx == NULL) return -1;
  int status = channel_nccl_p2p_group(ctx, recvbuf, recv_elems*16, recv_peer, sendbuf, send_elems*16, send_peer);
  if (status != 0) return status;
  status = channel_nccl_p2p_ticket(ctx, &ticket);
  if (status != 0) return status;
  channelError_t err = channelEventSynchronize(ctx->ticket_event[ticket]);
  return err == channelSuccess ? 0 : 100000 + (int)err;
}

int channel_nccl_p2p_flush(void *ctx_ptr) {
  channel_nccl_context *ctx = (channel_nccl_context *)ctx_ptr;
  if (ctx == NULL) return -1;
  return channel_nccl_p2p_group(ctx, NULL, 0, -1, NULL, 0, -1);
}

int channel_nccl_p2p_test(void *ctx_ptr, int ticket, int *done) {
  channel_nccl_context *ctx = (channel_nccl_context *)ctx_ptr;
  if (ctx == NULL || ticket < 0 || ticket >= ctx->nticket) return -4;
  channelError_t err = channelEventQuery(ctx->ticket_event[ticket]);
  if (err == channelSuccess) { *done = 1; return 0; }
  if (err == channelErrorNotReady) { *done = 0; return 0; }
  return 100000 + (int)err;
}

int channel_nccl_p2p_wait(void *ctx_ptr, int ticket) {
  channel_nccl_context *ctx = (channel_nccl_context *)ctx_ptr;
  if (ctx == NULL || ticket < 0 || ticket >= ctx->nticket) return -4;
  channelError_t err = channelEventSynchronize(ctx->ticket_event[ticket]);
  return err == channelSuccess ? 0 : 100000 + (int)err;
}

/* Opens the connections to the neighbours lo and hi (either may be -1) with
   a one-value exchange each way, and waits for it.  NCCL connects a pair of
   ranks at their first send/recv, and that blocks until both sides take part:
   a pair that first talks in the middle of the pipeline (a halo sent during
   the backward sweep but received at the end) would wait on each other.
   Every rank calls this at the same point, so the connections form here. */
int channel_nccl_p2p_connect(void *ctx_ptr, const void *sendbuf, void *recvbuf, int lo, int hi) {
  channel_nccl_context *ctx = (channel_nccl_context *)ctx_ptr;
  int status;
  if (ctx == NULL) return -1;
  status = channel_nccl_check(ncclGroupStart());
  if (status != 0) return status;
  if (lo >= 0) {
    status = channel_nccl_check(ncclSend(sendbuf, 16, ncclUint8, lo, ctx->comm, ctx->stream));
    if (status == 0) status = channel_nccl_check(ncclRecv(recvbuf, 16, ncclUint8, lo, ctx->comm, ctx->stream));
    if (status != 0) return status;
  }
  if (hi >= 0) {
    status = channel_nccl_check(ncclSend(sendbuf, 16, ncclUint8, hi, ctx->comm, ctx->stream));
    if (status == 0) status = channel_nccl_check(ncclRecv((char *)recvbuf + 16, 16, ncclUint8, hi, ctx->comm, ctx->stream));
    if (status != 0) return status;
  }
  status = channel_nccl_check(ncclGroupEnd());
  if (status != 0) return status;
  channelError_t err = channelStreamSynchronize(ctx->stream);
  return err == channelSuccess ? 0 : 100000 + (int)err;
}

/* Sends everything still queued and waits until the stream is idle; tickets
   are recycled afterwards. */
int channel_nccl_p2p_complete(void *ctx_ptr) {
  channel_nccl_context *ctx = (channel_nccl_context *)ctx_ptr;
  if (ctx == NULL) return -1;
  int status = channel_nccl_p2p_group(ctx, NULL, 0, -1, NULL, 0, -1);
  if (status != 0) return status;
  channelError_t err = channelStreamSynchronize(ctx->stream);
  if (err != channelSuccess) return 100000 + (int)err;
  ctx->nticket = 0;
  return 0;
}

int channel_nccl_context_destroy(void *ctx_ptr) {
  channel_nccl_context *ctx = (channel_nccl_context *)ctx_ptr;
  if (ctx == NULL) return 0;
  channelStreamSynchronize(ctx->stream);
  for (int s = 0; s < CHANNEL_NCCL_SLOTS; ++s) channelEventDestroy(ctx->slot_event[s]);
  for (int t = 0; t < ctx->cap_ticket; ++t) channelEventDestroy(ctx->ticket_event[t]);
  free(ctx->ticket_event);
  free(ctx->pending);
  channelStreamDestroy(ctx->stream);
  ncclCommDestroy(ctx->comm);
  free(ctx);
  return 0;
}

int channel_nccl_init(int nranks, int rank, const void *id_bytes) {
  return channel_nccl_context_create(nranks, rank, id_bytes, (void **)&channel_nccl_default_context);
}

int channel_nccl_sendrecv(const void *sendbuf, size_t send_elems, int send_peer,
                          void *recvbuf, size_t recv_elems, int recv_peer) {
  return channel_nccl_context_sendrecv(channel_nccl_default_context, sendbuf, send_elems, send_peer,
                                       recvbuf, recv_elems, recv_peer);
}

int channel_nccl_finalize(void) {
  int status = channel_nccl_context_destroy(channel_nccl_default_context);
  channel_nccl_default_context = NULL;
  return status;
}
#endif
