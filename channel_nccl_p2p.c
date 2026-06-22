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
#endif

typedef struct {
  ncclComm_t comm;
  channelStream_t stream;
  int nranks;
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

int channel_nccl_context_alltoall(void *ctx_ptr, const void *sendbuf, void *recvbuf, size_t count_elems) {
#if defined(HAVE_NCCL_ALLTOALL)
  int status;
  size_t count_bytes = count_elems*16;
  channel_nccl_context *ctx = (channel_nccl_context *)ctx_ptr;
  if (ctx == NULL) return -1;
  status = channel_nccl_check(ncclAlltoAll(sendbuf, recvbuf, count_bytes,
                                           ncclUint8, ctx->comm, ctx->stream));
  if (status != 0) return status;
  channelError_t stream_status = channelStreamSynchronize(ctx->stream);
  if (stream_status != channelSuccess) return 100000 + (int)stream_status;
  return 0;
#else
  int peer;
  int status;
  const char *send_bytes = (const char *)sendbuf;
  char *recv_bytes = (char *)recvbuf;
  size_t count_bytes = count_elems*16;
  channel_nccl_context *ctx = (channel_nccl_context *)ctx_ptr;
  if (ctx == NULL) return -1;
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
  status = channel_nccl_check(ncclGroupEnd());
  if (status != 0) return status;
  channelError_t stream_status = channelStreamSynchronize(ctx->stream);
  if (stream_status != channelSuccess) return 100000 + (int)stream_status;
  return 0;
#endif
}

int channel_nccl_context_destroy(void *ctx_ptr) {
  channel_nccl_context *ctx = (channel_nccl_context *)ctx_ptr;
  if (ctx == NULL) return 0;
  channelStreamSynchronize(ctx->stream);
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
