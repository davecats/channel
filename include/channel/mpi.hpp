#pragma once

#include "channel/types.hpp"

#include <algorithm>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <mpi.h>

namespace channel {

using MpiComm = MPI_Comm;
inline MpiComm world_comm() { return MPI_COMM_WORLD; }

int mpi_rank(MpiComm comm = world_comm());
int mpi_size(MpiComm comm = world_comm());
void mpi_barrier(MpiComm comm = world_comm());

MpiComm mpi_split(MpiComm parent, int color, int key);
void mpi_free(MpiComm& comm);

void alltoall_complex_device(const Complex* send, Complex* recv, int count, MpiComm comm, const std::string& label);
void alltoallv_complex_device(const Complex* send,
                               std::span<const int> send_counts,
                               std::span<const int> send_displs,
                               Complex* recv,
                               std::span<const int> recv_counts,
                               std::span<const int> recv_displs,
                               MpiComm comm,
                               const std::string& label);
void allgather_complex_device(const Complex* send, int send_count, Complex* recv, MpiComm comm, const std::string& label);
void allreduce_sum_complex_device(const Complex* send, Complex* recv, int count, MpiComm comm, const std::string& label);
void gatherv_complex_device(const Complex* send,
                            int send_count,
                            Complex* recv,
                            std::span<const int> recv_counts,
                            std::span<const int> recv_displs,
                            int root,
                            MpiComm comm,
                            const std::string& label);
void bcast_complex_device(Complex* values, int count, int root, MpiComm comm, const std::string& label);
void sendrecv_complex_device(const Complex* send,
                             int send_count,
                             int dest,
                             int send_tag,
                             Complex* recv,
                             int recv_count,
                             int source,
                             int recv_tag,
                             MpiComm comm,
                             const std::string& label);

template <class View>
void require_contiguous_complex_view(const View& view, std::size_t required_values, const std::string& label) {
  static_assert(std::is_same_v<typename View::non_const_value_type, Complex>,
                "MPI complex device helpers require Kokkos::View<channel::Complex...>");
  if (!view.span_is_contiguous()) {
    throw std::runtime_error("MPI buffer view is not contiguous in " + label);
  }
  if (view.span() < required_values) {
    throw std::runtime_error("MPI buffer view is too small in " + label);
  }
}

inline std::size_t checked_mpi_count(int count, const std::string& label) {
  if (count < 0) throw std::runtime_error("negative MPI count in " + label);
  return static_cast<std::size_t>(count);
}

inline std::size_t checked_mpi_count_span(std::span<const int> counts,
                                          std::span<const int> displs,
                                          const std::string& label) {
  std::size_t required = 0;
  for (std::size_t i = 0; i < counts.size(); ++i) {
    if (counts[i] < 0 || displs[i] < 0) {
      throw std::runtime_error("negative MPI count/displacement in " + label);
    }
    required = std::max(required, static_cast<std::size_t>(displs[i]) + static_cast<std::size_t>(counts[i]));
  }
  return required;
}

template <class SendView, class RecvView>
void alltoall_complex_device(const SendView& send,
                             const RecvView& recv,
                             int count,
                             MpiComm comm,
                             const std::string& label) {
  const auto peers = static_cast<std::size_t>(mpi_size(comm));
  const auto total = checked_mpi_count(count, label) * peers;
  require_contiguous_complex_view(send, total, label + "_send");
  require_contiguous_complex_view(recv, total, label + "_recv");
  alltoall_complex_device(static_cast<const Complex*>(send.data()),
                          static_cast<Complex*>(recv.data()),
                          count,
                          comm,
                          label);
}

template <class SendView, class RecvView>
void alltoallv_complex_device(const SendView& send,
                              std::span<const int> send_counts,
                              std::span<const int> send_displs,
                              const RecvView& recv,
                              std::span<const int> recv_counts,
                              std::span<const int> recv_displs,
                              MpiComm comm,
                              const std::string& label) {
  require_contiguous_complex_view(send, checked_mpi_count_span(send_counts, send_displs, label + "_send"), label + "_send");
  require_contiguous_complex_view(recv, checked_mpi_count_span(recv_counts, recv_displs, label + "_recv"), label + "_recv");
  alltoallv_complex_device(static_cast<const Complex*>(send.data()),
                           send_counts,
                           send_displs,
                           static_cast<Complex*>(recv.data()),
                           recv_counts,
                           recv_displs,
                           comm,
                           label);
}

template <class SendView, class RecvView>
void allgather_complex_device(const SendView& send,
                              int send_count,
                              const RecvView& recv,
                              MpiComm comm,
                              const std::string& label) {
  const auto send_values = checked_mpi_count(send_count, label);
  const auto recv_values = send_values * static_cast<std::size_t>(mpi_size(comm));
  require_contiguous_complex_view(send, send_values, label + "_send");
  require_contiguous_complex_view(recv, recv_values, label + "_recv");
  allgather_complex_device(static_cast<const Complex*>(send.data()),
                           send_count,
                           static_cast<Complex*>(recv.data()),
                           comm,
                           label);
}

template <class SendView, class RecvView>
void allreduce_sum_complex_device(const SendView& send,
                                  const RecvView& recv,
                                  int count,
                                  MpiComm comm,
                                  const std::string& label) {
  require_contiguous_complex_view(send, checked_mpi_count(count, label + "_send"), label + "_send");
  require_contiguous_complex_view(recv, checked_mpi_count(count, label + "_recv"), label + "_recv");
  allreduce_sum_complex_device(static_cast<const Complex*>(send.data()),
                               static_cast<Complex*>(recv.data()),
                               count,
                               comm,
                               label);
}

template <class SendView, class RecvView>
void gatherv_complex_device(const SendView& send,
                            int send_count,
                            const RecvView& recv,
                            std::span<const int> recv_counts,
                            std::span<const int> recv_displs,
                            int root,
                            MpiComm comm,
                            const std::string& label) {
  require_contiguous_complex_view(send, checked_mpi_count(send_count, label + "_send"), label + "_send");
  require_contiguous_complex_view(recv, checked_mpi_count_span(recv_counts, recv_displs, label + "_recv"), label + "_recv");
  gatherv_complex_device(static_cast<const Complex*>(send.data()),
                         send_count,
                         static_cast<Complex*>(recv.data()),
                         recv_counts,
                         recv_displs,
                         root,
                         comm,
                         label);
}

template <class View>
void bcast_complex_device(const View& values, int count, int root, MpiComm comm, const std::string& label) {
  require_contiguous_complex_view(values, checked_mpi_count(count, label), label);
  bcast_complex_device(static_cast<Complex*>(values.data()), count, root, comm, label);
}

template <class SendView, class RecvView>
void sendrecv_complex_device(const SendView& send,
                             int send_count,
                             int dest,
                             int send_tag,
                             const RecvView& recv,
                             int recv_count,
                             int source,
                             int recv_tag,
                             MpiComm comm,
                             const std::string& label) {
  require_contiguous_complex_view(send, checked_mpi_count(send_count, label + "_send"), label + "_send");
  require_contiguous_complex_view(recv, checked_mpi_count(recv_count, label + "_recv"), label + "_recv");
  sendrecv_complex_device(static_cast<const Complex*>(send.data()),
                          send_count,
                          dest,
                          send_tag,
                          static_cast<Complex*>(recv.data()),
                          recv_count,
                          source,
                          recv_tag,
                          comm,
                          label);
}

} // namespace channel
