#include "channel/mpi.hpp"

#include "channel/runtime.hpp"

#include <stdexcept>

namespace channel {

int mpi_rank(MpiComm comm) {
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  return rank;
}

int mpi_size(MpiComm comm) {
  int size = 1;
  MPI_Comm_size(comm, &size);
  return size;
}

void mpi_barrier(MpiComm comm) {
  MPI_Barrier(comm);
}

MpiComm mpi_split(MpiComm parent, int color, int key) {
  MpiComm child = MPI_COMM_NULL;
  const int err = MPI_Comm_split(parent, color, key, &child);
  if (err != MPI_SUCCESS) throw std::runtime_error("MPI_Comm_split failed");
  return child;
}

void mpi_free(MpiComm& comm) {
  if (comm != MPI_COMM_NULL && comm != MPI_COMM_WORLD) {
    MPI_Comm_free(&comm);
  }
  comm = MPI_COMM_NULL;
}

void alltoall_complex_device(const Complex* send, Complex* recv, int count, MpiComm comm, const std::string& label) {
  fence((label + "_before_alltoall").c_str());
  const int err = MPI_Alltoall(send, count, MPI_C_DOUBLE_COMPLEX,
                               recv, count, MPI_C_DOUBLE_COMPLEX, comm);
  if (err != MPI_SUCCESS) throw std::runtime_error("MPI_Alltoall failed in " + label);
  fence((label + "_after_alltoall").c_str());
}

void alltoallv_complex_device(const Complex* send,
                               std::span<const int> send_counts,
                               std::span<const int> send_displs,
                               Complex* recv,
                               std::span<const int> recv_counts,
                               std::span<const int> recv_displs,
                               MpiComm comm,
                               const std::string& label) {
  const int peers = mpi_size(comm);
  if (send_counts.size() != static_cast<std::size_t>(peers) ||
      send_displs.size() != static_cast<std::size_t>(peers) ||
      recv_counts.size() != static_cast<std::size_t>(peers) ||
      recv_displs.size() != static_cast<std::size_t>(peers)) {
    throw std::runtime_error("MPI_Alltoallv count/displacement size mismatch in " + label);
  }
  fence((label + "_before_alltoallv").c_str());
  const int err = MPI_Alltoallv(send, send_counts.data(), send_displs.data(), MPI_C_DOUBLE_COMPLEX,
                                recv, recv_counts.data(), recv_displs.data(), MPI_C_DOUBLE_COMPLEX, comm);
  if (err != MPI_SUCCESS) throw std::runtime_error("MPI_Alltoallv failed in " + label);
  fence((label + "_after_alltoallv").c_str());
}

void allgather_complex_device(const Complex* send, int send_count, Complex* recv, MpiComm comm, const std::string& label) {
  fence((label + "_before_allgather").c_str());
  const int err = MPI_Allgather(send, send_count, MPI_C_DOUBLE_COMPLEX,
                                recv, send_count, MPI_C_DOUBLE_COMPLEX, comm);
  if (err != MPI_SUCCESS) throw std::runtime_error("MPI_Allgather failed in " + label);
  fence((label + "_after_allgather").c_str());
}

void allreduce_sum_complex_device(const Complex* send,
                                  Complex* recv,
                                  int count,
                                  MpiComm comm,
                                  const std::string& label) {
  fence((label + "_before_allreduce").c_str());
  const int err = MPI_Allreduce(send, recv, count, MPI_C_DOUBLE_COMPLEX, MPI_SUM, comm);
  if (err != MPI_SUCCESS) throw std::runtime_error("MPI_Allreduce failed in " + label);
  fence((label + "_after_allreduce").c_str());
}

void gatherv_complex_device(const Complex* send,
                            int send_count,
                            Complex* recv,
                            std::span<const int> recv_counts,
                            std::span<const int> recv_displs,
                            int root,
                            MpiComm comm,
                            const std::string& label) {
  const int peers = mpi_size(comm);
  if (recv_counts.size() != static_cast<std::size_t>(peers) ||
      recv_displs.size() != static_cast<std::size_t>(peers)) {
    throw std::runtime_error("MPI_Gatherv count/displacement size mismatch in " + label);
  }
  fence((label + "_before_gatherv").c_str());
  const int err = MPI_Gatherv(send, send_count, MPI_C_DOUBLE_COMPLEX,
                              recv, recv_counts.data(), recv_displs.data(), MPI_C_DOUBLE_COMPLEX,
                              root, comm);
  if (err != MPI_SUCCESS) throw std::runtime_error("MPI_Gatherv failed in " + label);
  fence((label + "_after_gatherv").c_str());
}

void bcast_complex_device(Complex* values, int count, int root, MpiComm comm, const std::string& label) {
  fence((label + "_before_bcast").c_str());
  const int err = MPI_Bcast(values, count, MPI_C_DOUBLE_COMPLEX, root, comm);
  if (err != MPI_SUCCESS) throw std::runtime_error("MPI_Bcast failed in " + label);
  fence((label + "_after_bcast").c_str());
}

void sendrecv_complex_device(const Complex* send,
                             int send_count,
                             int dest,
                             int send_tag,
                             Complex* recv,
                             int recv_count,
                             int source,
                             int recv_tag,
                             MpiComm comm,
                             const std::string& label) {
  fence((label + "_before_sendrecv").c_str());
  const int err = MPI_Sendrecv(send,
                               send_count,
                               MPI_C_DOUBLE_COMPLEX,
                               dest,
                               send_tag,
                               recv,
                               recv_count,
                               MPI_C_DOUBLE_COMPLEX,
                               source,
                               recv_tag,
                               comm,
                               MPI_STATUS_IGNORE);
  if (err != MPI_SUCCESS) throw std::runtime_error("MPI_Sendrecv failed in " + label);
  fence((label + "_after_sendrecv").c_str());
}

} // namespace channel
