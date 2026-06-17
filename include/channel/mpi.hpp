#pragma once

#include "channel/types.hpp"

#include <span>
#include <string>

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
void gatherv_complex_device(const Complex* send,
                            int send_count,
                            Complex* recv,
                            std::span<const int> recv_counts,
                            std::span<const int> recv_displs,
                            int root,
                            MpiComm comm,
                            const std::string& label);
void bcast_complex_device(Complex* values, int count, int root, MpiComm comm, const std::string& label);

} // namespace channel
