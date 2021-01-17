//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file send_buffers.cpp
//  \brief packs and sends boundary buffers between MeshBlocks

#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "bvals/bvals.hpp"
#include "hydro/hydro.hpp"

//----------------------------------------------------------------------------------------
// \!fn void SendBuffers()
// \brief Pack boundary buffers and send to neighbors.
// This routine packs ALL the buffers on ALL the faces, edges, and corners simultaneously,
// for ALL the MeshBlocks.  This reduces the number of kernel launches when there are a
// large number of MeshBlocks per MPI rank.  Buffer data are then sent (via MPI) or copied
// directly for periodic or block boundaries.
//
// Input array must be 5D Kokkos View dimensioned (nmb, nvar, nx3, nx2, nx1)
//
// TODO: with AMR, buffer indices can be different on different MeshBlocks

TaskStatus SendBuffers(AthenaArray5D<Real> &a,
  std::vector<std::vector<BoundaryBuffer>> send_buf,
  std::vector<std::vector<BoundaryBuffer>> recv_buf, std::vector<MeshBlock> mblocks)
{
  // create local references for variables in kernel
  int nnghbr = send_buf.size();
  int nmb  = mblocks.size();
  int nvar = a.extent_int(1);  // 2nd index from L of input array must be NVAR
  int &my_rank = global_variable::my_rank;

  // load buffers, using 3 levels of hierarchical parallelism
  int scr_level = 0;
  size_t scr_size = 16*sizeof(int);
  par_for_outer("SendCC", DevExeSpace(), scr_size, scr_level, 0, (nmb-1), 0, (nnghbr-1),
    KOKKOS_LAMBDA(TeamMember_t tmember, const int m, const int n)
    {
      const int il = send_buf[m][n].index(0);
      const int iu = send_buf[m][n].index(1);
      const int jl = send_buf[m][n].index(2);
      const int ju = send_buf[m][n].index(3);
      const int kl = send_buf[m][n].index(4);
      const int ku = send_buf[m][n].index(5);
      const int ni = iu - il + 1;
      const int nj = ju - jl + 1;
      const int nk = ku - kl + 1;
      const int nkj  = nk*nj;
      const int nvkj = nvar*nk*nj;
      Kokkos::parallel_for(
        Kokkos::TeamThreadRange<>(tmember, nvkj), [&](const int idx) {
          int v = idx / nkj;
          int k = (idx - v * nkj) / nj;
          int j = idx - v * nkj - k * nj;
          k += kl;
          j += jl;
  
          // copy directly into recv buffer if MeshBlocks on same rank
          if (mblocks[m].nghbr[n].rank == my_rank) {
            int nn = mblocks[m].nghbr[n].destn; // index of recv'ing boundary buffer
            // index of recv;ing MB: assumes MB IDs are stored sequentially in mblocks[]
            int mm = mblocks[m].nghbr[n].gid - mblocks[0].mb_gid;
            Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(tmember, il, iu + 1), [&](const int i) {
                recv_buf[mm][nn].data(v, i-il + ni*(j-jl + nj*(k-kl))) = a(m, v, k, j, i);
              }
            );

          // else copy directly into send buffer for MPI communication below
          } else {
            Kokkos::parallel_for(
              Kokkos::ThreadVectorRange(tmember, il, iu + 1), [&](const int i) {
                send_buf[m][n].data(v, i-il + ni*(j-jl + nj*(k-kl))) = a(m, v, k, j, i);
              }
            );
          }
        }
      );
    }
  ); // end par_for_outer

  // Send boundary buffer to neighboring MeshBlocks using MPI or Kokkos::deep_copy if
  // neighbor is on same MPI rank.
  // Note send_buf[n] --> recv_buf[n + nghbr[n].dn]

#if MPI_PARALLEL_ENABLED
  using Kokkos::ALL;
  for (int n=0; n<nnghbr; ++n) {
    if (nghbr[n].gid >= 0) {  // ID of buffer != -1, so not a physical boundary
      if (nghbr[n].rank == global_variable::my_rank) {
//        Kokkos::deep_copy(pmb->exe_space,
//          pmesh_->FindMeshBlock(nghbr[n].gid)->phydro->recv_buf[n + nghbr[n].dn].data,
//          send_buf[n].data);
//        pmesh_->FindMeshBlock(nghbr[n].gid)->phydro->recv_buf[n + nghbr[n].dn].bcomm_stat
//          = BoundaryCommStatus::received;
      } else {
        // create tag using local ID and buffer index of *receiving* MeshBlock
        int lid = nghbr_x1face[n].gid - pmesh_->gidslist[nghbr_x1face[n].rank];
        int tag = CreateMPItag(lid, (1-n), key);
        void* send_ptr = sendbuf.data();
        int ierr = MPI_Isend(send_ptr, sendbuf.size(), MPI_ATHENA_REAL,
          nghbr_x1face[n].rank, tag, MPI_COMM_WORLD, &(pbb->send_rq_x1face[n]));
      }
    }
  }
#endif

  return TaskStatus::complete;
}