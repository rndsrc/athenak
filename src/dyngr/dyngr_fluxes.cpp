//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file dyngr_fluxes.cpp
//  \brief Calculate 3D fluxes for hydro

#include <iostream>
//#include <stdio.h>

#include "athena.hpp"
#include "athena_tensor.hpp"
#include "mesh/mesh.hpp"
#include "dyngr.hpp"
#include "dyngr_util.hpp"
#include "adm/adm.hpp"
#include "eos/eos.hpp"
#include "diffusion/viscosity.hpp"
#include "diffusion/conduction.hpp"
#include "reconstruct/dc.hpp"
#include "reconstruct/plm.hpp"
#include "reconstruct/ppm.hpp"
#include "reconstruct/wenoz.hpp"
// include inlined Riemann solvers (double yuck...)
#include "dyngr/rsolvers/llf_dyngrhyd.cpp" // NOLINT(build/include)
// include PrimitiveSolver stuff
#include "eos/primitive-solver/idealgas.hpp"
#include "eos/primitive-solver/reset_floor.hpp"

namespace dyngr {
//----------------------------------------------------------------------------------------
//! \fn  void Hydro::CalcFluxes
//! \brief Calls reconstruction and Riemann solver functions to compute hydro fluxes
//! Note this function is templated over RS for better performance on GPUs.

template<class EOSPolicy, class ErrorPolicy> template <DynGR_RSolver rsolver_method_>
TaskStatus DynGRPS<EOSPolicy, ErrorPolicy>::CalcFluxes(Driver *pdriver, int stage) {
  RegionIndcs indcs_ = pmy_pack->pmesh->mb_indcs;
  int is = indcs_.is, ie = indcs_.ie;
  int js = indcs_.js, je = indcs_.je;
  int ks = indcs_.ks, ke = indcs_.ke;
  int ncells1 = indcs_.nx1 + 2*(indcs_.ng);

  int nhyd = pmy_pack->phydro->nhydro;
  int nvars = pmy_pack->phydro->nhydro + pmy_pack->phydro->nscalars;
  int nmb1 = pmy_pack->nmb_thispack - 1;
  const auto recon_method_ = pmy_pack->phydro->recon_method;
  auto size_ = pmy_pack->pmb->mb_size;
  auto coord_ = pmy_pack->pcoord->coord_data;
  auto w0_ = pmy_pack->phydro->w0;
  auto &adm = pmy_pack->padm->adm;
  auto &eos_ = pmy_pack->phydro->peos->eos_data;
  auto &dyn_eos_ = eos;
  const auto rsolver_ = rsolver_method_;
  bool extrema = false;
  if (recon_method_ == ReconstructionMethod::ppmx) {
    extrema = true;
  }

  //--------------------------------------------------------------------------------------
  // i-direction

  size_t scr_size = ScrArray2D<Real>::shmem_size(nvars, ncells1) * 2 + ScrArray1D<Real>::shmem_size(ncells1)
                       + ScrArray2D<Real>::shmem_size(3, ncells1) + ScrArray2D<Real>::shmem_size(6, ncells1);
  int scr_level = scratch_level;
  auto flx1_ = pmy_pack->phydro->uflx.x1f;

  par_for_outer("dyngrflux_x1",DevExeSpace(), scr_size, scr_level, 0, nmb1, ks, ke, js, je,
  KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k, const int j) {
    ScrArray2D<Real> wl(member.team_scratch(scr_level), nvars, ncells1);
    ScrArray2D<Real> wr(member.team_scratch(scr_level), nvars, ncells1);

    // scratch memory for metric at faces
    AthenaScratchTensor<Real, TensorSymm::SYM2, 3, 2> gface1_dd;
    AthenaScratchTensor<Real, TensorSymm::NONE, 3, 1> betaface1_u;
    AthenaScratchTensor<Real, TensorSymm::NONE, 3, 0> alphaface1;

    gface1_dd.NewAthenaScratchTensor(member, scr_level, ncells1);
    betaface1_u.NewAthenaScratchTensor(member, scr_level, ncells1);
    alphaface1.NewAthenaScratchTensor(member, scr_level, ncells1);
    // Reconstruct qR[i] and qL[i+1]
    switch (recon_method_) {
      case ReconstructionMethod::dc:
        DonorCellX1(member, m, k, j, is-1, ie+1, w0_, wl, wr);
        break;
      case ReconstructionMethod::plm:
        PiecewiseLinearX1(member, m, k, j, is-1, ie+1, w0_, wl, wr);
        break;
      // JF: These higher-order reconstruction methods all need EOS_Data to calculate a floor.
      case ReconstructionMethod::ppm4:
      case ReconstructionMethod::ppmx:
        PiecewiseParabolicX1(member,eos_,extrema,false, m, k, j, is-1, ie+1, w0_, wl, wr);
        break;
      case ReconstructionMethod::wenoz:
        WENOZX1(member, eos_, false, m, k, j, is-1, ie+1, w0_, wl, wr);
        break;
      default:
        break;
    }
    // Sync all threads in the team so that scratch memory is consistent
    member.team_barrier();

    // Calculate metric at faces
    Face1Metric(member, m, k, j, is-1, ie+1, adm.g_dd, adm.beta_u, adm.alpha, gface1_dd, betaface1_u, alphaface1);

    // TODO do I need a member team barrier here?
    member.team_barrier();

    // compute fluxes over [is,ie+1]
    auto &dyn_eos = dyn_eos_;
    auto &indcs = indcs_;
    auto &size = size_;
    auto &coord = coord_;
    auto &flx1 = flx1_;
    auto &nhyd_ = nhyd;
    auto nscal_ = nvars - nhyd;
    const auto rsolver = rsolver_;
    int il = is; int iu = ie+1;
    if constexpr (rsolver == DynGR_RSolver::llf_dyngr) {
      LLF_DYNGR(member, dyn_eos, indcs, size, coord, m, k, j, il, iu, IVX, 
                wl, wr, nhyd_, nscal_, gface1_dd, betaface1_u, alphaface1, flx1);
    }
    //} else { other Riemann solvers here
    member.team_barrier();

    // Calculate fluxes of scalars (if any)
    if (nvars > nhyd) {
      for (int n=nhyd; n<nvars; ++n) {
        par_for_inner(member, is, ie+1, [&](const int i) {
          if (flx1(m,IDN,k,j,i) >= 0.0) {
            flx1(m,n,k,j,i) = flx1(m,IDN,k,j,i)*wl(n,i);
          } else {
            flx1(m,n,k,j,i) = flx1(m,IDN,k,j,i)*wr(n,i);
          }
        });
      }
    }
    member.team_barrier();
  });

  //--------------------------------------------------------------------------------------
  // j-direction

  if (pmy_pack->pmesh->multi_d) {
    scr_size = ScrArray2D<Real>::shmem_size(nvars, ncells1) * 3 + ScrArray1D<Real>::shmem_size(ncells1)
                         + ScrArray2D<Real>::shmem_size(3, ncells1) + ScrArray2D<Real>::shmem_size(6, ncells1);
    auto flx2_ = pmy_pack->phydro->uflx.x2f;

    par_for_outer("dyngrflux_x2",DevExeSpace(), scr_size, scr_level, 0, nmb1, ks, ke,
    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int k) {
      ScrArray2D<Real> scr1(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr2(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr3(member.team_scratch(scr_level), nvars, ncells1);

      // scratch memory for metric at faces
      AthenaScratchTensor<Real, TensorSymm::SYM2, 3, 2> gface2_dd;
      AthenaScratchTensor<Real, TensorSymm::NONE, 3, 1> betaface2_u;
      AthenaScratchTensor<Real, TensorSymm::NONE, 3, 0> alphaface2;

      gface2_dd.NewAthenaScratchTensor(member, scr_level, ncells1);
      betaface2_u.NewAthenaScratchTensor(member, scr_level, ncells1);
      alphaface2.NewAthenaScratchTensor(member, scr_level, ncells1);
      for (int j=js-1; j<=je+1; ++j) {
        // Permute scratch arrays.
        auto wl     = scr1;
        auto wl_jp1 = scr2;
        auto wr     = scr3;
        if ((j%2) == 0) {
          wl     = scr2;
          wl_jp1 = scr1;
        }

        // Reconstruct qR[j] and qL[j+1]
        switch (recon_method_) {
          case ReconstructionMethod::dc:
            DonorCellX2(member, m, k, j, is, ie, w0_, wl_jp1, wr);
            break;
          case ReconstructionMethod::plm:
            PiecewiseLinearX2(member, m, k, j, is, ie, w0_, wl_jp1, wr);
            break;
          // JF: These higher-order reconstruction methods all need EOS_Data to calculate a floor.
          case ReconstructionMethod::ppm4:
          case ReconstructionMethod::ppmx:
            PiecewiseParabolicX2(member,eos_,extrema,false, m, k, j, is, ie, w0_, wl_jp1, wr);
            break;
          case ReconstructionMethod::wenoz:
            WENOZX2(member, eos_, false, m, k, j, is-1, ie+1, w0_, wl_jp1, wr);
            break;
          default:
            break;
        }
        // Sync all threads in the team so that scratch memory is consistent
        member.team_barrier();

        // Calculate metric at faces
        Face2Metric(member, m, k, j, is, ie, adm.g_dd, adm.beta_u, adm.alpha, gface2_dd, betaface2_u, alphaface2);

        // TODO do I need a member team barrier here?
        member.team_barrier();

        // compute fluxes over [is,ie+1]
        auto &dyn_eos = dyn_eos_;
        auto &indcs = indcs_;
        auto &size = size_;
        auto &coord = coord_;
        auto &flx2 = flx2_;
        auto &nhyd_ = nhyd;
        auto nscal_ = nvars - nhyd;
        const auto rsolver = rsolver_;
        int il = is; int iu = ie;
        if (j>(js-1)) {
          if constexpr (rsolver == DynGR_RSolver::llf_dyngr) {
            LLF_DYNGR(member, dyn_eos, indcs, size, coord, m, k, j, il, iu, IVY, 
                      wl, wr, nhyd_, nscal_, gface2_dd, betaface2_u, alphaface2, flx2);
          }
          //} else { other Riemann solvers here
        }
        member.team_barrier();

        // Calculate fluxes of scalars (if any)
        if (nvars > nhyd) {
          for (int n=nhyd; n<nvars; ++n) {
            par_for_inner(member, is, ie, [&](const int i) {
              if (flx2(m,IDN,k,j,i) >= 0.0) {
                flx2(m,n,k,j,i) = flx2(m,IDN,k,j,i)*wl(n,i);
              } else {
                flx2(m,n,k,j,i) = flx2(m,IDN,k,j,i)*wr(n,i);
              }
            });
          }
        }
      } // end of loop over j
      member.team_barrier();
    });
  }

  //--------------------------------------------------------------------------------------
  // k-direction. Note order of k,j loops switched

  if (pmy_pack->pmesh->three_d) {
    scr_size = ScrArray2D<Real>::shmem_size(nvars, ncells1) * 3 + ScrArray1D<Real>::shmem_size(ncells1)
                         + ScrArray2D<Real>::shmem_size(3, ncells1) + ScrArray2D<Real>::shmem_size(6, ncells1);
    auto flx3_ = pmy_pack->phydro->uflx.x3f;

    par_for_outer("dyngrflux_x3",DevExeSpace(), scr_size, scr_level, 0, nmb1, js, je,
    KOKKOS_LAMBDA(TeamMember_t member, const int m, const int j) {
      ScrArray2D<Real> scr1(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr2(member.team_scratch(scr_level), nvars, ncells1);
      ScrArray2D<Real> scr3(member.team_scratch(scr_level), nvars, ncells1);

      // scratch memory for metric at faces
      AthenaScratchTensor<Real, TensorSymm::SYM2, 3, 2> gface3_dd;
      AthenaScratchTensor<Real, TensorSymm::NONE, 3, 1> betaface3_u;
      AthenaScratchTensor<Real, TensorSymm::NONE, 3, 0> alphaface3;

      gface3_dd.NewAthenaScratchTensor(member, scr_level, ncells1);
      betaface3_u.NewAthenaScratchTensor(member, scr_level, ncells1);
      alphaface3.NewAthenaScratchTensor(member, scr_level, ncells1);
      for (int k=ks-1; k<=ke+1; ++k) {
        // Permute scratch arrays.
        auto wl     = scr1;
        auto wl_kp1 = scr2;
        auto wr     = scr3;
        if ((k%2) == 0) {
          wl     = scr2;
          wl_kp1 = scr1;
        }

        // Reconstruct qR[j] and qL[j+1]
        switch (recon_method_) {
          case ReconstructionMethod::dc:
            DonorCellX3(member, m, k, j, is, ie, w0_, wl_kp1, wr);
            break;
          case ReconstructionMethod::plm:
            PiecewiseLinearX3(member, m, k, j, is, ie, w0_, wl_kp1, wr);
            break;
          // JF: These higher-order reconstruction methods all need EOS_Data to calculate a floor.
          case ReconstructionMethod::ppm4:
          case ReconstructionMethod::ppmx:
            PiecewiseParabolicX3(member,eos_,extrema,false, m, k, j, is, ie, w0_, wl_kp1, wr);
            break;
          case ReconstructionMethod::wenoz:
            WENOZX3(member, eos_, false, m, k, j, is-1, ie+1, w0_, wl_kp1, wr);
            break;
          default:
            break;
        }
        // Sync all threads in the team so that scratch memory is consistent
        member.team_barrier();

        // Calculate metric at faces
        Face3Metric(member, m, k, j, is, ie, adm.g_dd, adm.beta_u, adm.alpha, gface3_dd, betaface3_u, alphaface3);

        // TODO do I need a member team barrier here?
        member.team_barrier();

        // compute fluxes over [ks,ke+1]
        auto &dyn_eos = dyn_eos_;
        auto &indcs = indcs_;
        auto &size = size_;
        auto &coord = coord_;
        auto &flx3 = flx3_;
        auto &nhyd_ = nhyd;
        auto nscal_ = nvars - nhyd;
        const auto rsolver = rsolver_;
        int il = is; int iu = ie;
        if (k>(ks-1)) {
          if constexpr (rsolver == DynGR_RSolver::llf_dyngr) {
            LLF_DYNGR(member, dyn_eos, indcs, size, coord, m, k, j, il, iu, IVZ,
                      wl, wr, nhyd_, nscal_, gface3_dd, betaface3_u, alphaface3, flx3);
          }
          //} else { other Riemann solvers here
        }
        member.team_barrier();

        // Calculate fluxes of scalars (if any)
        if (nvars > nhyd) {
          for (int n=nhyd; n<nvars; ++n) {
            par_for_inner(member, is, ie, [&](const int i) {
              if (flx3(m,IDN,k,j,i) >= 0.0) {
                flx3(m,n,k,j,i) = flx3(m,IDN,k,j,i)*wl(n,i);
              } else {
                flx3(m,n,k,j,i) = flx3(m,IDN,k,j,i)*wr(n,i);
              }
            });
          }
        }
      } // end of loop over j
      member.team_barrier();
    });
    // TODO: handle excision masks
  }

  // excision masks
  /*if (coord_.bh_excise) {
    auto &fc_mask_ = pmy_pack->pcoord->fc_mask;

    auto &fcorr_x1 = pmy_pack->phydro->uflx.x1f;
    auto &fcorr_x2 = pmy_pack->phydro->uflx.x2f;
    auto &fcorr_x3 = pmy_pack->phydro->uflx.x3f;
    par_for("excise_flux", DevExeSpace(), 0, nmb1, ks, ke+1, js, je+1, is, ie+1,
    KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
      if (j<(je+1) && k<(ke+1)) {
        if (fc_mask_.x1f(m,k,j,i)) {
          Real wim1[NPRIM];
          ExtractPrimitives(w0_, wim1, dyn_eos_, nhyd, nvars - nhyd, m, k, j, i-1);

          Real wi[NPRIM];
          ExtractPrimitives(w0_, wi, dyn_eos_, nhyd, nvars - nhyd, m, k, j, i);

          Real g3d[NSPMETRIC], beta_u[3], alpha;
          Face1Metric(m, k, j, i, adm.g_dd, adm.beta_u, adm.alpha, g3d, beta_u, alpha);

          Real flux[NCONS];
          SingleStateLLF_DYNGR(dyn_eos_, wim1, wi, IVX, g3d, beta_u, alpha, flux);

          fcorr_x1(m, IDN, k, j, i) = flux[CDN];
          fcorr_x1(m, IM1, k, j, i) = flux[CSX];
          fcorr_x1(m, IM2, k, j, i) = flux[CSY];
          fcorr_x1(m, IM3, k, j, i) = flux[CSZ];
          fcorr_x1(m, IEN, k, j, i) = flux[CTA];
        }
      }

      if (i<(ie+1) && k<(ke+1)) {
        if (fc_mask_.x2f(m,k,j,i)) {
          Real wjm1[NPRIM];
          ExtractPrimitives(w0_, wjm1, dyn_eos_, nhyd, nvars - nhyd, m, k, j-1, i);

          Real wj[NPRIM];
          ExtractPrimitives(w0_, wj, dyn_eos_, nhyd, nvars - nhyd, m, k, j, i);

          Real g3d[NSPMETRIC], beta_u[3], alpha;
          Face2Metric(m, k, j, i, adm.g_dd, adm.beta_u, adm.alpha, g3d, beta_u, alpha);

          Real flux[NCONS];
          SingleStateLLF_DYNGR(dyn_eos_, wjm1, wj, IVY, g3d, beta_u, alpha, flux);

          fcorr_x2(m, IDN, k, j, i) = flux[CDN];
          fcorr_x2(m, IM1, k, j, i) = flux[CSX];
          fcorr_x2(m, IM2, k, j, i) = flux[CSY];
          fcorr_x2(m, IM3, k, j, i) = flux[CSZ];
          fcorr_x2(m, IEN, k, j, i) = flux[CTA];
        }
      }

      if (i<(ie+1) && j<(je+1)) {
        if (fc_mask_.x3f(m,k,j,i)) {
          Real wkm1[NPRIM];
          ExtractPrimitives(w0_, wkm1, dyn_eos_, nhyd, nvars - nhyd, m, k-1, j, i);

          Real wk[NPRIM];
          ExtractPrimitives(w0_, wk, dyn_eos_, nhyd, nvars - nhyd, m, k, j, i);
          
          Real g3d[NSPMETRIC], beta_u[3], alpha;
          Face3Metric(m, k, j, i, adm.g_dd, adm.beta_u, adm.alpha, g3d, beta_u, alpha);

          Real flux[NCONS];
          SingleStateLLF_DYNGR(dyn_eos_, wkm1, wk, IVZ, g3d, beta_u, alpha, flux);

          fcorr_x3(m, IDN, k, j, i) = flux[CDN];
          fcorr_x3(m, IM1, k, j, i) = flux[CSX];
          fcorr_x3(m, IM2, k, j, i) = flux[CSY];
          fcorr_x3(m, IM3, k, j, i) = flux[CSZ];
          fcorr_x3(m, IEN, k, j, i) = flux[CTA];
        }
      }
    });
  }*/
  
  return TaskStatus::complete;
}

// function definitions for each template parameter
// Macro for instantiating every flux function for each Riemann solver
#define INSTANTIATE_CALC_FLUXES(EOSPolicy, ErrorPolicy) \
template \
TaskStatus DynGRPS<EOSPolicy, ErrorPolicy>::CalcFluxes<DynGR_RSolver::llf_dyngr>(Driver *pdriver, int stage);

INSTANTIATE_CALC_FLUXES(Primitive::IdealGas, Primitive::ResetFloor)

} // namespace dyngr