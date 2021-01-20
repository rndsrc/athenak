//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file shock_tube.cpp
//  \brief Problem generator for shock tube problems.
//
// Problem generator for shock tube (1-D Riemann) problems. Initializes plane-parallel
// shock along x1 (in 1D, 2D, 3D), along x2 (in 2D, 3D), and along x3 (in 3D).

#include <iostream>
#include <sstream>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "hydro/eos/hydro_eos.hpp"
#include "hydro/hydro.hpp"
#include "utils/grid_locations.hpp"
#include "pgen.hpp"

//----------------------------------------------------------------------------------------
//! \fn
//  \brief Problem Generator for the shock tube (Riemann problem) tests

void ProblemGenerator::ShockTube_(MeshBlockPack *pmbp, ParameterInput *pin)
{
  using namespace hydro;

  // parse shock direction: {1,2,3} -> {x1,x2,x3}
  int shk_dir = pin->GetInteger("problem","shock_dir");

  // parse shock location (must be inside grid)
  Real xshock = pin->GetReal("problem","xshock");
  if (shk_dir == 1 && (xshock < pmesh_->mesh_size.x1min ||
                       xshock > pmesh_->mesh_size.x1max)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "xshock=" << xshock << " lies outside x1 domain" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (shk_dir == 2 && (xshock < pmesh_->mesh_size.x2min ||
                       xshock > pmesh_->mesh_size.x2max)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "xshock=" << xshock << " lies outside x2 domain" << std::endl;
    exit(EXIT_FAILURE);
  }
  if (shk_dir == 3 && (xshock < pmesh_->mesh_size.x3min ||
                       xshock > pmesh_->mesh_size.x3max)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "xshock=" << xshock << " lies outside x3 domain" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Parse left state read from input file: dl,ul,vl,wl,[pl]
  Real wl[5];
  wl[IDN] = pin->GetReal("problem","dl");
  wl[IVX] = pin->GetReal("problem","ul");
  wl[IVY] = pin->GetReal("problem","vl");
  wl[IVZ] = pin->GetReal("problem","wl");
  wl[IPR] = pin->GetReal("problem","pl");

  // Parse right state read from input file: dr,ur,vr,wr,[pr]
  Real wr[5];
  wr[IDN] = pin->GetReal("problem","dr");
  wr[IVX] = pin->GetReal("problem","ur");
  wr[IVY] = pin->GetReal("problem","vr");
  wr[IVZ] = pin->GetReal("problem","wr");
  wr[IPR] = pin->GetReal("problem","pr");

  // Initialize the discontinuity in the Hydro variables ---------------------------------

 // capture variables for the kernel
  Real gm1 = pmbp->phydro->peos->eos_data.gamma - 1.0;
  int &is = pmbp->mb_cells.is, &ie = pmbp->mb_cells.ie;
  int &js = pmbp->mb_cells.js, &je = pmbp->mb_cells.je;
  int &ks = pmbp->mb_cells.ks, &ke = pmbp->mb_cells.ke;
  int &nx1 = pmbp->mb_cells.nx1;
  int &nx2 = pmbp->mb_cells.nx2;
  int &nx3 = pmbp->mb_cells.nx3;
  auto &u0 = pmbp->phydro->u0;
  auto size = pmbp->pmb->mbsize;

  switch(shk_dir) {

    //--- shock in 1-direction
    case 1:
      par_for("pgen_shock1", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
        KOKKOS_LAMBDA(int m,int k, int j, int i)
        {
          Real x1 = CellCenterX(i-is, nx1, size.x1min.d_view(m), size.x1max.d_view(m));
          if (x1 < xshock) {
            u0(m,IDN,k,j,i) = wl[IDN];
            u0(m,IM1,k,j,i) = wl[IVX]*wl[IDN];
            u0(m,IM2,k,j,i) = wl[IVY]*wl[IDN];
            u0(m,IM3,k,j,i) = wl[IVZ]*wl[IDN];
            u0(m,IEN,k,j,i) = wl[IPR]/gm1 +
               0.5*wl[IDN]*(SQR(wl[IVX]) + SQR(wl[IVY]) + SQR(wl[IVZ]));
          } else {
            u0(m,IDN,k,j,i) = wr[IDN];
            u0(m,IM1,k,j,i) = wr[IVX]*wr[IDN];
            u0(m,IM2,k,j,i) = wr[IVY]*wr[IDN];
            u0(m,IM3,k,j,i) = wr[IVZ]*wr[IDN];
            u0(m,IEN,k,j,i) = wr[IPR]/gm1 +
               0.5*wr[IDN]*(SQR(wr[IVX]) + SQR(wr[IVY]) + SQR(wr[IVZ]));
          }
        }
      );
      break;

    //--- shock in 2-direction
    case 2:
      par_for("pgen_shock2", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
        KOKKOS_LAMBDA(int m,int k, int j, int i)
        {
          Real x2 = CellCenterX(j-js, nx2, size.x2min.d_view(m), size.x2max.d_view(m));
          if (x2 < xshock) {
            u0(m,IDN,k,j,i) = wl[IDN];
            u0(m,IM2,k,j,i) = wl[IVX]*wl[IDN];
            u0(m,IM3,k,j,i) = wl[IVY]*wl[IDN];
            u0(m,IM1,k,j,i) = wl[IVZ]*wl[IDN];
            u0(m,IEN,k,j,i) = wl[IPR]/gm1 +
               0.5*wl[IDN]*(SQR(wl[IVX]) + SQR(wl[IVY]) + SQR(wl[IVZ]));
          } else {
            u0(m,IDN,k,j,i) = wr[IDN];
            u0(m,IM2,k,j,i) = wr[IVX]*wr[IDN];
            u0(m,IM3,k,j,i) = wr[IVY]*wr[IDN];
            u0(m,IM1,k,j,i) = wr[IVZ]*wr[IDN];
            u0(m,IEN,k,j,i) = wr[IPR]/gm1 +
               0.5*wr[IDN]*(SQR(wr[IVX]) + SQR(wr[IVY]) + SQR(wr[IVZ]));
          }
        }
      );
      break;

    //--- shock in 3-direction
    case 3:
      par_for("pgen_shock3", DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
        KOKKOS_LAMBDA(int m,int k, int j, int i)
        {
          Real x3 = CellCenterX(k-ks, nx3, size.x3min.d_view(m), size.x3max.d_view(m));
          if (x3 < xshock) {
            u0(m,IDN,k,j,i) = wl[IDN];
            u0(m,IM3,k,j,i) = wl[IVX]*wl[IDN];
            u0(m,IM1,k,j,i) = wl[IVY]*wl[IDN];
            u0(m,IM2,k,j,i) = wl[IVZ]*wl[IDN];
            u0(m,IEN,k,j,i) = wl[IPR]/gm1 +
               0.5*wl[IDN]*(SQR(wl[IVX]) + SQR(wl[IVY]) + SQR(wl[IVZ]));
          } else {
            u0(m,IDN,k,j,i) = wr[IDN];
            u0(m,IM3,k,j,i) = wr[IVX]*wr[IDN];
            u0(m,IM1,k,j,i) = wr[IVY]*wr[IDN];
            u0(m,IM2,k,j,i) = wr[IVZ]*wr[IDN];
            u0(m,IEN,k,j,i) = wr[IPR]/gm1 +
               0.5*wr[IDN]*(SQR(wr[IVX]) + SQR(wr[IVY]) + SQR(wr[IVZ]));
          }
        }
      );
      break;

    //--- Invaild input value for shk_dir
    default:
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "shock_dir=" <<shk_dir<< " must be either 1,2, or 3" << std::endl;
      exit(EXIT_FAILURE);
  }

  return;
}
