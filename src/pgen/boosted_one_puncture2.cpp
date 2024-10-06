//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file z4c_one_puncture.cpp
//  \brief Problem generator for a single puncture placed at the origin of the domain

#include <algorithm>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <iostream>   // endl
#include <limits>     // numeric_limits::max()
#include <memory>
#include <string>     // c_str(), string
#include <vector>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "z4c/z4c.hpp"
#include "adm/adm.hpp"
#include "coordinates/cell_locations.hpp"

void ADMOnePunctureBoosted(MeshBlockPack *pmbp, ParameterInput *pin);

KOKKOS_INLINE_FUNCTION
AthenaScratchTensor<Real, TensorSymm::SYM2, 4, 2> 
inverse(AthenaScratchTensor<Real, TensorSymm::SYM2, 4, 2> matrix);

KOKKOS_INLINE_FUNCTION
void LorentzBoost(Real vx1, Real vx2, Real vx3, AthenaScratchTensor<Real, TensorSymm::SYM2, 4, 2> &lambda);

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::UserProblem_()
//! \brief Problem Generator for single puncture
void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  auto &indcs = pmy_mesh_->mb_indcs;

  if (pmbp->pz4c == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "One Puncture test can only be run in Z4c, but no <z4c> block "
              << "in input file" << std::endl;
    exit(EXIT_FAILURE);
  }

  ADMOnePunctureBoosted(pmbp, pin);
  // pmbp->pz4c->GaugePreCollapsedLapse(pmbp, pin);
  switch (indcs.ng) {
    case 2: pmbp->pz4c->ADMToZ4c<2>(pmbp, pin);
            break;
    case 3: pmbp->pz4c->ADMToZ4c<3>(pmbp, pin);
            break;
    case 4: pmbp->pz4c->ADMToZ4c<4>(pmbp, pin);
            break;
  }
  pmbp->pz4c->Z4cToADM(pmbp);
  switch (indcs.ng) {
    case 2: pmbp->pz4c->ADMConstraints<2>(pmbp);
            break;
    case 3: pmbp->pz4c->ADMConstraints<3>(pmbp);
            break;
    case 4: pmbp->pz4c->ADMConstraints<4>(pmbp);
            break;
  }

  std::cout<<"OnePuncture initialized."<<std::endl;


  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ADMOnePuncture(MeshBlockPack *pmbp, ParameterInput *pin)
//! \brief Initialize ADM vars to single boosted puncture (no spin), based on 1909.02997

void ADMOnePunctureBoosted(MeshBlockPack *pmbp, ParameterInput *pin) {
  // capture variables for the kernel
  auto &indcs = pmbp->pmesh->mb_indcs;
  auto &size = pmbp->pmb->mb_size;
  int &is = indcs.is; int &ie = indcs.ie;
  int &js = indcs.js; int &je = indcs.je;
  int &ks = indcs.ks; int &ke = indcs.ke;
  // For GLOOPS
  int isg = is-indcs.ng; int ieg = ie+indcs.ng;
  int jsg = js-indcs.ng; int jeg = je+indcs.ng;
  int ksg = ks-indcs.ng; int keg = ke+indcs.ng;
  int nmb = pmbp->nmb_thispack;
  Real ADM_mass = pin->GetOrAddReal("problem", "punc_ADM_mass", 1.);
  Real center_x1 = pin->GetOrAddReal("problem", "punc_center_x1", 0.);
  Real center_x2 = pin->GetOrAddReal("problem", "punc_center_x2", 0.);
  Real center_x3 = pin->GetOrAddReal("problem", "punc_center_x3", 0.);
  Real vx1 = pin->GetOrAddReal("problem", "punc_velocity_x1", 0.);
  Real vx2 = pin->GetOrAddReal("problem", "punc_velocity_x2", 0.);
  Real vx3 = pin->GetOrAddReal("problem", "punc_velocity_x3", 0.);

  adm::ADM::ADM_vars &adm = pmbp->padm->adm;

  par_for("pgen one puncture",
  DevExeSpace(),0,nmb-1,ksg,keg,jsg,jeg,isg,ieg,
  KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {
    Real &x1min = size.d_view(m).x1min;
    Real &x1max = size.d_view(m).x1max;
    int nx1 = indcs.nx1;
    Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

    Real &x2min = size.d_view(m).x2min;
    Real &x2max = size.d_view(m).x2max;
    int nx2 = indcs.nx2;
    Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

    Real &x3min = size.d_view(m).x3min;
    Real &x3max = size.d_view(m).x3max;
    int nx3 = indcs.nx3;
    Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

    x1v -= center_x1;
    x2v -= center_x2;
    x3v -= center_x3;

    // construct lorentz boost
    AthenaScratchTensor<Real, TensorSymm::SYM2, 4, 2> lambda;
    LorentzBoost(vx1, vx2, vx3, lambda);

    // inverse Lorentz boost
    AthenaScratchTensor<Real, TensorSymm::SYM2, 4, 2> lambda_inv;
    LorentzBoost(-vx1, -vx2, -vx3, lambda_inv);

    // coordinate in the boosted frame (p or primed)
    Real xp[4];
    Real xinit[4] = {0,x1v,x2v,x3v};

    for(int a = 0; a < 4; ++a) {
      xp[a] = 0;
      for (int b = 0; b < 4; ++b) {
        xp[a] += lambda(a,b)*xinit[b];
      }
    }
    // radial coordinate in boosted frame
    Real r = std::sqrt(std::pow(xp[1],2) + std::pow(xp[2],2) + std::pow(xp[3],2));

    // metric in boosted frame
    AthenaScratchTensor<Real, TensorSymm::SYM2, 4, 2> g_p_dd;
    g_p_dd.ZeroClear();
    // conformal factor and omega (whatever it's called)
    Real psi = 1.0 + 0.25*ADM_mass/r;
    Real omega = 1.0 - 0.25*ADM_mass/r;

    g_p_dd(0,0) = - std::pow(omega,2) / std::pow(psi,2);
    for(int a = 1; a < 4; ++a) {
      g_p_dd(a,a) = std::pow(psi,4);
    }

    // metric partial derivative in boosted frame
    AthenaScratchTensor<Real, TensorSymm::SYM2, 4, 3> dg_p_ddd;
    dg_p_ddd.ZeroClear();
    // spatial partial of psi^4
    AthenaScratchTensor<Real, TensorSymm::NONE, 4, 1> dpsi4;
    for(int a = 1; a < 4; ++a) {
      dpsi4(a) = -ADM_mass*xp[a]*std::pow(psi,3)/std::pow(r,1.5);
    }

    // put together partial derivatives
    for(int a = 1; a < 4; ++a) {
      dg_p_ddd(a,0,0) = 8*xp[a]*SQR(ADM_mass-4*r)/(r*std::pow(ADM_mass+4*r,3))
              + 8*xp[a]*(ADM_mass-4*r)/(r*SQR(ADM_mass+4*r));
      for(int b = 1; b < 4; ++b) {
        dg_p_ddd(a,b,b) = dpsi4(a);
      }
    }

    // inverse metric partial derivative in boosted frame
    AthenaScratchTensor<Real, TensorSymm::SYM2, 4, 3> dg_p_duu;
    dg_p_duu.ZeroClear();
    // since the metric is diagonal, each component is simply g^aa=1/g_aa
    // the partials are just then g^bb,a = -1/g_bb^2 g_bb,a
    // put together partial derivatives
    for(int a = 1; a < 4; ++a) {
      for(int b = 0; b < 4; ++b) {
        dg_p_duu(a,b,b) = -1/SQR(g_p_dd(b,b))*dg_p_ddd(a,b,b);
      }
    }

    // Construct the boosted metric
    AthenaScratchTensor<Real, TensorSymm::SYM2, 4, 2> g_boosted_dd;
    g_boosted_dd.ZeroClear();
    for(int a = 0; a < 4; ++a)
    for(int b = a; b < 4; ++b)
    for(int c = 0; c < 4; ++c)
    for(int d = 0; d < 4; ++d) {
      g_boosted_dd(a,b) += lambda(c,a)*lambda(d,b)*g_p_dd(c,d);
    }

    // Construct the boosted metric derivative
    AthenaScratchTensor<Real, TensorSymm::SYM2, 4, 3> dg_boosted_ddd;
    dg_boosted_ddd.ZeroClear();

    for(int a = 0; a < 4; ++a)
    for(int b = a; b < 4; ++b)
    for(int c = 0; c < 4; ++c)
    for(int d = 0; d < 4; ++d)
    for(int e = 0; e < 4; ++e)
    for(int f = 0; f < 4; ++f) {
      dg_boosted_ddd(e,a,b) += lambda(f,e)*lambda(c,a)*lambda(d,b)*dg_p_ddd(f,c,d);
    }

    // Construct the boosted inverse metric derivative
    AthenaScratchTensor<Real, TensorSymm::SYM2, 4, 3> dg_boosted_duu;
    dg_boosted_duu.ZeroClear();

    for(int a = 0; a < 4; ++a)
    for(int b = a; b < 4; ++b)
    for(int c = 0; c < 4; ++c)
    for(int d = 0; d < 4; ++d)
    for(int e = 0; e < 4; ++e)
    for(int f = 0; f < 4; ++f) {
      dg_boosted_duu(e,a,b) += lambda(f,e)*lambda_inv(c,a)*lambda_inv(d,b)*dg_p_duu(f,c,d);
    }

    // Lastly calculate the ADM variable
    // First the inverted spacetime metric
    AthenaScratchTensor<Real, TensorSymm::SYM2, 4, 2> g_boosted_uu = inverse(g_boosted_dd);

    // Gauge variables
    adm.alpha(m,k,j,i) = 1/sqrt(-g_boosted_uu(0,0));
    for(int a = 0; a < 3; ++a) {
      adm.beta_u(m,a,k,j,i) = g_boosted_uu(0,a+1)*SQR(adm.alpha(m,k,j,i));
    }

    // adm metric
    for(int a = 0; a < 3; ++a)
    for(int b = a; b < 3; ++b) {
      adm.g_dd(m,a,b,k,j,i) = g_boosted_dd(a+1,b+1);
    }

    // partials for the shift vector
    AthenaScratchTensor<Real, TensorSymm::NONE, 3, 2> dbeta_du;
    for(int a = 0; a < 3; ++a)
    for(int b = 0; b < 3; ++b) {
      dbeta_du(a,b) = SQR(adm.alpha(m,k,j,i))*(dg_boosted_duu(a+1,0,b+1)
                      +dg_boosted_duu(a+1,0,0)*adm.beta_u(m,b,k,j,i));
    }

    // extrinsic curvature
    for(int a = 0; a < 3; ++a)
    for(int b = a; b < 3; ++b) {
      adm.vK_dd(m,a,b,k,j,i) = dg_boosted_ddd(0,a,b);
      for(int c = 0; c < 3; ++c) {
        adm.vK_dd(m,a,b,k,j,i) += - adm.beta_u(m,c,k,j,i)*dg_boosted_ddd(c,a,b)
                                  - g_boosted_dd(c,b)*dbeta_du(a,c)
                                  - g_boosted_dd(a,c)*dbeta_du(b,c);
      }
      adm.vK_dd(m,a,b,k,j,i) /= - 2*adm.alpha(m,k,j,i);
    }
  });
}

KOKKOS_INLINE_FUNCTION
void LorentzBoost(Real vx1, Real vx2, Real vx3, AthenaScratchTensor<Real, TensorSymm::SYM2, 4, 2> &lambda) {
  lambda.ZeroClear();
  // Velocity magnitude squared
  Real vsq = SQR(vx1)+SQR(vx2)+SQR(vx3);
  if (vsq == 0) {
    for (int a = 0; a < 4; ++a) {
     lambda(a,a) = 1;
    }
  } else {
    // Calculate Lorentz factor
    Real gamma = 1/sqrt(1-vsq);

    lambda(0,0) = gamma;
    lambda(0,1) = -gamma*vx1;
    lambda(0,2) = -gamma*vx2;
    lambda(0,3) = -gamma*vx3;

    lambda(1,1) = 1 + (gamma-1)*SQR(vx1)/vsq;
    lambda(2,2) = 1 + (gamma-1)*SQR(vx2)/vsq;
    lambda(3,3) = 1 + (gamma-1)*SQR(vx3)/vsq;
    lambda(1,2) = (gamma-1)*vx1*vx2/vsq;
    lambda(1,3) = (gamma-1)*vx1*vx3/vsq;
    lambda(2,3) = (gamma-1)*vx2*vx3/vsq;
  }
}

KOKKOS_INLINE_FUNCTION
AthenaScratchTensor<Real, TensorSymm::SYM2, 4, 2> 
inverse(AthenaScratchTensor<Real, TensorSymm::SYM2, 4, 2> matrix) {
    AthenaScratchTensor<Real, TensorSymm::SYM2, 4, 2> inv;
    Real det = 0;

    // Calculate the determinant using the formula for a 4x4 matrix
    det = matrix(0,0) * (matrix(1,1) * (matrix(2,2) * matrix(3,3) - matrix(2,3) * matrix(3,2)) -
                          matrix(1,2) * (matrix(2,1) * matrix(3,3) - matrix(2,3) * matrix(3,1)) +
                          matrix(1,3) * (matrix(2,1) * matrix(3,2) - matrix(2,2) * matrix(3,1))) -
          matrix(0,1) * (matrix(1,0) * (matrix(2,2) * matrix(3,3) - matrix(2,3) * matrix(3,2)) -
                          matrix(1,2) * (matrix(2,0) * matrix(3,3) - matrix(2,3) * matrix(3,0)) +
                          matrix(1,3) * (matrix(2,0) * matrix(3,2) - matrix(2,2) * matrix(3,0))) +
          matrix(0,2) * (matrix(1,0) * (matrix(2,1) * matrix(3,3) - matrix(2,3) * matrix(3,1)) -
                          matrix(1,1) * (matrix(2,0) * matrix(3,3) - matrix(2,3) * matrix(3,0)) +
                          matrix(1,3) * (matrix(2,0) * matrix(3,1) - matrix(2,1) * matrix(3,0))) -
          matrix(0,3) * (matrix(1,0) * (matrix(2,1) * matrix(3,2) - matrix(2,2) * matrix(3,1)) -
                          matrix(1,1) * (matrix(2,0) * matrix(3,2) - matrix(2,2) * matrix(3,0)) +
                          matrix(1,2) * (matrix(2,0) * matrix(3,1) - matrix(2,1) * matrix(3,0)));

    // Calculate the inverse using the formula for a 4x4 matrix
    inv(0,0) = (matrix(1,1) * (matrix(2,2) * matrix(3,3) - matrix(2,3) * matrix(3,2)) -
                     matrix(1,2) * (matrix(2,1) * matrix(3,3) - matrix(2,3) * matrix(3,1)) +
                     matrix(1,3) * (matrix(2,1) * matrix(3,2) - matrix(2,2) * matrix(3,1))) /
                    det;
    inv(0,1) = -(matrix(0,1) * (matrix(2,2) * matrix(3,3) - matrix(2,3) * matrix(3,2)) -
                      matrix(0,2) * (matrix(2,1) * matrix(3,3) - matrix(2,3) * matrix(3,1)) +
                      matrix(0,3) * (matrix(2,1) * matrix(3,2) - matrix(2,2) * matrix(3,1))) /
                    det;
    inv(0,2) = (matrix(0,1) * (matrix(1,2) * matrix(3,3) - matrix(1,3) * matrix(3,2)) -
                     matrix(0,2) * (matrix(1,1) * matrix(3,3) - matrix(1,3) * matrix(3,1)) +
                     matrix(0,3) * (matrix(1,1) * matrix(3,2) - matrix(1,2) * matrix(3,1))) /
                    det;
    inv(0,3) = -(matrix(0,1) * (matrix(1,2) * matrix(2,3) - matrix(1,3) * matrix(2,2)) -
                     matrix(0,2) * (matrix(1,1) * matrix(2,3) - matrix(1,3) * matrix(2,1)) +
                     matrix(0,3) * (matrix(1,1) * matrix(2,2) - matrix(1,2) * matrix(2,1))) /
                    det;

    inv(1,1) = (matrix(0,0) * (matrix(2,2) * matrix(3,3) - matrix(2,3) * matrix(3,2)) -
                     matrix(0,2) * (matrix(2,0) * matrix(3,3) - matrix(2,3) * matrix(3,0)) +
                     matrix(0,3) * (matrix(2,0) * matrix(3,2) - matrix(2,2) * matrix(3,0))) /
                    det;
    inv(1,2) = -(matrix(0,0) * (matrix(1,2) * matrix(3,3) - matrix(1,3) * matrix(3,2)) -
                      matrix(0,2) * (matrix(1,0) * matrix(3,3) - matrix(1,3) * matrix(3,0)) +
                      matrix(0,3) * (matrix(1,0) * matrix(3,2) - matrix(1,2) * matrix(3,0))) /
                    det;
    inv(1,3) = (matrix(0,0) * (matrix(1,2) * matrix(2,3) - matrix(1,3) * matrix(2,2)) -
                     matrix(0,2) * (matrix(1,0) * matrix(2,3) - matrix(1,3) * matrix(2,0)) +
                     matrix(0,3) * (matrix(1,0) * matrix(2,2) - matrix(1,2) * matrix(2,0))) /
                    det;

    inv(2,2) = (matrix(0,0) * (matrix(1,1) * matrix(3,3) - matrix(1,3) * matrix(3,1)) -
                     matrix(0,1) * (matrix(1,0) * matrix(3,3) - matrix(1,3) * matrix(3,0)) +
                     matrix(0,3) * (matrix(1,0) * matrix(3,1) - matrix(1,1) * matrix(3,0))) /
                    det;
    inv(2,3) = -(matrix(0,0) * (matrix(1,1) * matrix(2,3) - matrix(1,3) * matrix(2,1)) -
                      matrix(0,1) * (matrix(1,0) * matrix(2,3) - matrix(1,3) * matrix(2,0)) +
                      matrix(0,3) * (matrix(1,0) * matrix(2,1) - matrix(1,1) * matrix(2,0))) /
                    det;

    inv(3,3) = (matrix(0,0) * (matrix(1,1) * matrix(2,2) - matrix(1,2) * matrix(2,1)) -
                     matrix(0,1) * (matrix(1,0) * matrix(2,2) - matrix(1,2) * matrix(2,0)) +
                     matrix(0,2) * (matrix(1,0) * matrix(2,1) - matrix(1,1) * matrix(2,0))) /
                    det;
    return inv;
}
