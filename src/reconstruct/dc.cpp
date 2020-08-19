//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file dc.cpp
//  \brief piecewise constant (donor cell) reconstruction implemented in a derived class

#include <iostream>
#include "athena.hpp"
#include "athena_arrays.hpp"
#include "mesh/mesh.hpp"
#include "reconstruct.hpp"

//----------------------------------------------------------------------------------------
// DonorCell constructor

DonorCell::DonorCell(std::unique_ptr<ParameterInput> &pin) : Reconstruction(pin) {
  
}

//----------------------------------------------------------------------------------------
//! \fn DonorCell::ReconstructX1()
//  \brief reconstruct L/R surfaces of the i-th cells

void DonorCell::ReconstructX1(const int il, const int iu, const AthenaArray<Real> &q,
                              AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  int nvar = q.GetDim(2);
  for (int n=0; n<nvar; ++n) {
    for (int i=il; i<=iu; ++i) {
      ql(n,i+1) = q(n,i);
      qr(n,i  ) = q(n,i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn DonorCell::ReconstructX2()
//  \brief


void DonorCell::ReconstructX2(const int il, const int iu, const AthenaArray<Real> &q,
                              AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  int nvar = q.GetDim(2);
  for (int n=0; n<nvar; ++n) {
    for (int i=il; i<=iu; ++i) {
      ql(n,i) = q(n,i);
      qr(n,i) = q(n,i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn DonorCell::ReconstructX3()
//  \brief

void DonorCell::ReconstructX3(const int il, const int iu, const AthenaArray<Real> &q,
                                 AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  int nvar = q.GetDim(2);
  for (int n=0; n<nvar; ++n) {
    for (int i=il; i<=iu; ++i) {
      ql(n,i) = q(n,i);
      qr(n,i) = q(n,i);
    }
  }
  return;
}