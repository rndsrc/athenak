#ifndef MHD_MHD_HPP_
#define MHD_MHD_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd.hpp
//  \brief definitions for MHD class

#include "athena.hpp"
#include "parameter_input.hpp"
#include "tasklist/task_list.hpp"
#include "bvals/bvals.hpp"

// forward declarations
class Driver;
class EquationOfState;

// constants that enumerate MHD Riemann Solver options
enum MHD_RSolver {advect, llf, hlld, roe};

namespace mhd {

//----------------------------------------------------------------------------------------
//! \class MHD

class MHD
{
 public:
  MHD(MeshBlockPack *ppack, ParameterInput *pin);
  ~MHD();

  // data
  EquationOfState *peos;    // object that implements chosen EOS

  int nhydro;             // number of hydro variables (5/4 for adiabatic/isothermal)
  int nscalars;           // number of passive scalars
  DvceArray5D<Real> u0;   // conserved variables
  DvceArray5D<Real> w0;   // primitive variables
  FaceArray4D<Real> b0;   // face-centered magnetic fields

  // Object containing boundary communication buffers and routines
  BoundaryValues *pbvals;

  // following only used for time-evolving flow
  DvceArray5D<Real> u1;           // conserved variables at intermediate step 
  DvceArray5D<Real> divf;         // divergence of fluxes
  FaceArray4D<Real> b1;           // face-centered magnetic fields at intermediate step
  DvceArray3D<Real> uflx_x1face;  // fluxes on x1-faces
  DvceArray3D<Real> uflx_x2face;  // fluxes on x2-faces
  DvceArray3D<Real> uflx_x3face;  // fluxes on x3-faces
  Real dtnew;

  // functions
  void MHDStageStartTasks(TaskList &tl, TaskID start, std::vector<TaskID> &added);
  void MHDStageRunTasks(TaskList &tl, TaskID start, std::vector<TaskID> &added);
  void MHDStageEndTasks(TaskList &tl, TaskID start, std::vector<TaskID> &added);
  TaskStatus MHDInitRecv(Driver *d, int stage);
  TaskStatus MHDClearRecv(Driver *d, int stage);
  TaskStatus MHDClearSend(Driver *d, int stage);
  TaskStatus MHDCopyCons(Driver *d, int stage);
  TaskStatus MHDDivFlux(Driver *d, int stage);
  TaskStatus MHDUpdate(Driver *d, int stage);
  TaskStatus MHDSend(Driver *d, int stage); 
  TaskStatus MHDReceive(Driver *d, int stage); 
  TaskStatus ConToPrim(Driver *d, int stage);
  TaskStatus NewTimeStep(Driver *d, int stage);
  TaskStatus MHDApplyPhysicalBCs(Driver* pdrive, int stage);

  // functions to set physical BCs for Hydro conserved variables, applied to single MB
  // specified by argument 'm'. 
  void ReflectInnerX1(int m);
  void ReflectOuterX1(int m);
  void ReflectInnerX2(int m);
  void ReflectOuterX2(int m);
  void ReflectInnerX3(int m);
  void ReflectOuterX3(int m);
  void OutflowInnerX1(int m);
  void OutflowOuterX1(int m);
  void OutflowInnerX2(int m);
  void OutflowOuterX2(int m);
  void OutflowInnerX3(int m);
  void OutflowOuterX3(int m);

 private:
  MeshBlockPack* pmy_pack;   // ptr to MeshBlockPack containing this MHD
  ReconstructionMethod recon_method_;
  MHD_RSolver rsolver_method_;
};

} // namespace mhd
#endif // MHD_MHD_HPP_
