//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================

#include <assert.h>
#include <unistd.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

#include "horizon_dump.hpp"
#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "parameter_input.hpp"
#include "utils/cart_grid.hpp"
#include "coordinates/adm.hpp"
#include "mhd/mhd.hpp"
#include "z4c/z4c.hpp"

//----------------------------------------------------------------------------------------
HorizonDump::HorizonDump(MeshBlockPack *pmbp, ParameterInput *pin, int n, bool common_horizon):
              common_horizon{false}, pos{NAN, NAN, NAN},
              pmbp{pmbp} {
  std::string nstr = std::to_string(n);

  pos[0] = pin->GetOrAddReal("z4c", "co_" + nstr + "_x", 0.0);
  pos[1] = pin->GetOrAddReal("z4c", "co_" + nstr + "_y", 0.0);
  pos[2] = pin->GetOrAddReal("z4c", "co_" + nstr + "_z", 0.0);

  radius = pin->GetOrAddReal("z4c", "co_" + nstr + "_radius", 0.0);
  horizon_nx = pin->GetOrAddInteger("z4c", "horizon_"
                              + std::to_string(n)+"_Nx",100);
  Real extend[3] = {radius,radius,radius};
  int Nx[3] = {horizon_nx,horizon_nx,horizon_nx};
  pcat_grid = new CartesianGrid(pmbp, pos, extend, Nx);
}

//----------------------------------------------------------------------------------------
HorizonDump::~HorizonDump() {}

void HorizonDump::SetGridAndInterpolate(Real center[NDIM]) {
  // update center location
  pos[0] = center[0];
  pos[1] = center[1];
  pos[2] = center[2];

  pcat_grid->ResetCenter(pos);

  Real data_out[horizon_nx][horizon_nx][horizon_nx][16];
  for(int nvar=0; nvar<16; nvar++) {
    // Interpolate here
    for (int nx = 0; nx < horizon_nx; nx ++)
    for (int ny = 0; ny < horizon_nx; ny ++)
    for (int nz = 0; nz < horizon_nx; nz ++) {
      // fill data_out here
    }
  }

  // MPI reduce here

  // Then write output file
  WriteFile();
}

void HorizonDump::WriteFile() {

}
