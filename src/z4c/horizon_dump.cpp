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
#include <cstdio>

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
HorizonDump::HorizonDump(MeshBlockPack *pmbp, ParameterInput *pin, int n, int is_common):
              common_horizon{is_common}, pos{NAN, NAN, NAN},
              pmbp{pmbp} {
  std::string nstr = std::to_string(n);

  pos[0] = pin->GetOrAddReal("z4c", "co_" + nstr + "_x", 0.0);
  pos[1] = pin->GetOrAddReal("z4c", "co_" + nstr + "_y", 0.0);
  pos[2] = pin->GetOrAddReal("z4c", "co_" + nstr + "_z", 0.0);

  horizon_extent = pin->GetOrAddReal("z4c", "co_" + nstr + "_dump_radius", 2.0);
  horizon_nx = pin->GetOrAddInteger("z4c", "horizon_"
                              + std::to_string(n)+"_Nx",10);
  Real extend[3] = {horizon_extent,horizon_extent,horizon_extent};
  int Nx[3] = {horizon_nx,horizon_nx,horizon_nx};
  pcat_grid = new CartesianGrid(pmbp, pos, extend, Nx);
  
  // Initializing variables that will be dumped
  // The order is alpha, betax, betay, betaz,
  // gxx, gxy, gxz, gyy, gyz, gzz
  // Kxx, Kxy, Kxz, Kyy, Kyz, Kzz
  variable_to_dump.push_back(std::make_pair(pmbp->pz4c->I_Z4C_ALPHA, true));
  variable_to_dump.push_back(std::make_pair(pmbp->pz4c->I_Z4C_BETAX, true));
  variable_to_dump.push_back(std::make_pair(pmbp->pz4c->I_Z4C_BETAY, true));
  variable_to_dump.push_back(std::make_pair(pmbp->pz4c->I_Z4C_BETAZ, true));
  variable_to_dump.push_back(std::make_pair(pmbp->padm->I_ADM_GXX, false));
  variable_to_dump.push_back(std::make_pair(pmbp->padm->I_ADM_GXY, false));
  variable_to_dump.push_back(std::make_pair(pmbp->padm->I_ADM_GXZ, false));
  variable_to_dump.push_back(std::make_pair(pmbp->padm->I_ADM_GYY, false));
  variable_to_dump.push_back(std::make_pair(pmbp->padm->I_ADM_GYZ, false));
  variable_to_dump.push_back(std::make_pair(pmbp->padm->I_ADM_GZZ, false));
  variable_to_dump.push_back(std::make_pair(pmbp->padm->I_ADM_KXX, false));
  variable_to_dump.push_back(std::make_pair(pmbp->padm->I_ADM_KXY, false));
  variable_to_dump.push_back(std::make_pair(pmbp->padm->I_ADM_KXZ, false));
  variable_to_dump.push_back(std::make_pair(pmbp->padm->I_ADM_KYY, false));
  variable_to_dump.push_back(std::make_pair(pmbp->padm->I_ADM_KYZ, false));
  variable_to_dump.push_back(std::make_pair(pmbp->padm->I_ADM_KZZ, false));
}

//----------------------------------------------------------------------------------------
HorizonDump::~HorizonDump() {}

void HorizonDump::SetGridAndInterpolate(Real center[NDIM]) {
  // update center location
  pos[0] = center[0];
  pos[1] = center[1];
  pos[2] = center[2];

  pcat_grid->ResetCenter(pos);
  // Real* data_out = new Real [];
  // swap out to 1d array
  Real data_out[horizon_nx][horizon_nx][horizon_nx][16];
  for(int nvar=0; nvar<16; nvar++) {
    // Interpolate here
    if (variable_to_dump[nvar].second) {
      pcat_grid->InterpolateToGrid(variable_to_dump[nvar].first,pmbp->pz4c->u0);
    } else {
      pcat_grid->InterpolateToGrid(variable_to_dump[nvar].first,pmbp->padm->u_adm);
    }
    for (int nx = 0; nx < horizon_nx; nx ++)
    for (int ny = 0; ny < horizon_nx; ny ++)
    for (int nz = 0; nz < horizon_nx; nz ++) {
      // fill data_out here
      data_out[nx][ny][nz][nvar] = pcat_grid->interp_vals.h_view(nx,ny,nz);
    }
  }

  // MPI reduce here
  // Reduction to the master rank for data_out
  int count = 16*horizon_nx*horizon_nx*horizon_nx;
  #if MPI_PARALLEL_ENABLED
  if (0 == global_variable::my_rank) {
    MPI_Reduce(MPI_IN_PLACE, &data_out[0][0][0][0], count, MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
  } else {
    MPI_Reduce(&data_out[0][0][0][0], &data_out[0][0][0][0], count, MPI_ATHENA_REAL, MPI_SUM, 0, MPI_COMM_WORLD);
  }
  #endif
  // Then write output file
  // Open the file in binary write mode
  if (0 == global_variable::my_rank) {
    FILE* etk_output_file = fopen("horizons/etk_output_file.dat", "wb");
    if (etk_output_file == nullptr) {
      perror("Error opening file");
      return;
    }
    // fwrite(&common_horizon, sizeof(int), 1, etk_output_file);
    // fwrite(&pmbp->pmesh->time, sizeof(Real), 1, etk_output_file);
    // Write the 4D array to the binary file
    size_t elementsWritten = fwrite(data_out, sizeof(Real), count, etk_output_file);
    if (elementsWritten != count) {
      perror("Error writing to file");
    }
    // Close the file
    fclose(etk_output_file);
  }
}
