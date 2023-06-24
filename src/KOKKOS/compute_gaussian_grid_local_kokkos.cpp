/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Drew Rohskopf (SNL)
------------------------------------------------------------------------- */

#include "compute_gaussian_grid_local_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor_kokkos.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeGaussianGridLocalKokkos<DeviceType>::ComputeGaussianGridLocalKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeGaussianGridLocal(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeGaussianGridLocalKokkos<DeviceType>::~ComputeGaussianGridLocalKokkos()
{
  if (copymode) return;

  //memoryKK->destroy_kokkos(k_result,result);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeGaussianGridLocalKokkos<DeviceType>::init()
{
  ComputeGaussianGridLocal::init();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeGaussianGridLocalKokkos<DeviceType>::compute_local()
{

  printf(">>> compute_local Kokkos\n");

}

namespace LAMMPS_NS {
template class ComputeGaussianGridLocalKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeGaussianGridLocalKokkos<LMPHostType>;
#endif
}