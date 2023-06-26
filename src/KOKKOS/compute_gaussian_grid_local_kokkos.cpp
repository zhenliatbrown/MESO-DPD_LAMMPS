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

  k_cutsq = tdual_fparams("ComputeSNAGridKokkos::cutsq",atom->ntypes+1,atom->ntypes+1);
  auto d_cutsq = k_cutsq.template view<DeviceType>();
  rnd_cutsq = d_cutsq;

  host_flag = (execution_space == Host);

  // TODO: Extract cutsq in double loop below, no need for cutsq_tmp

  //cutsq_tmp = cutsq[1][1];

  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = 1; j <= atom->ntypes; j++){
      k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutsq[i][j]; //cutsq_tmp;
      k_cutsq.template modify<LMPHostType>();
    }
  }
  //printf(">>> 1\n");
  // Set up element lists
  MemKK::realloc_kokkos(d_radelem,"ComputeSNAGridKokkos::radelem",nelements);
  int n = atom->ntypes;
  MemKK::realloc_kokkos(d_map,"ComputeSNAGridKokkos::map",n+1);
  //printf(">>> 2\n");
  auto h_radelem = Kokkos::create_mirror_view(d_radelem);
  auto h_map = Kokkos::create_mirror_view(d_map);
  //printf(">>> 3\n");
  // start from index 1 because of how compute sna/grid is
  for (int i = 1; i <= atom->ntypes; i++) {
    h_radelem(i-1) = radelem[i];
  }
  //printf(">>> 4\n");
  // In pair snap some things like `map` get allocated regardless of chem flag.
  // In this compute, however, map does not get allocated in parent classes.
  /*
  for (int i = 1; i <= atom->ntypes; i++) {
    h_map(i) = map[i];
  }
  */
  //printf(">>> 5\n");
  Kokkos::deep_copy(d_radelem,h_radelem);
  Kokkos::deep_copy(d_map,h_map);
  //printf(">>> 6\n");

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeGaussianGridLocalKokkos<DeviceType>::~ComputeGaussianGridLocalKokkos()
{
  printf(">>> ComputeGaussianGridLocalKokkos destruct begin, copymode %d\n", copymode);
  if (copymode) return;

  memoryKK->destroy_kokkos(k_cutsq,cutsq);
  memoryKK->destroy_kokkos(k_alocal,alocal);
  //gridlocal_allocated = 0;

  printf(">>> ComputeGaussianGridLocalKokkos end\n");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeGaussianGridLocalKokkos<DeviceType>::setup()
{

  // Do not call ComputeGrid::setup(), we don't wanna allocate the grid array there.
  // Instead, call ComputeGrid::set_grid_global and set_grid_local to set the n indices.

  //ComputeGrid::set_grid_global();
  //ComputeGrid::set_grid_local();
  ComputeGridLocal::setup();

  // allocate arrays
  printf(">>> rows cols kokkos init: %d %d\n", size_local_rows, size_local_cols);
  memoryKK->create_kokkos(k_alocal, alocal, size_local_rows, size_local_cols, "grid:alocal");

  //gridlocal_allocated = 1;
  //array = gridall;

  d_alocal = k_alocal.template view<DeviceType>();
  //d_grid = k_grid.template view<DeviceType>();
  //d_gridall = k_gridall.template view<DeviceType>();

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
  printf(">>> compute_local Kokkos begin\n");

  if (host_flag) {
    return;
  }

  invoked_local = update->ntimestep;

  copymode = 1;

  zlen = nzhi-nzlo+1;
  ylen = nyhi-nylo+1;
  xlen = nxhi-nxlo+1;
  total_range = (nzhi-nzlo+1)*(nyhi-nylo+1)*(nxhi-nxlo+1);

  atomKK->sync(execution_space,X_MASK|F_MASK|TYPE_MASK);
  x = atomKK->k_x.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  k_cutsq.template sync<DeviceType>();

  // max_neighs is defined here - think of more elaborate methods.
  max_neighs = 100;

  // Pair snap/kk uses grow_ij with some max number of neighs but compute sna/grid uses total 
  // number of atoms.

  ntotal = atomKK->nlocal + atomKK->nghost;
  // Allocate view for number of neighbors per grid point
  MemKK::realloc_kokkos(d_ninside,"ComputeSNAGridKokkos:ninside",total_range);

  // "chunksize" variable is default 32768 in compute_sna_grid.cpp, and set by user 
  // `total_range` is the number of grid points which may be larger than chunk size.
  //printf(">>> total_range: %d\n", total_range);
  chunksize = 32768;
  chunk_size = MIN(chunksize, total_range);
  chunk_offset = 0;

  int vector_length_default = 1;
  int team_size_default = 1;
  if (!host_flag)
    team_size_default = 32;//max_neighs;

  while (chunk_offset < total_range) { // chunk up loop to prevent running out of memory

    if (chunk_size > total_range - chunk_offset)
      chunk_size = total_range - chunk_offset;

    //Neigh
    {
      int vector_length = vector_length_default;
      int team_size = team_size_default;
      check_team_size_for<TagComputeGaussianGridLocalNeigh>(chunk_size,team_size,vector_length);
      printf(">>> Check 1 %d %d %d\n", chunk_size, team_size, vector_length);
      typename Kokkos::TeamPolicy<DeviceType, TagComputeGaussianGridLocalNeigh> policy_neigh(chunk_size,team_size,vector_length);
      printf(">>> Check 2\n");
      Kokkos::parallel_for("ComputeGaussianGridLocalNeigh",policy_neigh,*this);
    }

    // Proceed to the next chunk.
    chunk_offset += chunk_size;
  } // end while

  copymode = 0;

  k_alocal.template modify<DeviceType>();
  k_alocal.template sync<LMPHostType>();

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void ComputeGaussianGridLocalKokkos<DeviceType>::operator() (TagComputeGaussianGridLocalNeigh,const typename Kokkos::TeamPolicy<DeviceType, TagComputeGaussianGridLocalNeigh>::member_type& team) const
{
  const int ii = team.league_rank();
  //printf("%d\n", ii);
}

/* ----------------------------------------------------------------------
   check max team size
------------------------------------------------------------------------- */

template<class DeviceType>
template<class TagStyle>
void ComputeGaussianGridLocalKokkos<DeviceType>::check_team_size_for(int inum, int &team_size, int vector_length) {
  int team_size_max;

  team_size_max = Kokkos::TeamPolicy<DeviceType,TagStyle>(inum,Kokkos::AUTO).team_size_max(*this,Kokkos::ParallelForTag());

  if (team_size*vector_length > team_size_max)
    team_size = team_size_max/vector_length;
}

namespace LAMMPS_NS {
template class ComputeGaussianGridLocalKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeGaussianGridLocalKokkos<LMPHostType>;
#endif
}