// clang-format off
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

#include "compute_sna_grid_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor_kokkos.h"
#include "sna_kokkos.h"
#include "update.h"

using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
ComputeSNAGridKokkos<DeviceType, real_type, vector_length>::ComputeSNAGridKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeSNAGrid(lmp, narg, arg)
{

  printf("^^^ inside ComputeSNAGridKokkos constructor\n");
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  host_flag = (execution_space == Host);

}

/* ---------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
ComputeSNAGridKokkos<DeviceType, real_type, vector_length>::~ComputeSNAGridKokkos()
{
  if (copymode) return;


}

/* ---------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
void ComputeSNAGridKokkos<DeviceType, real_type, vector_length>::init()
{

  printf("^^^ beginning of ComputeSNAGridKokkos init()\n");

  // init non-kk compute
  // this calls snaptr->init(), we probably want to init the kokkos snaptr?
  // let's copy pair_snap_kokkos by making a snaKK in header
  ComputeSNAGrid::init();

  // adjust neighbor list request for KOKKOS

  // taken from compute_coord_atom_kokkos
  // this segfaults
  /*
  printf("^^^ before neigh request\n");
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same<DeviceType,LMPHostType>::value &&
                           !std::is_same<DeviceType,LMPDeviceType>::value);
  request->set_kokkos_device(std::is_same<DeviceType,LMPDeviceType>::value);
  */


  // taken from pair_snap_kokkos init
  // compile errors with:
  // error: pointer to incomplete class type "LAMMPS_NS::KokkosLMP" is not allowed"
  /*
  if (host_flag) {
    if (lmp->kokkos->nthreads > 1)
      error->all(FLERR,"compute sna grid can currently only run on a single "
                         "CPU thread");

    // this calls snaptr->init()
    // we probably wanna call init of kokkos snaptr
    ComputeSNAGrid::init();
    return;
  }

  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style SNAP requires newton pair on");

  // neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;

  auto request = neighbor->add_request(this, NeighConst::REQ_FULL);
  request->set_kokkos_host(std::is_same<DeviceType,LMPHostType>::value &&
                           !std::is_same<DeviceType,LMPDeviceType>::value);
  request->set_kokkos_device(std::is_same<DeviceType,LMPDeviceType>::value);
  if (neighflag == FULL)
    error->all(FLERR,"Must use half neighbor list style with pair snap/kk");
  */

  // Overall, I think maybe this compute does not need a neighlist request because the original
  // compute_sna_grid.cpp does not have one.
}

namespace LAMMPS_NS {
template class ComputeSNAGridKokkos<LMPDeviceType, real_type, vector_length>;
#ifdef LMP_KOKKOS_GPU
template class ComputeSNAGridKokkos<LMPHostType, real_type, vector_length>;
#endif
}
