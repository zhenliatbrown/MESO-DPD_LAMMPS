// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
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
   Contributing authors: Christian Trott (SNL), Stan Moore (SNL),
                         Evan Weinberg (NVIDIA)
------------------------------------------------------------------------- */

#include "compute_sna_grid_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neighbor_kokkos.h"
#include "neigh_request.h"
#include "sna.h"

#include <cmath>
#include <cstdlib>
#include <cstring>

#define MAXLINE 1024
#define MAXWORD 3

namespace LAMMPS_NS {

/* ---------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
ComputeSNAGridKokkos<DeviceType, real_type, vector_length>::ComputeSNAGridKokkos(LAMMPS *lmp, int narg, char **arg) : ComputeSNAGrid (lmp, narg, arg)
{
  //respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  //datamask_read = EMPTY_MASK;
  //datamask_modify = EMPTY_MASK;

  //host_flag = (execution_space == Host);
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

  printf("^^^ inside ComputeSNAGridKokkos init\n");
  // from pair_snap_kokkos_impl.h :
  /*
  if (host_flag) {
    if (lmp->kokkos->nthreads > 1)
      error->all(FLERR,"Pair style snap/kk can currently only run on a single "
                         "CPU thread");

    PairSNAP::init_style();
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
}

/* ----------------------------------------------------------------------
   routines used by template reference classes
------------------------------------------------------------------------- */

template<class DeviceType>
ComputeSNAGridKokkosDevice<DeviceType>::ComputeSNAGridKokkosDevice(class LAMMPS *lmp, int narg, char **arg)
   : ComputeSNAGridKokkos<DeviceType, SNAP_KOKKOS_REAL, SNAP_KOKKOS_DEVICE_VECLEN>(lmp, narg, arg) { ; }

template<class DeviceType>
void ComputeSNAGridKokkosDevice<DeviceType>::init()
{
  Base::init();
}

#ifdef LMP_KOKKOS_GPU
template<class DeviceType>
ComputeSNAGridKokkosHost<DeviceType>::ComputeSNAGridKokkosHost(class LAMMPS *lmp, int narg, char **arg)
   : ComputeSNAGridKokkos<DeviceType, SNAP_KOKKOS_REAL, SNAP_KOKKOS_HOST_VECLEN>(lmp, narg, arg) { ; }

template<class DeviceType>
void ComputeSNAGridKokkosHost<DeviceType>::init()
{
  Base::init();
}
#endif

}
