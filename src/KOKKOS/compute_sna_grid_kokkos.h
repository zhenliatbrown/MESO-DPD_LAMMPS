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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(sna/grid/kk,ComputeSNAGridKokkosDevice<LMPDeviceType>);
ComputeStyle(sna/grid/kk/device,ComputeSNAGridKokkosDevice<LMPDeviceType>);
#ifdef LMP_KOKKOS_GPU
ComputeStyle(sna/grid/kk/host,ComputeSNAGridKokkosHost<LMPHostType>);
#else
ComputeStyle(sna/grid/kk/host,ComputeSNAGridKokkosDevice<LMPHostType>);
#endif
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_SNA_GRID_KOKKOS_H
#define LMP_COMPUTE_SNA_GRID_KOKKOS_H

#include "compute_sna_grid.h"
#include "kokkos_type.h"
//#include "neigh_list_kokkos.h"
#include "sna_kokkos.h"
//#include "pair_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType, typename real_type_, int vector_length_>
class ComputeSNAGridKokkos : public ComputeSNAGrid {
 public:
  //enum {EnabledNeighFlags=FULL|HALF|HALFTHREAD};
  //enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  static constexpr int vector_length = vector_length_;
  using real_type = real_type_;
  //using complex = SNAComplex<real_type>;

  ComputeSNAGridKokkos(class LAMMPS *, int, char **);
  ~ComputeSNAGridKokkos() override;

  void init() override;
  //void compute_array(int, int) override;
  //double memory_usage() override;
  
 protected:
  SNAKokkos<DeviceType, real_type, vector_length> snaKK;


  using KKDeviceType = typename KKDevice<DeviceType>::value;

};

// These wrapper classes exist to make the compute style factory happy/avoid having
// to extend the compute style factory to support Compute classes w/an arbitrary number
// of extra template parameters

template <class DeviceType>
class ComputeSNAGridKokkosDevice : public ComputeSNAGridKokkos<DeviceType, SNAP_KOKKOS_REAL, SNAP_KOKKOS_DEVICE_VECLEN> {

 private:
  using Base = ComputeSNAGridKokkos<DeviceType, SNAP_KOKKOS_REAL, SNAP_KOKKOS_DEVICE_VECLEN>;

 public:

  ComputeSNAGridKokkosDevice(class LAMMPS *, int, char **);
  //ComputeSNAGridKokkosDevice(class LAMMPS *);

  void init() override;
  //double memory_usage() override;

};

#ifdef LMP_KOKKOS_GPU
template <class DeviceType>
class ComputeSNAGridKokkosHost : public ComputeSNAGridKokkos<DeviceType, SNAP_KOKKOS_REAL, SNAP_KOKKOS_HOST_VECLEN> {

 private:
  using Base = ComputeSNAGridKokkos<DeviceType, SNAP_KOKKOS_REAL, SNAP_KOKKOS_HOST_VECLEN>;

 public:

  ComputeSNAGridKokkosHost(class LAMMPS *, int, char **);
  //ComputeSNAGridKokkosHost(class LAMMPS *);

  void init();
  //double memory_usage();

};
#endif

}

#endif
#endif
