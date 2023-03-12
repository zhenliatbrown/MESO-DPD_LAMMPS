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
ComputeStyle(sna/grid/kk,ComputeSNAGridKokkos<LMPDeviceType>);
//ComputeStyle(sna/grid/kk/device,ComputeSNAGridKokkos<LMPDeviceType>);
//ComputeStyle(sna/grid/kk/host,ComputeSNAGridKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_SNA_GRID_KOKKOS_H
#define LMP_COMPUTE_SNA_GRID_KOKKOS_H

#include "compute_sna_grid.h"
#include "kokkos_type.h"
#include "sna_kokkos.h"

namespace LAMMPS_NS {

//template<int CSTYLE, int NCOL>
//struct TagComputeCoordAtom{};

// copying pair_snap_kokkos, template args are real_type and vector_length
template<class DeviceType, typename real_type_, int vector_length_>
class ComputeSNAGridKokkos : public ComputeSNAGrid {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  static constexpr int vector_length = vector_length_;
  using real_type = real_type_;

  ComputeSNAGridKokkos(class LAMMPS *, int, char **);
  ~ComputeSNAGridKokkos() override;
  void init() override;
  //void compute_peratom() override;
  //enum {NONE,CUTOFF,ORIENT};

  //template<int CSTYLE, int NCOL>
  //KOKKOS_INLINE_FUNCTION
  //void operator()(TagComputeCoordAtom<CSTYLE,NCOL>, const int&) const;
  
 protected:

  // these are used by pair_snap_kokkos
  // neighflag gets set in init()
  // what about host_flag?
  // dunno... commented these out for now
  int host_flag, neighflag;

  SNAKokkos<DeviceType, real_type, vector_length> snaKK;

 private:


  /*
  int inum;

  typename AT::t_x_array_randomread x;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread type;
  typename ArrayTypes<DeviceType>::t_int_1d mask;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  typename AT::t_int_1d d_typelo;
  typename AT::t_int_1d d_typehi;

  DAT::tdual_float_1d k_cvec;
  typename AT::t_float_1d d_cvec;
  DAT::tdual_float_2d k_carray;
  typename AT::t_float_2d d_carray;

  typename AT::t_float_2d d_normv;
  */
};

}

#endif
#endif

