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

#ifdef FIX_CLASS
// clang-format off
FixStyle(mdr/radius/update,FixMDRradiusUpdate);
// clang-format on
#else

#ifndef LMP_FIX_MDR_RADIUS_UPDATE_H
#define LMP_FIX_MDR_RADIUS_UPDATE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMDRradiusUpdate : public Fix {
 public:
  double * Ro;
  double * Vgeo;
  double * Velas;
  double * Vcaps;
  double * eps_bar;
  double * dRnumerator; 
  double * dRdenominator;
  double * Acon0;
  double * Acon1;
  double * Atot;
  double * Atot_sum;
  double * ddelta_bar;
  double * psi;
  double * psi_b;
  double * sigmaxx;
  double * sigmayy;
  double * sigmazz;
  double * history_setup_flag; 
  double * contacts;
  double * adhesive_length;
  
  FixMDRradiusUpdate(class LAMMPS *, int, char **);
  int setmask() override;
  void setup(int) override;
  void setup_pre_force(int) override;
  void pre_force(int) override;
  void end_of_step() override; // FOR MDR
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 private:

};

}    // namespace LAMMPS_NS

#endif
#endif
