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
FixStyle(mdr/mean/surf/disp,FixMDRmeanSurfDisp);
// clang-format on
#else

#ifndef LMP_FIX_MDR_MEAN_SURF_DISP_H
#define LMP_FIX_MDR_MEAN_SURF_DISP_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMDRmeanSurfDisp : public Fix {
 public:           
  double * Acon0; 
  double * ddelta_bar;
  
  FixMDRmeanSurfDisp(class LAMMPS *, int, char **);
  int setmask() override;
  void setup_pre_force(int) override;
  void pre_force(int) override; // FOR MDR
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 private:

};

}    // namespace LAMMPS_NS

#endif
#endif
