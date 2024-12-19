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
FixStyle(GRANULAR/MDR,FixGranularMDR);
// clang-format on
#else

#ifndef LMP_FIX_GRANULAR_MDR_H
#define LMP_FIX_GRANULAR_MDR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGranularMDR : public Fix {
 public:
  FixGranularMDR(class LAMMPS *, int, char **);
  ~FixGranularMDR() override;
  int setmask() override;
  void post_constructor() override;
  void setup(int) override;
  void setup_pre_force(int) override;
  void pre_force(int) override;
  void end_of_step() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 private:
  int comm_stage;
  char *id_fix;

  void radius_update();
  void mean_surf_disp();
  void calculate_contact_penalty();
  void update_fix_gran_wall(Fix*);

  double *Ro;                 // initial radius
  double *Vgeo;               // geometric particle volume of apparent particle afterremoving spherical cap volume
  double *Velas;              // particle volume from linear elasticity
  double *Vcaps;              // spherical cap volume from intersection of apparentradius particle and contact planes
  double *eps_bar;            // volume-averaged infinitesimal strain tensor
  double *dRnumerator;        // summation of numerator terms in calculation of dR
  double *dRdenominator;      // summation of denominator terms in calculation of dR
  double *Acon0;              // total area involved in contacts: Acon^{n}
  double *Acon1;              // total area involved in contacts: Acon^{n+1}
  double *Atot;               // total particle area
  double *Atot_sum;           // running sum of contact area minus cap area
  double *ddelta_bar;         // change in mean surface displacement
  double *psi;                // ratio of free surface area to total surface area
  double *psi_b;              // TEMPORARY, SINCE PSI_B IS ALREADY DEFINED IN THEINPUT SCRIPT
  double *sigmaxx;            // xx-component of the stress tensor, not necessary forforce calculation
  double *sigmayy;            // yy-component of the stress tensor, not necessary forforce calculation
  double *sigmazz;            // zz-component of the stress tensor, not necessary forforce calculation
  double *history_setup_flag; // flag to check if history variables have beeninitialized
  double *contacts;           // total contacts on particle
  double *adhesive_length;    // total length of adhesive contact on a particle
};

}    // namespace LAMMPS_NS

#endif
#endif
