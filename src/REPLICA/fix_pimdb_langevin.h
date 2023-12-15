/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(pimdb/langevin,FixPIMDBLangevin);

#else

#ifndef FIX_PIMDB_LANGEVIN_H
#define FIX_PIMDB_LANGEVIN_H

#include "fix_pimd_langevin.h"
#include "bosonic_exchange.h"

namespace LAMMPS_NS {

class FixPIMDBLangevin : public FixPIMDLangevin {
 public:
    FixPIMDBLangevin(class LAMMPS *, int, char **);
    ~FixPIMDBLangevin();

protected:
    void spring_force() override;

private:
  int nbosons;
  BosonicExchange bosonic_exchange;
  double** f_tag_order;

  void compute_kinetic_energy();
};


}

#endif
#endif
