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
FixStyle(pimdb/nvt,FixPIMDBNVT);
// clang-format on
#else

#ifndef FIX_PIMDB_NVT_H
#define FIX_PIMDB_NVT_H

#include "fix_pimd_nvt.h"
#include "bosonic_exchange.h"

namespace LAMMPS_NS {

class FixPIMDBNVT : public FixPIMDNVT {
 public:
    FixPIMDBNVT(class LAMMPS *, int, char **);
    ~FixPIMDBNVT();
   //  void post_force(int) override;
    double compute_vector(int) override;

 protected:
    void spring_force() override;
    void kinetic_estimators() override;

 private:
    BosonicExchange bosonic_exchange;
    double prim;
};

}    // namespace LAMMPS_NS

#endif
#endif
