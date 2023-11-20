/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Package      FixPIMDBLangevin
   Purpose      TODO
   Copyright    TODO
   Authors      TODO

   Updated      TODO
   Version      1.0
------------------------------------------------------------------------- */

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <fstream>
#include "fix_pimdb_langevin.h"
#include "universe.h"
#include "comm.h"
#include "force.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "memory.h"
#include "error.h"
#include <algorithm> 

using namespace LAMMPS_NS;
using namespace FixConst;

enum{PIMD,NMPIMD,CMD};

/* ---------------------------------------------------------------------- */

FixPIMDBLangevin::FixPIMDBLangevin(LAMMPS *lmp, int narg, char **arg) :
    FixPIMDLangevin(lmp, narg, arg),
    bosonic_exchange(lmp, atom->nlocal, np, universe->me,
                     false)
{

    // TODO: update, ensure consistency with PIMDB

  if (method == CMD) {
    error->universe_all(FLERR, "Method cmd not supported in fix pimdb");
  }

  nbosons    = atom->nlocal;
}

/* ---------------------------------------------------------------------- */

FixPIMDBLangevin::~FixPIMDBLangevin() {
}

/* ---------------------------------------------------------------------- */

void FixPIMDBLangevin::spring_force() {
    double ff = fbond * atom->mass[atom->type[0]]; // TODO: ensure that all masses are the same
    // TODO: verify these pointers
    int nlocal = atom->nlocal;
    double* me_bead_positions = *(atom->x);
    double* last_bead_positions = bufsortedall[x_last * nlocal];
    double* next_bead_positions = bufsortedall[x_next * nlocal];
    bosonic_exchange.prepare_with_coordinates(me_bead_positions,
                                              last_bead_positions, next_bead_positions,
                                              beta, ff);

    // TODO: translate order back from tag order
    // TODO virial
    double virial = bosonic_exchange.spring_force(atom->f);

    // TODO spring force
    // if (universe->me == np - 1) {
    //     spring_energy = bosonic_exchange.get_potential();
    // }
}

/* ---------------------------------------------------------------------- */