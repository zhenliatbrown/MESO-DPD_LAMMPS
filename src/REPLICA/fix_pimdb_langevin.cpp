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
  if (method != PIMD) {
    error->universe_all(FLERR, "Method not supported in fix pimdb/langevin; only method PIMD");
  }

  nbosons    = atom->nlocal;

  memory->create(f_tag_order, nbosons, 3, "FixPIMDBLangevin:f_tag_order");
}

/* ---------------------------------------------------------------------- */

FixPIMDBLangevin::~FixPIMDBLangevin() {
    memory->destroy(f_tag_order);
}

/* ---------------------------------------------------------------------- */

void FixPIMDBLangevin::spring_force() {
    double ff = fbond * atom->mass[atom->type[0]]; // TODO: ensure that all masses are the same
    int nlocal = atom->nlocal;
    double* me_bead_positions = *(atom->x);
    double* last_bead_positions = &bufsortedall[x_last * nlocal][0];
    double* next_bead_positions = &bufsortedall[x_next * nlocal][0];

    bosonic_exchange.prepare_with_coordinates(me_bead_positions,
                                              last_bead_positions, next_bead_positions,
                                              beta, ff);

    for (int i = 0; i < nbosons; i++) {
        f_tag_order[i][0] = 0.0;
        f_tag_order[i][1] = 0.0;
        f_tag_order[i][2] = 0.0;
    }
    bosonic_exchange.spring_force(f_tag_order);

    double** f = atom->f;
    tagint* tag = atom->tag;
    for (int i = 0; i < nbosons; i++) {
        f[i][0] += f_tag_order[tag[i] - 1][0];
        f[i][1] += f_tag_order[tag[i] - 1][1];
        f[i][2] += f_tag_order[tag[i] - 1][2];
    }

     total_spring_energy = bosonic_exchange.get_potential();
}