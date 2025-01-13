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

   Purpose      Path Integral Molecular Dynamics of Bosons with Langevin Thermostat
   Copyright    Hirshberg lab @ Tel Aviv University
   Authors      Yotam M. Y. Feldman, Ofir Blumer

   Updated      Jan-06-2025
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

/* ---------------------------------------------------------------------- */

FixPIMDBLangevin::FixPIMDBLangevin(LAMMPS *lmp, int narg, char **arg) :
    FixPIMDLangevin(lmp, narg, arg),
    bosonic_exchange(lmp, atom->nlocal, np, universe->me, false, true)
{
    for (int i = 3; i < narg - 1; i += 2) {
        if ((strcmp(arg[i], "method") == 0) && (strcmp(arg[i+1], "pimd") != 0)) {
            error->universe_all(FLERR, "Method not supported in fix pimdb/langevin; only method PIMD");
        }
        else if (strcmp(arg[i], "scale") == 0) {
            error->universe_all(FLERR, "The scale parameter of the PILE_L thermostat is not supported for pimdb, and should be removed.");
        }
        else if ((strcmp(arg[i], "iso") == 0) || (strcmp(arg[i], "aniso") == 0) || (strcmp(arg[i], "barostat") == 0) || (strcmp(arg[i], "taup") == 0)) {
            error->universe_all(FLERR, "Barostat parameters are not available for pimdb.");
        }
    }

    if (fmmode != PHYSICAL) {
        error->universe_all(FLERR, "The only available fmmode for pimdb is physical, please remove the fmmode keyword.");
    }
    if (ensemble != NVE && ensemble != NVT) {
        error->universe_all(FLERR, "The only available ensembles for pimdb are nve and nvt, please choose one of these ensembles.");
    }

    method = PIMD;     

    size_vector = 6;

    nbosons = atom->nlocal;

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
                                              beta_np, 1 / beta, ff);

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
}

/* ---------------------------------------------------------------------- */

void FixPIMDBLangevin::compute_spring_energy() {
    total_spring_energy = bosonic_exchange.get_potential();
    se_bead = (universe->iworld == np - 1 ? total_spring_energy : 0.0);
}

/* ---------------------------------------------------------------------- */

void FixPIMDBLangevin::compute_t_prim()
{
    if (universe->iworld == 0)
        t_prim = bosonic_exchange.prim_estimator();
}

/* ---------------------------------------------------------------------- */

double FixPIMDBLangevin::compute_vector(int n)
{
    if (0 <= n && n < 6) {
        return FixPIMDLangevin::compute_vector(n);
    }
    return 0.0;
}
