/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Package      FixPIMDNVT
   Purpose      Quantum Path Integral Algorithm for Quantum Chemistry
   Copyright    Voth Group @ University of Chicago
   Authors      Chris Knight & Yuxing Peng (yuxing at uchicago.edu)

   Updated      Oct-01-2011
   Version      1.0
------------------------------------------------------------------------- */

#include "fix_pimdb_nvt.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "math_special.h"
#include "memory.h"
#include "universe.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

using MathSpecial::powint;

enum{PIMD,NMPIMD,CMD};

/* ---------------------------------------------------------------------- */

FixPIMDBNVT::FixPIMDBNVT(LAMMPS *lmp, int narg, char **arg) :
    FixPIMDNVT(lmp, narg, arg),
    bosonic_exchange(lmp, atom->nlocal, np, universe->me,false)
{
  beta = 1.0 / force->boltz / nhc_temp;
  t_prim = 0.;
  virial = virial2 = 0.;
  size_vector = 4;
  if (method != PIMD) {
    error->universe_all(FLERR, "Method not supported in fix pimdb/nvt; only method PIMD");
  }
}

/* ---------------------------------------------------------------------- */

FixPIMDBNVT::~FixPIMDBNVT() {
}
void FixPIMDBNVT::post_force(int /*flag*/)
{
  for (int i = 0; i < atom->nlocal; i++)
    for (int j = 0; j < 3; j++) atom->f[i][j] /= np;

  comm_exec(atom->x);
  virial = bosonic_exchange.vir_estimator(atom->x, atom->f);
  spring_force();
  virial2 = bosonic_exchange.vir_estimator(atom->x, atom->f);
}

/* ---------------------------------------------------------------------- */

void FixPIMDBNVT::spring_force()
{
  double **x = atom->x;
  double **f = atom->f;
  double *xlast = buf_beads[x_last];
  double *xnext = buf_beads[x_next];
  double ff = fbond * atom->mass[atom->type[0]]; 
  
  bosonic_exchange.prepare_with_coordinates(*x, xlast, xnext, beta, -ff);
  bosonic_exchange.spring_force(f);
}

double FixPIMDBNVT::compute_vector(int n)
{
  if (n == 0) return bosonic_exchange.get_potential();
  if (n == 1) return t_sys;
  if (n == 2) return virial;
  if (n == 3) return virial2; // bosonic_exchange.prim_estimator();
  return 0.0;
}