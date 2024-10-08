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
  virial = 0.;
  prim = 0.;
  spring_energy = 0.;
  size_vector = 4;
  if (method != PIMD && method != NMPIMD) {
    error->universe_all(FLERR, "Method not supported in fix pimdb/nvt; only methods PIMD and NMPIMD");
  }
}

/* ---------------------------------------------------------------------- */

FixPIMDBNVT::~FixPIMDBNVT() {
}
void FixPIMDBNVT::post_force(int /*flag*/)
{
  double **x = atom->x;
  double **f = atom->f;
  for (int i = 0; i < atom->nlocal; i++)
    for (int j = 0; j < 3; j++) atom->f[i][j] /= np;

  comm_exec(atom->x);
  virial = bosonic_exchange.vir_estimator(x, f);
  if (0 == universe->me)
  {
    prim = bosonic_exchange.prim_estimator();
    spring_energy = bosonic_exchange.get_potential();
  }
  else {
    prim = -bosonic_exchange.get_spring_energy();
    spring_energy = bosonic_exchange.get_spring_energy();
  }
  spring_force(x, f);

  if (method == NMPIMD) 
  {
  /* forward comm for the force on ghost atoms */

  nmpimd_fill(atom->f);

  /* inter-partition comm */

  comm_exec(atom->f);

  /* normal-mode transform */

  nmpimd_transform(buf_beads, atom->f, M_f2fp[universe->iworld]);
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDBNVT::spring_force(double **x, double **f)
{
  double *xlast = buf_beads[x_last];
  double *xnext = buf_beads[x_next];
  double ff = fbond * atom->mass[atom->type[0]]; 
  
  bosonic_exchange.prepare_with_coordinates(*x, xlast, xnext, beta, 1 / beta, -ff);
  bosonic_exchange.spring_force(f);
}

/* ---------------------------------------------------------------------- */
double FixPIMDBNVT::compute_vector(int n)
{
  if (n == 0) return spring_energy;
  if (n == 1) return t_sys;
  if (n == 2) return virial;
  if (n == 3) return prim;
  return 0.0;
}