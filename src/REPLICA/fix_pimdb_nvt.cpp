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
#include "error.h"
#include "force.h"
#include "universe.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixPIMDBNVT::FixPIMDBNVT(LAMMPS *lmp, int narg, char **arg) :
    FixPIMDNVT(lmp, narg, arg),
    // CR: apply mic (compatible with previous implementation, and with pimd_nvt)
    // OB: Did you mean apply minimum image convention? Because "bosonic_exchange" takes care of it...
    bosonic_exchange(lmp, atom->nlocal, np, universe->me, false, false)
{
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

/* ---------------------------------------------------------------------- */

void FixPIMDBNVT::estimate_energies()
{
  vir_estimator();
  if (universe->me == 0)
  {
    prim = bosonic_exchange.prim_estimator();
    spring_energy = bosonic_exchange.get_potential();
  }
  else {
    double total_spring_energy_for_bead = bosonic_exchange.get_total_spring_energy_for_bead();
    prim = -total_spring_energy_for_bead;
    spring_energy = total_spring_energy_for_bead;
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMDBNVT::spring_force()
{
  double **x = atom->x;
  double **f = atom->f;

  double *xlast = buf_beads[x_last];
  double *xnext = buf_beads[x_next];
  double ff = fbond * atom->mass[atom->type[0]]; 
  
  bosonic_exchange.prepare_with_coordinates(*x, xlast, xnext, beta, 1 / beta, -ff);
  bosonic_exchange.spring_force(f);
}

/* ---------------------------------------------------------------------- */

double FixPIMDBNVT::compute_vector(int n)
{
    if (0 <= n && n < 3) {
        return FixPIMDNVT::compute_vector(n);
    }
  // CR: needs to be added also to the documentation.
  // CR: Reminds that we need to add documentation about the entire bosonic fix
  if (n == 3) return prim;
  return 0.0;
}
