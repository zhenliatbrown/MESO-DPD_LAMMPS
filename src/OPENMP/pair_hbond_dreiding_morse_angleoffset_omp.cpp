// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U), Don Xu/EiPi Fun
------------------------------------------------------------------------- */

#include "pair_hbond_dreiding_morse_angleoffset_omp.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "math_const.h"
#include "math_special.h"
#include "molecule.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathSpecial;

static constexpr double SMALL = 0.001;

/* ---------------------------------------------------------------------- */

PairHbondDreidingMorseAngleoffsetOMP::PairHbondDreidingMorseAngleoffsetOMP(LAMMPS *lmp) :
  PairHbondDreidingMorseOMP(lmp) {
  angle_offset_flag = 1;
}
