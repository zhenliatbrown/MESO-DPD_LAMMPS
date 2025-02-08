/* -*- c++ -*--------------------------------------------------------------
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
   Contributing author:  Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "region2vmd.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "region.h"

#include "region_block.h"
#include "region_cone.h"
#include "region_cylinder.h"
#include "region_sphere.h"

#include <cmath>

using namespace LAMMPS_NS;

static constexpr double SMALL = 1.0e-10;

/* ---------------------------------------------------------------------- */

void Region2VMD::command(int narg, char **arg)
{
  FILE *fp = nullptr;

  if (narg < 2) utils::missing_cmd_args(FLERR, "region2vmd", error);

  if (comm->me == 0) {
    fp = fopen(arg[0], "w");
    if (fp == nullptr) {
      error->one(FLERR, Error::ARGZERO, "Cannot open file {} for writing: {}", arg[0],
                 utils::getsyserror());
    } else {
      utils::logmesg(lmp, "Writing region visualizations to VMD Tcl script file {}:\n", arg[0]);
      fputs("# save old top molecule index\nset oldtop [molinfo top]\n", fp);
    }
  }

  for (int iarg = 1; iarg < narg; ++iarg) {
    auto *region = domain->get_region_by_id(arg[iarg]);
    if (!region) {
      if (fp) fclose(fp);
      error->all(FLERR, iarg, "Region {} does not exist", arg[iarg]);
    } else {
      write_region(fp, region);
    }
  }

  // done, close file
  if (comm->me == 0) {
    // reset views and restore previous top molecule
    fputs("display resetview\nif {$oldtop >= 0} {mol top $oldtop}\n", fp);
    fclose(fp);
  }
}

/* ----------------------------------------------------------------------
  write out one region using VMD graphics primitives
  fp is non-NULL only on MPI rank 0
   ---------------------------------------------------------------------- */

void Region2VMD::write_region(FILE *fp, Region *region)
{
  if (!fp || !region) return;

  if (region->dynamic_check()) {
    utils::logmesg(lmp, "Cannot (yet) handle moving or rotating region {}. Skipping... ",
                   region->id);
    return;
  }

  // create a new (empty) VMD molecule to hold the graphics for this specific region.

  utils::logmesg(lmp, " writing region {} ...", region->id);
  utils::print(fp, "\n# region {} of style {}\n", region->id, region->style);
  fputs("# create new empty VMD molecule to store the graphics primitives\n"
        "set gfxmol [mol new]\nmol top $gfxmol\n",
        fp);
  utils::print(fp, "mol rename $gfxmol {{LAMMPS region {}}}\n", region->id);
  fputs("# set color to desired value. Use 'silver' by default\n"
        "graphics $gfxmol color silver\n",
        fp);
  fputs("# set material to desired choice\n"
        "graphics $gfxmol material Transparent\n",
        fp);

  // translate compatible regions to VMD graphics primitives, skip others.

  const std::string regstyle = region->style;
  if (regstyle == "block") {
    const auto block = dynamic_cast<RegBlock *>(region);
    if (!block) {
      error->one(FLERR, Error::NOLASTLINE, "Region {} is not of style 'block'", region->id);
    } else {

      // a block is represented by 12 triangles

      utils::print(fp,
                   "draw triangle {{{0} {2} {4}}} {{{0} {2} {5}}} {{{0} {3} {4}}}\n"
                   "draw triangle {{{0} {3} {4}}} {{{0} {3} {5}}} {{{0} {2} {5}}}\n"
                   "draw triangle {{{0} {2} {4}}} {{{0} {2} {5}}} {{{1} {2} {4}}}\n"
                   "draw triangle {{{1} {2} {4}}} {{{1} {2} {5}}} {{{0} {2} {5}}}\n"
                   "draw triangle {{{1} {2} {4}}} {{{1} {2} {5}}} {{{1} {3} {4}}}\n"
                   "draw triangle {{{1} {3} {4}}} {{{1} {3} {5}}} {{{1} {2} {5}}}\n"
                   "draw triangle {{{0} {2} {4}}} {{{1} {2} {4}}} {{{1} {3} {4}}}\n"
                   "draw triangle {{{0} {2} {4}}} {{{0} {3} {4}}} {{{1} {3} {4}}}\n"
                   "draw triangle {{{0} {2} {5}}} {{{1} {2} {5}}} {{{1} {3} {5}}}\n"
                   "draw triangle {{{0} {2} {5}}} {{{0} {3} {5}}} {{{1} {3} {5}}}\n"
                   "draw triangle {{{0} {3} {4}}} {{{0} {3} {5}}} {{{1} {3} {5}}}\n"
                   "draw triangle {{{0} {3} {4}}} {{{1} {3} {4}}} {{{1} {3} {5}}}\n",
                   block->xlo, block->xhi, block->ylo, block->yhi, block->zlo, block->zhi);
    }

  } else if (regstyle == "cone") {
    const auto cone = dynamic_cast<RegCone *>(region);
    if (!cone) {
      error->one(FLERR, Error::NOLASTLINE, "Region {} is not of style 'cone'", region->id);
    } else {

      // The VMD cone primitive requires one radius set to zero
      if (cone->radiuslo < SMALL) {
        // a cone uses a single cone primitive
        if (cone->axis == 'x') {
          utils::print(fp, "draw cone {{{1} {2} {3}}} {{{0} {2} {3}}} radius {4} resolution 20\n",
                       cone->lo, cone->hi, cone->c1, cone->c2, cone->radiushi);
        } else if (cone->axis == 'y') {
          utils::print(fp, "draw cone {{{2} {1} {3}}} {{{2} {0} {3}}} radius {4} resolution 20\n",
                       cone->lo, cone->hi, cone->c1, cone->c2, cone->radiushi);
        } else if (cone->axis == 'z') {
          utils::print(fp, "draw cone {{{2} {3} {1}}} {{{2} {3} {0}}} radius {4} resolution 20\n",
                       cone->lo, cone->hi, cone->c1, cone->c2, cone->radiushi);
        }
      } else if (cone->radiushi < SMALL) {
        // a cone uses a single cone primitive
        if (cone->axis == 'x') {
          utils::print(fp, "draw cone {{{0} {2} {3}}} {{{1} {2} {3}}} radius {4} resolution 20\n",
                       cone->lo, cone->hi, cone->c1, cone->c2, cone->radiuslo);
        } else if (cone->axis == 'y') {
          utils::print(fp, "draw cone {{{2} {0} {3}}} {{{2} {1} {3}}} radius {4} resolution 20\n",
                       cone->lo, cone->hi, cone->c1, cone->c2, cone->radiuslo);
        } else if (cone->axis == 'z') {
          utils::print(fp, "draw cone {{{2} {3} {0}}} {{{2} {3} {1}}} radius {4} resolution 20\n",
                       cone->lo, cone->hi, cone->c1, cone->c2, cone->radiuslo);
        }
      } else {
        utils::logmesg(lmp,
                       "Cannot (yet) translate a truncated cone to VMD graphics. Skipping...\n");
      }
    }
  } else if (regstyle == "cylinder") {
    const auto cylinder = dynamic_cast<RegCylinder *>(region);
    if (!cylinder) {
      error->one(FLERR, Error::NOLASTLINE, "Region {} is not of style 'cylinder'", region->id);
    } else {
      // a cylinder uses a single cylinder primitive
      if (cylinder->axis == 'x') {
        utils::print(fp, "draw cylinder {{{0} {2} {3}}} {{{1} {2} {3}}} radius {4} resolution 20\n",
                     cylinder->lo, cylinder->hi, cylinder->c1, cylinder->c2, cylinder->radius);
      } else if (cylinder->axis == 'y') {
        utils::print(fp, "draw cylinder {{{2} {0} {3}}} {{{2} {1} {3}}} radius {4} resolution 20\n",
                     cylinder->lo, cylinder->hi, cylinder->c1, cylinder->c2, cylinder->radius);
      } else if (cylinder->axis == 'z') {
        utils::print(fp, "draw cylinder {{{2} {3} {0}}} {{{2} {3} {1}}} radius {4} resolution 20\n",
                     cylinder->lo, cylinder->hi, cylinder->c1, cylinder->c2, cylinder->radius);
      }
    }
  } else if (regstyle == "sphere") {
    const auto sphere = dynamic_cast<RegSphere *>(region);
    if (!sphere) {
      error->one(FLERR, Error::NOLASTLINE, "Region {} is not of style 'sphere'", region->id);
    } else {
      // a sphere uses a single sphere primitive
      utils::print(fp, "draw sphere {{{} {} {}}} radius {} resolution 20\n", sphere->xc, sphere->yc,
                   sphere->zc, sphere->radius);
    }
  } else {
    utils::logmesg(lmp,
                   "Cannot (yet) translate region {} of style {} to VMD graphics. Skipping... ",
                   region->id, region->style);
  }
  utils::logmesg(lmp, " done\n");
}
