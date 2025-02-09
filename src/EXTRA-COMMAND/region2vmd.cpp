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

#include <cstring>
#include <unordered_set>

using namespace LAMMPS_NS;

static constexpr double SMALL = 1.0e-10;

static const std::unordered_set<std::string> vmdcolors{
    "blue",    "red",      "gray",   "orange", "yellow",  "tan",    "silver",  "green",  "white",
    "pink",    "cyan",     "purple", "lime",   "mauve",   "ochre",  "iceblue", "black",  "yellow2",
    "yellow3", "green2",   "green3", "cyan2",  "cyan3",   "blue2",  "blue3",   "violet", "violet2",
    "magenta", "magenta2", "red2",   "red3",   "orange2", "orange3"};

static const std::unordered_set<std::string> vmdmaterials{
    "Opaque",      "Transparent", "BrushedMetal", "Diffuse",     "Ghost",          "Glass1",
    "Glass2",      "Glass3",      "Glossy",       "HardPlastic", "MetallicPastel", "Steel",
    "Translucent", "Edgy",        "EdgyShiny",    "EdgyGlass",   "Goodsell",       "AOShiny",
    "AOChalky",    "AOEdgy",      "BlownGlass",   "GlassBubble", "RTChrome"};

// class that "owns" the file pointer and closes it when going out of scope.
// this avoids a lot of redundant checks and calls.
class AutoClose {
 public:
  AutoClose() = delete;
  AutoClose(const AutoClose &) = delete;
  AutoClose(const AutoClose &&) = delete;
  explicit AutoClose(FILE *_fp) : fp(_fp) {};
  ~AutoClose()
  {
    if (fp) fclose(fp);
  }

 private:
  FILE *fp;
};

/* ---------------------------------------------------------------------- */

void Region2VMD::command(int narg, char **arg)
{
  if (narg < 3) utils::missing_cmd_args(FLERR, "region2vmd", error);

  FILE *fp = nullptr;
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

  // automatically close fp when fpowner goes out of scope
  AutoClose fpowner(fp);

  // defaults
  std::string color = "silver";
  std::string material = "Transparent";

  int iarg = 1;
  std::string thisarg = arg[iarg];
  while (iarg < narg) {
    if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "region2vmd", error);
    thisarg = arg[iarg];
    ++iarg;

    if (thisarg == "color") {
      color = arg[iarg];
     if (const auto &search = vmdcolors.find(color); search == vmdcolors.end())
        error->all(FLERR, iarg, "Color {} is not a known VMD color", color);

    } else if (thisarg == "material") {
      material = arg[iarg];
      if (const auto &search = vmdmaterials.find(material); search == vmdmaterials.end())
        error->all(FLERR, iarg, "Material {} is not a known VMD material", material);

    } else if (thisarg == "command") {
      if (fp) {
        fputs("\n# custom command\n", fp);
        fputs(arg[iarg], fp);
        fputs("\n", fp);
      }

    } else if (thisarg == "region") {
      auto *region = domain->get_region_by_id(arg[iarg]);
      if (!region) {
        error->all(FLERR, iarg, "Region {} does not exist", arg[iarg]);
      } else {
        if (fp) {
          utils::logmesg(lmp, " writing region {} ...", region->id);
          utils::print(fp, "\n# region {} of style {}\n", region->id, region->style);

          fputs("# create new empty VMD molecule to store the graphics primitives\n"
                "set gfxmol [mol new]\nmol top $gfxmol\n",
                fp);
          utils::print(fp, "mol rename $gfxmol {{LAMMPS region {}}}\n", region->id);
          fputs("# set color and material\n", fp);
          utils::print(fp, "graphics $gfxmol color {}\n", color);
          utils::print(fp, "graphics $gfxmol material {}\n", material);
          write_region(fp, region);
        }
      }

    } else {
      error->all(FLERR, iarg - 1, "Unknown region2vmd keyword {}", thisarg);
    }
    ++iarg;
  }

  // done. file will be close automatically
  if (fp) {
    // reset views and restore previous top molecule
    fputs("after idle {if {$oldtop >= 0} {mol top $oldtop}; display resetview}\n", fp);
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
    utils::logmesg(lmp, "Cannot (yet) handle moving or rotating regions {}. Skipping... ",
                   region->id);
    return;
  }

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
        // a VMD cone uses a single cone primitive
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
        // a VMD cone uses a single cone primitive
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
