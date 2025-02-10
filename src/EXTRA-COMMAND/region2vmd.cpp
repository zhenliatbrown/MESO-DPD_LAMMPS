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
#include "region_prism.h"
#include "region_sphere.h"

#include <cstring>
#include <unordered_set>

using namespace LAMMPS_NS;

static constexpr double SMALL = 1.0e-10;
static constexpr double DELTA = 1.0e-5;

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

  if (region->rotateflag) {
    utils::logmesg(lmp, "Cannot (yet) handle rotating region {}. Skipping... ", region->id);
    return;
  }

  // compute position offset for moving regions

  double dx = 0.0;
  double dy = 0.0;
  double dz = 0.0;
  if (region->moveflag) {
    dx = region->dx;
    dy = region->dy;
    dz = region->dz;
  }

  // translate compatible regions to VMD graphics primitives, skip others.

  const std::string regstyle = region->style;
  if (regstyle == "block") {
    const auto block = dynamic_cast<RegBlock *>(region);
    if (!block) {
      error->one(FLERR, Error::NOLASTLINE, "Region {} is not of style 'block'", region->id);
    } else {
      // a block is represented by 12 triangles
      // comment out VMD command when side is open
      utils::print(fp,
                   "{6}draw triangle {{{0} {2} {4}}} {{{0} {2} {5}}} {{{0} {3} {4}}}\n"
                   "{6}draw triangle {{{0} {3} {4}}} {{{0} {3} {5}}} {{{0} {2} {5}}}\n"
                   "{8}draw triangle {{{0} {2} {4}}} {{{0} {2} {5}}} {{{1} {2} {4}}}\n"
                   "{8}draw triangle {{{1} {2} {4}}} {{{1} {2} {5}}} {{{0} {2} {5}}}\n"
                   "{7}draw triangle {{{1} {2} {4}}} {{{1} {2} {5}}} {{{1} {3} {4}}}\n"
                   "{7}draw triangle {{{1} {3} {4}}} {{{1} {3} {5}}} {{{1} {2} {5}}}\n"
                   "{10}draw triangle {{{0} {2} {4}}} {{{1} {2} {4}}} {{{1} {3} {4}}}\n"
                   "{10}draw triangle {{{0} {2} {4}}} {{{0} {3} {4}}} {{{1} {3} {4}}}\n"
                   "{11}draw triangle {{{0} {2} {5}}} {{{1} {2} {5}}} {{{1} {3} {5}}}\n"
                   "{11}draw triangle {{{0} {2} {5}}} {{{0} {3} {5}}} {{{1} {3} {5}}}\n"
                   "{9}draw triangle {{{0} {3} {4}}} {{{0} {3} {5}}} {{{1} {3} {5}}}\n"
                   "{9}draw triangle {{{0} {3} {4}}} {{{1} {3} {4}}} {{{1} {3} {5}}}\n",
                   block->xlo + dx, block->xhi + dx, block->ylo + dy, block->yhi + dy,
                   block->zlo + dz, block->zhi + dz, block->open_faces[0] ? "# " : "",
                   block->open_faces[1] ? "# " : "", block->open_faces[2] ? "# " : "",
                   block->open_faces[3] ? "# " : "", block->open_faces[4] ? "# " : "",
                   block->open_faces[5] ? "# " : "");
    }

  } else if (regstyle == "cone") {
    const auto cone = dynamic_cast<RegCone *>(region);
    if (!cone) {
      error->one(FLERR, Error::NOLASTLINE, "Region {} is not of style 'cone'", region->id);
    } else {
      if (cone->open_faces[0] || cone->open_faces[1])
        error->warning(FLERR, "Drawing open-faced cones is not supported");
      // The VMD cone primitive requires one radius set to zero
      if (cone->radiuslo < SMALL) {
        // a VMD cone uses a single cone primitive
        if (cone->axis == 'x') {
          utils::print(fp, "draw cone {{{1} {2} {3}}} {{{0} {2} {3}}} radius {4} resolution 20\n",
                       cone->lo + dx, cone->hi + dx, cone->c1 + dy, cone->c2 + dz, cone->radiushi);
        } else if (cone->axis == 'y') {
          utils::print(fp, "draw cone {{{2} {1} {3}}} {{{2} {0} {3}}} radius {4} resolution 20\n",
                       cone->lo + dy, cone->hi + dy, cone->c1 + dx, cone->c2 + dz, cone->radiushi);
        } else if (cone->axis == 'z') {
          utils::print(fp, "draw cone {{{2} {3} {1}}} {{{2} {3} {0}}} radius {4} resolution 20\n",
                       cone->lo + dz, cone->hi + dz, cone->c1 + dx, cone->c2 + dy, cone->radiushi);
        }
      } else if (cone->radiushi < SMALL) {
        // a VMD cone uses a single cone primitive
        if (cone->axis == 'x') {
          utils::print(fp, "draw cone {{{0} {2} {3}}} {{{1} {2} {3}}} radius {4} resolution 20\n",
                       cone->lo + dx, cone->hi + dx, cone->c1 + dy, cone->c2 + dz, cone->radiuslo);
        } else if (cone->axis == 'y') {
          utils::print(fp, "draw cone {{{2} {0} {3}}} {{{2} {1} {3}}} radius {4} resolution 20\n",
                       cone->lo + dy, cone->hi + dy, cone->c1 + dx, cone->c2 + dz, cone->radiuslo);
        } else if (cone->axis == 'z') {
          utils::print(fp, "draw cone {{{2} {3} {0}}} {{{2} {3} {1}}} radius {4} resolution 20\n",
                       cone->lo + dz, cone->hi + dz, cone->c1 + dx, cone->c2 + dy, cone->radiuslo);
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
      std::string filled = "yes";
      if (cylinder->open_faces[0] && cylinder->open_faces[1]) {
        filled = "no";
      } else if (cylinder->open_faces[0] != cylinder->open_faces[1]) {
        filled = "no";
        // we put a single "lid" on an open cylinder by adding a filled cylinder of zero height
        double lid = cylinder->lo;
        if (cylinder->open_faces[0]) lid = cylinder->hi;
        if (cylinder->axis == 'x') {
          utils::print(fp,
                       "draw cylinder {{{0} {2} {3}}} {{{1:.15} {2} {3}}} radius {4} resolution 20 "
                       "filled yes\n",
                       lid + dx, lid + dx + DELTA, cylinder->c1 + dy, cylinder->c2 + dz,
                       cylinder->radius);
        } else if (cylinder->axis == 'y') {
          utils::print(fp,
                       "draw cylinder {{{2} {0} {3}}} {{{2} {1:.15} {3}}} radius {4} resolution 20 "
                       "filled yes\n",
                       lid + dy, lid + dy + DELTA, cylinder->c1 + dx, cylinder->c2 + dz,
                       cylinder->radius);
        } else if (cylinder->axis == 'z') {
          utils::print(fp,
                       "draw cylinder {{{2} {3} {0}}} {{{2} {3} {1:.15}}} radius {4} resolution 20 "
                       "filled yes\n",
                       lid + dz, lid + dz + DELTA, cylinder->c1 + dx, cylinder->c2 + dy,
                       cylinder->radius);
        }
      }
      if (cylinder->open_faces[2]) {
        // need to handle two lids case only. Single lid is already done
        if (!cylinder->open_faces[0] && !cylinder->open_faces[1]) {
          if (cylinder->axis == 'x') {
            utils::print(
                fp,
                "draw cylinder {{{0} {2} {3}}} {{{1:.15} {2} {3}}} radius {4} resolution 20 "
                "filled yes\n",
                cylinder->lo + dx, cylinder->lo + dx + DELTA, cylinder->c1 + dy, cylinder->c2 + dz,
                cylinder->radius);
            utils::print(
                fp,
                "draw cylinder {{{0} {2} {3}}} {{{1:.15} {2} {3}}} radius {4} resolution 20 "
                "filled yes\n",
                cylinder->hi + dx, cylinder->hi + dx + DELTA, cylinder->c1 + dy, cylinder->c2 + dz,
                cylinder->radius);
          } else if (cylinder->axis == 'y') {
            utils::print(
                fp,
                "draw cylinder {{{2} {0} {3}}} {{{2} {1:.15} {3}}} radius {4} resolution 20 "
                "filled yes\n",
                cylinder->lo + dy, cylinder->lo + dy + DELTA, cylinder->c1 + dx, cylinder->c2 + dz,
                cylinder->radius);
            utils::print(
                fp,
                "draw cylinder {{{2} {0} {3}}} {{{2} {1:.15} {3}}} radius {4} resolution 20 "
                "filled yes\n",
                cylinder->hi + dy, cylinder->hi + dy + DELTA, cylinder->c1 + dx, cylinder->c2 + dz,
                cylinder->radius);
          } else if (cylinder->axis == 'z') {
            utils::print(
                fp,
                "draw cylinder {{{2} {3} {0}}} {{{2} {3} {1:.15}}} radius {4} resolution 20 "
                "filled yes\n",
                cylinder->lo + dz, cylinder->lo + dz + DELTA, cylinder->c1 + dx, cylinder->c2 + dy,
                cylinder->radius);
            utils::print(
                fp,
                "draw cylinder {{{2} {3} {0}}} {{{2} {3} {1:.15}}} radius {4} resolution 20 "
                "filled yes\n",
                cylinder->hi + dz, cylinder->hi + dz + DELTA, cylinder->c1 + dx, cylinder->c2 + dy,
                cylinder->radius);
          }
        }
      } else {
        // a cylinder uses a single cylinder primitive and possibly a single "lid"
        if (cylinder->axis == 'x') {
          utils::print(
              fp,
              "draw cylinder {{{0} {2} {3}}} {{{1} {2} {3}}} radius {4} resolution 20 filled {5}\n",
              cylinder->lo + dx, cylinder->hi + dx, cylinder->c1 + dy, cylinder->c2 + dz,
              cylinder->radius, filled);
        } else if (cylinder->axis == 'y') {
          utils::print(
              fp,
              "draw cylinder {{{2} {0} {3}}} {{{2} {1} {3}}} radius {4} resolution 20 filled {5}\n",
              cylinder->lo + dy, cylinder->hi + dy, cylinder->c1 + dx, cylinder->c2 + dz,
              cylinder->radius, filled);
        } else if (cylinder->axis == 'z') {
          utils::print(
              fp,
              "draw cylinder {{{2} {3} {0}}} {{{2} {3} {1}}} radius {4} resolution 20 filled {5}\n",
              cylinder->lo + dz, cylinder->hi + dz, cylinder->c1 + dx, cylinder->c2 + dy,
              cylinder->radius, filled);
        }
      }
    }

  } else if (regstyle == "prism") {
    const auto prism = dynamic_cast<RegPrism *>(region);
    if (!prism) {
      error->one(FLERR, Error::NOLASTLINE, "Region {} is not of style 'prism'", region->id);
    } else {
      // a prism is represented by 12 triangles
      // comment out VMD command when side is open
      utils::print(fp,
                   "{14}draw triangle {{{0} {2} {4}}} {{{1} {2} {4}}} {{{7} {3} {4}}}\n"
                   "{14}draw triangle {{{0} {2} {4}}} {{{6} {3} {4}}} {{{7} {3} {4}}}\n"
                   "{16}draw triangle {{{0} {2} {4}}} {{{1} {2} {4}}} {{{8} {9} {5}}}\n"
                   "{16}draw triangle {{{0} {2} {4}}} {{{10} {9} {5}}} {{{8} {9} {5}}}\n"
                   "{19}draw triangle {{{1} {2} {4}}} {{{8} {9} {5}}} {{{7} {3} {4}}}\n"
                   "{19}draw triangle {{{11} {12} {5}}} {{{8} {9} {5}}} {{{7} {3} {4}}}\n"
                   "{18}draw triangle {{{0} {2} {4}}} {{{6} {3} {4}}} {{{13} {12} {5}}}\n"
                   "{18}draw triangle {{{0} {2} {4}}} {{{10} {9} {5}}} {{{13} {12} {5}}}\n"
                   "{15}draw triangle {{{10} {9} {5}}} {{{8} {9} {5}}} {{{11} {12} {5}}}\n"
                   "{15}draw triangle {{{10} {9} {5}}} {{{13} {12} {5}}} {{{11} {12} {5}}}\n"
                   "{17}draw triangle {{{6} {3} {4}}} {{{7} {3} {4}}} {{{11} {12} {5}}}\n"
                   "{17}draw triangle {{{6} {3} {4}}} {{{13} {12} {5}}} {{{11} {12} {5}}}\n",
                   prism->xlo, prism->xhi, prism->ylo, prism->yhi, prism->zlo, prism->zhi,
                   prism->xlo + prism->xy, prism->xhi + prism->xy, prism->xhi + prism->xz,
                   prism->ylo + prism->yz, prism->xlo + prism->xz,
                   prism->xhi + prism->xy + prism->xz, prism->yhi + prism->yz,
                   prism->xlo + prism->xy + prism->xz, prism->open_faces[0] ? "# " : "",
                   prism->open_faces[1] ? "# " : "", prism->open_faces[2] ? "# " : "",
                   prism->open_faces[3] ? "# " : "", prism->open_faces[4] ? "# " : "",
                   prism->open_faces[5] ? "# " : "");
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
