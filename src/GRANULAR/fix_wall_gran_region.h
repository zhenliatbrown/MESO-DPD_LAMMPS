/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(wall/gran/region,FixWallGranRegion);
// clang-format on
#else

#ifndef LMP_FIX_WALL_GRAN_REGION_H
#define LMP_FIX_WALL_GRAN_REGION_H

#include "fix_wall_gran.h"

namespace LAMMPS_NS {

class FixWallGranRegion : public FixWallGran {
 public:
  FixWallGranRegion(class LAMMPS *, int, char **);
  ~FixWallGranRegion() override;
  void post_force(int) override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  void init() override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;

  class Region *region; // MOVED FROM PRIVATE TO PUBLIC FOR MDR MODEL
  void update_contacts(int, int); // MOVED FROM PRIVATE TO PUBLIC FOR MDR MODEL
  int *ncontact;             // MOVED FROM PRIVATE TO PUBLIC FOR MDR MODEL
  int **walls;               // MOVED FROM PRIVATE TO PUBLIC FOR MDR MODEL
  int *c2r;                  // MOVED FROM PRIVATE TO PUBLIC FOR MDR MODEL
  double ***history_many;    // MOVED FROM PRIVATE TO PUBLIC FOR MDR MODEL
  int tmax;                  // MOVED FROM PRIVATE TO PUBLIC FOR MDR MODEL

 private:

  int nregion;

  // shear history for multiple contacts per particle




                             // c2r[i] = index of Ith contact in
                             //   region-contact[] list of contacts
  int motion_resetflag;      // used by restart to indicate that region
                             //    vel info is to be reset


};

}    // namespace LAMMPS_NS

#endif
#endif
