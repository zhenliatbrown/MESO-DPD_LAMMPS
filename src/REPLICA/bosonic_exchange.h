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

#ifndef BOSONIC_EXCHANGE_H
#define BOSONIC_EXCHANGE_H

#include "pointers.h"

namespace LAMMPS_NS {

class BosonicExchange : protected Pointers {
 public:
  BosonicExchange(class LAMMPS *, int nbosons, int np, int bead_num, bool mic,
                  bool beta_convention);
  ~BosonicExchange();

  void prepare_with_coordinates(const double *x, const double *x_prev, const double *x_next,
                                double beta, double spring_constant);

  double get_potential() const;
  double get_bead_spring_energy() const;

  void spring_force(double **f) const;

  double prim_estimator();

 private:
  void evaluate_cycle_energies();
  void diff_two_beads(const double *x1, int l1, const double *x2, int l2, double diff[3]) const;
  double get_interior_bead_spring_energy() const;
  double distance_squared_two_beads(const double *x1, int l1, const double *x2, int l2) const;
  double get_Enk(int m, int k) const;
  void set_Enk(int m, int k, double val);
  void evaluate_connection_probabilities();
  void spring_force_last_bead(double **f) const;
  void spring_force_first_bead(double **f) const;
  void spring_force_interior_bead(double **f) const;
  void Evaluate_VBn();
  void Evaluate_V_backwards();

  const int nbosons;
  const int np;
  const int bead_num;
  const bool apply_minimum_image;
  // In the "reduced-beta convention" [e.g. in J. Chem. Phys. 133, 124104 (2010); also J. Chem. Phys. 74, 4078-4095 (1981)],
  // the Boltzmann exponents have the form exp[-(beta/P)H], where H is the classical Hamiltonian of the
  // ring polymers. This results in a canonical distribution at P times the physical temperature.
  // In contrast, the "physical-beta convention" [e.g. in J. Chem. Phys. 99, 2796-2808 (1993)] uses weights of the form exp(-beta*H),
  // such that the temperature of the canonical ensemble coincides with the physical temperature.
  // Notably, the classical Hamiltonians of the two conventions differ, with the spring constant
  // in the reduced-beta convention being P times larger than that in the physical-beta convention. Additionally, the reduced-beta convention
  // lacks a 1/P prefactor in front of the external potential. The Hamiltonians of the two conventions are related through
  // H_physical = H_reduced / P. Note however that the expressions for the various estimators are unaffected by this choice,
  // so as the algorithm for bosonic exchange. The code below was designed to be compatible with both conventions,
  // and the choice of convention only affects a single calculation within it.
  // Setting the following boolian variable to false amounts to adopting the reduced-beta convention.
  const bool physical_beta_convention;

  double spring_constant;
  double beta;
  const double *x;
  const double *x_prev;
  const double *x_next;

  double *E_kn;
  double *V;
  double *V_backwards;
  double *connection_probabilities;

  double *temp_nbosons_array;
};
}    // namespace LAMMPS_NS
#endif
