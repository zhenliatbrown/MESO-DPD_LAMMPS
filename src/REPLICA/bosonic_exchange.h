/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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
        BosonicExchange(class LAMMPS *, int nbosons, int np, int bead_num, bool mic, bool iPyConvention);
        ~BosonicExchange();

        void prepare_with_coordinates(const double* x, const double* x_prev, const double* x_next,
                                      double beta, double kT, double spring_constant);

        double get_potential() const;
        double get_total_spring_energy_for_bead();
        double get_Vn(int n) const;

        // CR: actually this function can also be const? (probably my bad)
        // OB: I don't think it can (see below)
        void spring_force(double** f);

        double prim_estimator();
        // OB: I left **x and **f not const such that they'll be compatible with atom->x and atom->f
        double vir_estimator(double **x, double **f) const;

    private:
        void evaluate_cycle_energies();
        void diff_two_beads(const double* x1, int l1, const double* x2, int l2, double diff[3]);
        double distance_squared_two_beads(const double* x1, int l1, const double* x2, int l2);
        double get_Enk(int m, int k);
        void set_Enk(int m, int k, double val);
        void evaluate_connection_probabilities();
        // CR: actually these force functions can also be const? (probably my bad)
        // OB: as written now they cannot because they use "diff_two_beads" - it gave me an error...
        void spring_force_last_bead(double** f);
        void spring_force_first_bead(double** f);
        void spring_force_interior_bead(double** f);
        void Evaluate_VBn();
        void Evaluate_V_backwards();
        void calc_total_spring_energy_for_bead();

        const int nbosons;
        const int np;
        const int bead_num;
        const bool apply_minimum_image;

        double spring_constant;
        double beta;
        const double* x;
        const double* x_prev;
        const double* x_next;

        double* E_kn;
        double* V;
        double* V_backwards;
        double* connection_probabilities;

        double* temp_nbosons_array;
        double spring_energy_for_bead;

        bool iPyConvention;
    };


}

#endif