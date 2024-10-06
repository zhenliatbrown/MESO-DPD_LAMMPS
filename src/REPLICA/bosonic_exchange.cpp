#include "bosonic_exchange.h"

#include "domain.h"
#include "universe.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

BosonicExchange::BosonicExchange(LAMMPS *lmp, int nbosons, int np, int bead_num, bool mic) :
        Pointers(lmp),
        nbosons(nbosons), np(np), bead_num(bead_num), apply_minimum_image(mic) {
    memory->create(temp_nbosons_array, nbosons, "BosonicExchange: temp_nbosons_array");
    memory->create(separate_atom_spring, nbosons, "BosonicExchange: separate_atom_spring");
    memory->create(E_kn, (nbosons * (nbosons + 1) / 2), "BosonicExchange: E_kn");
    memory->create(V, nbosons + 1, "BosonicExchange: V");
    memory->create(V_backwards, nbosons + 1, "BosonicExchange: V_backwards");
    memory->create(connection_probabilities, nbosons * nbosons, "BosonicExchange: connection probabilities");
    memory->create(prim_est, nbosons + 1, "BosonicExchange: prim_est");
}

void BosonicExchange::prepare_with_coordinates(const double* x, const double* x_prev, const double* x_next,
                                               double beta, double kT, double spring_constant) {
    this->x = x;
    this->x_prev = x_prev;
    this->x_next = x_next;
    this->beta = beta;
    this->kT = kT;
    this->spring_constant = spring_constant;

    evaluate_cycle_energies();
    if (bead_num == 0 || bead_num == np - 1) {
        // exterior beads
        Evaluate_VBn();
        Evaluate_V_backwards();
        evaluate_connection_probabilities();
    }
}

/* ---------------------------------------------------------------------- */

BosonicExchange::~BosonicExchange() {
    memory->destroy(prim_est);
    memory->destroy(connection_probabilities);
    memory->destroy(V_backwards);
    memory->destroy(V);
    memory->destroy(E_kn);
    memory->destroy(separate_atom_spring);
    memory->destroy(temp_nbosons_array);
}

/* ---------------------------------------------------------------------- */

void BosonicExchange::diff_two_beads(const double* x1, int l1, const double* x2, int l2,
                              double diff[3]) {
    l1 = l1 % nbosons;
    l2 = l2 % nbosons;
    double delx2 = x2[3 * l2 + 0] - x1[3 * l1 + 0];
    double dely2 = x2[3 * l2 + 1] - x1[3 * l1 + 1];
    double delz2 = x2[3 * l2 + 2] - x1[3 * l1 + 2];
    if (apply_minimum_image) {
        domain->minimum_image(delx2, dely2, delz2);
    }

    diff[0] = delx2;
    diff[1] = dely2;
    diff[2] = delz2;
}

/* ---------------------------------------------------------------------- */

double BosonicExchange::distance_squared_two_beads(const double* x1, int l1, const double* x2, int l2) {
    double diff[3];
    diff_two_beads(x1, l1, x2, l2, diff);
    return diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
}

/* ---------------------------------------------------------------------- */

void BosonicExchange::evaluate_cycle_energies() // VVVVV
{ 
    for (int i = 0; i < nbosons; i++) {
        temp_nbosons_array[i] = distance_squared_two_beads(x, i, x_next, i);
    }
    // Reduce the result and send to bead_num=0
    MPI_Reduce(temp_nbosons_array, separate_atom_spring, nbosons,
                  MPI_DOUBLE, MPI_SUM, 0, universe->uworld);

    if (bead_num == 0 || bead_num == np - 1) {
        const double* x_first_bead;
        const double* x_last_bead;
        if (bead_num == 0) {
            // Send to bead_num=np-1
            MPI_Send(separate_atom_spring, nbosons, MPI_DOUBLE, np - 1, 0, universe->uworld);

            x_first_bead = x;
            x_last_bead = x_prev;
        } else {
            // Receive at bead_num=np-1 from bead_num=0
            MPI_Recv(separate_atom_spring, nbosons, MPI_DOUBLE, 0, 0, universe->uworld, MPI_STATUS_IGNORE);
            
            x_first_bead = x_next;
            x_last_bead = x;
        }

        for (int v = 0; v < nbosons; v++) {
            set_Enk(v + 1, 1,
                    0.5 * spring_constant * (separate_atom_spring[v]));

            for (int u = v - 1; u >= 0; u--) {
                double val = get_Enk(v + 1, v - u) +
                             0.5 * spring_constant * (
                                     // Eint(u)
                                     separate_atom_spring[u] - distance_squared_two_beads(x_first_bead, u, x_last_bead, u)
                                     // connect u to u+1
                                     + distance_squared_two_beads(x_last_bead, u, x_first_bead, u + 1)
                                     // break cycle [u+1,v]
                                     - distance_squared_two_beads(x_first_bead, u + 1, x_last_bead, v)
                                     // close cycle from v to u
                                     + distance_squared_two_beads(x_first_bead, u, x_last_bead, v));

                set_Enk(v + 1, v - u + 1, val);
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

double BosonicExchange::get_Enk(int m, int k) {
    int end_of_m = m * (m + 1) / 2;
    return E_kn[end_of_m - k];
}

/* ---------------------------------------------------------------------- */

void BosonicExchange::set_Enk(int m, int k, double val) {
    int end_of_m = m * (m + 1) / 2;
    E_kn[end_of_m - k] = val;
}

/* ---------------------------------------------------------------------- */

void BosonicExchange::Evaluate_VBn()
{
    V[0] = 0.0;

    for (int m = 1; m < nbosons + 1; m++) {
        double Elongest = std::numeric_limits<double>::max();

        for (int k = m; k > 0; k--) {
            double val = get_Enk(m,k) + V[m-k];
            Elongest = std::min(Elongest, val);
            temp_nbosons_array[k - 1] = val;
        }

        double sig_denom = 0.0;
        for (int k = m; k > 0; k--) {
            sig_denom += exp(-beta * (temp_nbosons_array[k - 1] - Elongest));
        }
        V[m] = Elongest - (1.0 / beta) * log(sig_denom / (double)m);

        if (!std::isfinite(V[m])) {
            error->universe_one(
                    FLERR,
                    fmt::format("Invalid sig_denom {} with Elongest {} in bosonic exchange potential",
                                sig_denom, Elongest));
        }
    }
}

void BosonicExchange::Evaluate_V_backwards() {
    V_backwards[nbosons] = 0.0;

    for (int l = nbosons - 1; l > 0; l--) {
        double Elongest = std::numeric_limits<double>::max();
        for (int p = l; p < nbosons; p++) {
            double val = get_Enk(p + 1, p - l + 1) + V_backwards[p + 1];
            Elongest = std::min(Elongest, val);
            temp_nbosons_array[p] = val;
        }

        double sig_denom = 0.0;
        for (int p = l; p < nbosons; p++) {
            sig_denom += 1.0 / (p + 1) * exp(-beta *
                                             (temp_nbosons_array[p]
                                              - Elongest)
            );
        }

        V_backwards[l] = Elongest - log(sig_denom) / beta;

        if (!std::isfinite(V_backwards[l])) {
            error->universe_one(
                    FLERR,
                    fmt::format("Invalid sig_denom {} with Elongest {} in bosonic exchange potential backwards",
                                sig_denom, Elongest));
        }
    }

    V_backwards[0] = V[nbosons];
}


/* ---------------------------------------------------------------------- */

double BosonicExchange::get_potential() const {
    return V[nbosons];
}

/* ---------------------------------------------------------------------- */

double BosonicExchange::get_Vn(int n) const {
    return V[n];
}

/* ---------------------------------------------------------------------- */

double BosonicExchange::get_E_kn_serial_order(int i) const {
    return E_kn[i];
}

/* ---------------------------------------------------------------------- */

void BosonicExchange::spring_force(double** f) {
    if (bead_num == np - 1) {
        spring_force_last_bead(f);
    } else if (bead_num == 0) {
        spring_force_first_bead(f);
    } else {
        spring_force_interior_bead(f);
    }
}

/* ---------------------------------------------------------------------- */

void BosonicExchange::evaluate_connection_probabilities() {
    for (int l = 0; l < nbosons - 1; l++) {
        double direct_link_probability = 1.0 - (exp(-beta *
                                                    (V[l + 1] + V_backwards[l + 1] -
                                                     V[nbosons])));
        connection_probabilities[nbosons * l + (l + 1)] = direct_link_probability;
    }
    for (int u = 0; u < nbosons; u++) {
        for (int l = u; l < nbosons; l++) {
            double close_cycle_probability = 1.0 / (l + 1) *
                                             exp(-beta * (V[u] + get_Enk(l + 1, l - u + 1) + V_backwards[l + 1]
                                                          - V[nbosons]));
            connection_probabilities[nbosons * l + u] = close_cycle_probability;
        }
    }
}

/* ---------------------------------------------------------------------- */

void BosonicExchange::spring_force_last_bead(double** f)
{
    const double* x_first_bead = x_next;
    const double* x_last_bead = x;

    for (int l = 0; l < nbosons; l++) {
        double sum_x = 0.0;
        double sum_y = 0.0;
        double sum_z = 0.0;
        for (int next_l = 0; next_l <= l + 1 && next_l < nbosons; next_l++) {
            double diff_next[3];

            diff_two_beads(x_last_bead, l, x_first_bead, next_l, diff_next);

            double prob = connection_probabilities[nbosons * l + next_l];

            sum_x += prob * diff_next[0];
            sum_y += prob * diff_next[1];
            sum_z += prob * diff_next[2];
        }

        double diff_prev[3];
        diff_two_beads(x_last_bead, l, x_prev, l, diff_prev);
        sum_x += diff_prev[0];
        sum_y += diff_prev[1];
        sum_z += diff_prev[2];

        f[l][0] += sum_x * spring_constant;
        f[l][1] += sum_y * spring_constant;
        f[l][2] += sum_z * spring_constant;
    }
}

/* ---------------------------------------------------------------------- */

void BosonicExchange::spring_force_first_bead(double** f)
{
    const double* x_first_bead = x;
    const double* x_last_bead = x_prev;

    for (int l = 0; l < nbosons; l++) {
        double sum_x = 0.0;
        double sum_y = 0.0;
        double sum_z = 0.0;
        for (int prev_l = std::max(0, l - 1); prev_l < nbosons; prev_l++) {
            double diff_prev[3];

            diff_two_beads(x_first_bead, l, x_last_bead, prev_l, diff_prev);

            double prob = connection_probabilities[nbosons * prev_l + l];

            sum_x += prob * diff_prev[0];
            sum_y += prob * diff_prev[1];
            sum_z += prob * diff_prev[2];
        }

        double diff_next[3];
        diff_two_beads(x_first_bead, l, x_next, l, diff_next);
        sum_x += diff_next[0];
        sum_y += diff_next[1];
        sum_z += diff_next[2];

        f[l][0] += sum_x * spring_constant;
        f[l][1] += sum_y * spring_constant;
        f[l][2] += sum_z * spring_constant;
    }
}

/* ---------------------------------------------------------------------- */

void BosonicExchange::spring_force_interior_bead(double **f) // VVVVVV
{
    for (int l = 0; l < nbosons; l++) {
        double sum_x = 0.0;
        double sum_y = 0.0;
        double sum_z = 0.0;

        double diff_prev[3];
        diff_two_beads(x, l, x_prev, l, diff_prev);
        sum_x += diff_prev[0];
        sum_y += diff_prev[1];
        sum_z += diff_prev[2];

        double diff_next[3];
        diff_two_beads(x, l, x_next, l, diff_next);
        sum_x += diff_next[0];
        sum_y += diff_next[1];
        sum_z += diff_next[2];

        f[l][0] += sum_x * spring_constant;
        f[l][1] += sum_y * spring_constant;
        f[l][2] += sum_z * spring_constant;
    }
}

/* ---------------------------------------------------------------------- */

double BosonicExchange::prim_estimator()
{
  prim_est[0] = 0.0;

  for (int m = 1; m < nbosons + 1; ++m) {
    double sig = 0.0;

    // Numerical stability (Xiong & Xiong method)
    double Elongest = std::numeric_limits<double>::max();

    for (int k = m; k > 0; k--) {
      Elongest = std::min(Elongest, get_Enk(m, k) + V[m - k]);
    }
    
    for (int k = m; k > 0; --k) {
      double E_kn_val = get_Enk(m, k);

      sig += (prim_est[m - k] - E_kn_val) * exp(-beta * (E_kn_val + V[m - k] - Elongest));
    }

    double sig_denom_m = m * exp(-beta * (V[m] - Elongest));

    prim_est[m] = sig / sig_denom_m;
  }

  return 0.5 * domain->dimension * nbosons * np * kT + prim_est[nbosons];
}

/* ---------------------------------------------------------------------- */

double BosonicExchange::vir_estimator(double **x, double **f)
{
  double virial = 0;
  for (int i = 0; i < nbosons; i++) {
      virial += -0.5 * (x[i][0] * f[i][0] + x[i][1] * f[i][1] + x[i][2] * f[i][2]);
  }
  return virial;
}