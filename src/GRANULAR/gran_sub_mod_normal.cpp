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

#include "gran_sub_mod_normal.h"

#include "atom.h"
#include "error.h"
#include "citeme.h"
#include "fix_granular_mdr.h"
#include "granular_model.h"
#include "math_const.h"
#include "modify.h"
#include "update.h"

#include <cmath>
#include <iomanip>
#include <sstream>

using namespace LAMMPS_NS;
using namespace Granular_NS;
using namespace MathConst;

static constexpr double PI27SQ = 266.47931882941264802866;      // 27*PI**2
static constexpr double THREEROOT3 = 5.19615242270663202362;    // 3*sqrt(3)
static constexpr double SIXROOT6 = 14.69693845669906728801;     // 6*sqrt(6)
static constexpr double INVROOT6 = 0.40824829046386307274;      // 1/sqrt(6)
static constexpr double FOURTHIRDS = (4.0 / 3.0);               // 4/3
static constexpr double FIVETHIRDS = (5.0 / 3.0);               // 5/3
static constexpr double JKRPREFIX = 1.2277228507842888;         // cbrt(3*PI**2/16)
static constexpr int MDR_MAX_IT = 100;                          // Newton-Raphson for MDR
static constexpr double MDR_EPSILON1 = 1e-10;                   // Newton-Raphson for MDR
static constexpr double MDR_EPSILON2 = 1e-16;                   // Newton-Raphson for MDR

static const char cite_mdr[] =
    "MDR contact model command: \n\n"
    "@Article\n"
    " author = \n"
    " title =   \n"
    " journal = ,\n"
    " year =    ,\n"
    " volume =  ,\n"
    " pages =  \n"
    "}\n\n";

/* ----------------------------------------------------------------------
   Default normal model
------------------------------------------------------------------------- */

GranSubModNormal::GranSubModNormal(GranularModel *gm, LAMMPS *lmp) : GranSubMod(gm, lmp)
{
  material_properties = 0;
  cohesive_flag = 0;
}

/* ---------------------------------------------------------------------- */

bool GranSubModNormal::touch()
{
  bool touchflag = (gm->rsq < gm->radsum * gm->radsum);
  return touchflag;
}

/* ---------------------------------------------------------------------- */

double GranSubModNormal::pulloff_distance(double /*radi*/, double /*radj*/)
{
  // called outside of compute(), do not assume correct geometry defined in contact
  return 0.0;
}

/* ---------------------------------------------------------------------- */

double GranSubModNormal::calculate_contact_radius()
{
  return sqrt(gm->dR);
}

/* ---------------------------------------------------------------------- */

void GranSubModNormal::set_fncrit()
{
  Fncrit = fabs(gm->Fntot);
}

/* ----------------------------------------------------------------------
   No model
------------------------------------------------------------------------- */

GranSubModNormalNone::GranSubModNormalNone(GranularModel *gm, LAMMPS *lmp) :
    GranSubModNormal(gm, lmp)
{
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalNone::calculate_forces()
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   Hookean normal force
------------------------------------------------------------------------- */

GranSubModNormalHooke::GranSubModNormalHooke(GranularModel *gm, LAMMPS *lmp) :
    GranSubModNormal(gm, lmp)
{
  num_coeffs = 2;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalHooke::coeffs_to_local()
{
  k = coeffs[0];
  damp = coeffs[1];

  if (k < 0.0 || damp < 0.0) error->all(FLERR, "Illegal Hooke normal model");
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalHooke::calculate_forces()
{
  return k * gm->delta;
}

/* ----------------------------------------------------------------------
   Hertzian normal force
------------------------------------------------------------------------- */

GranSubModNormalHertz::GranSubModNormalHertz(GranularModel *gm, LAMMPS *lmp) :
    GranSubModNormal(gm, lmp)
{
  num_coeffs = 2;
  contact_radius_flag = 1;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalHertz::coeffs_to_local()
{
  k = coeffs[0];
  damp = coeffs[1];

  if (k < 0.0 || damp < 0.0) error->all(FLERR, "Illegal Hertz normal model");
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalHertz::calculate_forces()
{
  return k * gm->contact_radius * gm->delta;
}

/* ----------------------------------------------------------------------
   Hertzian normal force with material properties
------------------------------------------------------------------------- */

GranSubModNormalHertzMaterial::GranSubModNormalHertzMaterial(GranularModel *gm, LAMMPS *lmp) :
    GranSubModNormalHertz(gm, lmp)
{
  material_properties = 1;
  num_coeffs = 3;
  contact_radius_flag = 1;
  mixed_coefficients = 0;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalHertzMaterial::coeffs_to_local()
{
  Emod = coeffs[0];
  damp = coeffs[1];
  poiss = coeffs[2];
  if (!mixed_coefficients) {
    if (gm->contact_type == PAIR) {
      k = FOURTHIRDS * mix_stiffnessE(Emod, Emod, poiss, poiss);
    } else {
      k = FOURTHIRDS * mix_stiffnessE_wall(Emod, poiss);
    }
  }

  if (Emod < 0.0 || damp < 0.0) error->all(FLERR, "Illegal Hertz material normal model");
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalHertzMaterial::mix_coeffs(double *icoeffs, double *jcoeffs)
{
  coeffs[0] = mix_stiffnessE(icoeffs[0], jcoeffs[0], icoeffs[2], jcoeffs[2]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);

  k = FOURTHIRDS * coeffs[0];
  mixed_coefficients = 1;

  coeffs_to_local();
}

/* ----------------------------------------------------------------------
   DMT normal force
------------------------------------------------------------------------- */

GranSubModNormalDMT::GranSubModNormalDMT(GranularModel *gm, LAMMPS *lmp) : GranSubModNormal(gm, lmp)
{
  material_properties = 1;
  cohesive_flag = 1;
  num_coeffs = 4;
  contact_radius_flag = 1;
  mixed_coefficients = 0;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalDMT::coeffs_to_local()
{
  Emod = coeffs[0];
  damp = coeffs[1];
  poiss = coeffs[2];
  cohesion = coeffs[3];

  if (!mixed_coefficients) {
    if (gm->contact_type == PAIR) {
      k = FOURTHIRDS * mix_stiffnessE(Emod, Emod, poiss, poiss);
    } else {
      k = FOURTHIRDS * mix_stiffnessE_wall(Emod, poiss);
    }
  }

  if (Emod < 0.0 || damp < 0.0) error->all(FLERR, "Illegal DMT normal model");
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalDMT::mix_coeffs(double *icoeffs, double *jcoeffs)
{
  coeffs[0] = mix_stiffnessE(icoeffs[0], jcoeffs[0], icoeffs[2], jcoeffs[2]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs[3] = mix_geom(icoeffs[3], jcoeffs[3]);

  k = FOURTHIRDS * coeffs[0];
  mixed_coefficients = 1;

  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalDMT::calculate_forces()
{
  Fne = k * gm->contact_radius * gm->delta;
  F_pulloff = 4.0 * MY_PI * cohesion * gm->Reff;
  Fne -= F_pulloff;
  return Fne;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalDMT::set_fncrit()
{
  Fncrit = fabs(Fne + 2.0 * F_pulloff);
}

/* ----------------------------------------------------------------------
   JKR normal force
------------------------------------------------------------------------- */

GranSubModNormalJKR::GranSubModNormalJKR(GranularModel *gm, LAMMPS *lmp) : GranSubModNormal(gm, lmp)
{
  material_properties = 1;
  cohesive_flag = 1;
  beyond_contact = 1;
  num_coeffs = 4;
  contact_radius_flag = 1;
  mixed_coefficients = 0;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalJKR::coeffs_to_local()
{
  Emod = coeffs[0];
  damp = coeffs[1];
  poiss = coeffs[2];
  cohesion = coeffs[3];

  if (!mixed_coefficients) {
    if (gm->contact_type == PAIR) {
      Emix = mix_stiffnessE(Emod, Emod, poiss, poiss);
    } else {
      Emix = mix_stiffnessE_wall(Emod, poiss);
    }
  }

  k = FOURTHIRDS * Emix;

  if (Emod < 0.0 || damp < 0.0) error->all(FLERR, "Illegal JKR normal model");
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalJKR::mix_coeffs(double *icoeffs, double *jcoeffs)
{
  coeffs[0] = mix_stiffnessE(icoeffs[0], jcoeffs[0], icoeffs[2], jcoeffs[2]);
  coeffs[1] = mix_geom(icoeffs[1], jcoeffs[1]);
  coeffs[2] = mix_geom(icoeffs[2], jcoeffs[2]);
  coeffs[3] = mix_geom(icoeffs[3], jcoeffs[3]);

  Emix = coeffs[0];
  mixed_coefficients = 1;

  coeffs_to_local();
}

/* ---------------------------------------------------------------------- */

bool GranSubModNormalJKR::touch()
{
  double delta_pulloff, dist_pulloff;
  bool touchflag;

  if (gm->touch) {
    // delta_pulloff defined as positive so center-to-center separation is > radsum
    delta_pulloff = JKRPREFIX * cbrt(gm->Reff * cohesion * cohesion / (Emix * Emix));
    dist_pulloff = gm->radsum + delta_pulloff;
    touchflag = gm->rsq < (dist_pulloff * dist_pulloff);
  } else {
    touchflag = gm->rsq < (gm->radsum * gm->radsum);
  }

  return touchflag;
}

/* ----------------------------------------------------------------------
  called outside of compute(), do not assume geometry defined in contact
------------------------------------------------------------------------- */

double GranSubModNormalJKR::pulloff_distance(double radi, double radj)
{
  double Reff_tmp;

  Reff_tmp = radi * radj / (radi + radj);    // May not be defined
  if (Reff_tmp <= 0) return 0;
  // Defined as positive so center-to-center separation is > radsum
  return JKRPREFIX * cbrt(Reff_tmp * cohesion * cohesion / (Emix * Emix));
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalJKR::calculate_contact_radius()
{
  double R2, dR2, t0, t1, t2, t3, t4, t5, t6;
  double sqrt1, sqrt2, sqrt3;

  R2 = gm->Reff * gm->Reff;
  dR2 = gm->dR * gm->dR;
  t0 = cohesion * cohesion * R2 * R2 * Emix;
  t1 = PI27SQ * t0;
  t2 = 8.0 * gm->dR * dR2 * Emix * Emix * Emix;
  t3 = 4.0 * dR2 * Emix;

  // in case sqrt(0) < 0 due to precision issues
  sqrt1 = MAX(0, t0 * (t1 + 2.0 * t2));
  t4 = cbrt(t1 + t2 + THREEROOT3 * MY_PI * sqrt(sqrt1));
  t5 = t3 / t4 + t4 / Emix;
  sqrt2 = MAX(0, 2.0 * gm->dR + t5);
  t6 = sqrt(sqrt2);
  sqrt3 = MAX(0, 4.0 * gm->dR - t5 + SIXROOT6 * cohesion * MY_PI * R2 / (Emix * t6));

  return INVROOT6 * (t6 + sqrt(sqrt3));
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalJKR::calculate_forces()
{
  double a2;
  a2 = gm->contact_radius * gm->contact_radius;
  Fne = k * gm->contact_radius * a2 / gm->Reff -
      MY_2PI * a2 * sqrt(4.0 * cohesion * Emix / (MY_PI * gm->contact_radius));
  F_pulloff = 3.0 * MY_PI * cohesion * gm->Reff;

  return Fne;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalJKR::set_fncrit()
{
  Fncrit = fabs(Fne + 2.0 * F_pulloff);
}

/* ----------------------------------------------------------------------
   MDR contact model

   Contributing authors:
   William Zunker (MIT), Sachith Dunatunga (MIT),
   Dan Bolintineanu (SNL), Joel Clemmer (SNL)
------------------------------------------------------------------------- */

GranSubModNormalMDR::GranSubModNormalMDR(GranularModel *gm, LAMMPS *lmp) :
    GranSubModNormal(gm, lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_mdr);

  num_coeffs = 6;
  contact_radius_flag = 1;
  size_history = 26;
  fix_mdr_flag = 0;
  id_fix = nullptr;

  nondefault_history_transfer = 1;
  transfer_history_factor = new double[size_history];
  for (int i = 0; i < size_history; i++) {
    transfer_history_factor[i] = +1;
  }
}

/* ---------------------------------------------------------------------- */

GranSubModNormalMDR::~GranSubModNormalMDR()
{
  if (id_fix && modify->nfix)
    modify->delete_fix(id_fix);
  delete[] id_fix;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalMDR::coeffs_to_local()
{
  E = coeffs[0];      // Young's modulus
  nu = coeffs[1];     // Poisson's ratio
  Y = coeffs[2];      // yield stress
  gamma = coeffs[3];  // effective surface energy
  psi_b = coeffs[4];  // bulk response trigger based on ratio of remaining free area: A_{free}/A_{total}
  CoR = coeffs[5];    // coefficent of restitution

  if (E <= 0.0) error->all(FLERR, "Illegal MDR normal model, Young's modulus must be greater than 0");
  if (nu < 0.0 || nu > 0.5) error->all(FLERR, "Illegal MDR normal model, Poisson's ratio must be between 0 and 0.5");
  if (Y < 0.0) error->all(FLERR, "Illegal MDR normal model, yield stress must be greater than or equal to 0");
  if (gamma < 0.0) error->all(FLERR, "Illegal MDR normal model, effective surface energy must be greater than or equal to 0");
  if (psi_b < 0.0 || psi_b > 1.0) error->all(FLERR, "Illegal MDR normal model, psi_b must be between 0 and 1.0");
  if (CoR < 0.0 || CoR > 1.0) error->all(FLERR, "Illegal MDR normal model, coefficent of restitution must be between 0 and 1.0");

  G = E / (2.0 * (1.0 + nu));            // shear modulus
  kappa = E / (3.0 * (1.0 - 2.0 * nu));  // bulk modulus
  Eeff = E / (1.0 - pow(nu, 2.0));       // composite plane strain modulus
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalMDR::init()
{
  if (!fix_mdr_flag) {
    if (modify->get_fix_by_style("GRANULAR/MDR").size() == 0) {
      id_fix = utils::strdup("MDR");
      modify->add_fix(fmt::format("{} all GRANULAR/MDR", id_fix));
    }
    fix_mdr_flag = 1;
  }

  // initialize particle history variables
  int tmp1, tmp2;
  index_Ro = atom->find_custom("Ro", tmp1, tmp2);                       // initial radius
  index_Vcaps = atom->find_custom("Vcaps", tmp1, tmp2);                 // spherical cap volume from intersection of apparent radius particle and contact planes
  index_Vgeo = atom->find_custom("Vgeo", tmp1, tmp2);                   // geometric particle volume of apparent particle after removing spherical cap volume
  index_Velas = atom->find_custom("Velas", tmp1, tmp2);                 // particle volume from linear elasticity
  index_eps_bar = atom->find_custom("eps_bar", tmp1, tmp2);             // volume-averaged infinitesimal strain tensor
  index_dRnumerator = atom->find_custom("dRnumerator", tmp1, tmp2);     // summation of numerator terms in calculation of dR
  index_dRdenominator = atom->find_custom("dRdenominator", tmp1, tmp2); // summation of denominator terms in calculation of dR
  index_Acon0 = atom->find_custom("Acon0", tmp1, tmp2);                 // total area involved in contacts: Acon^{n}
  index_Acon1 = atom->find_custom("Acon1", tmp1, tmp2);                 // total area involved in contacts: Acon^{n+1}
  index_Atot = atom->find_custom("Atot", tmp1, tmp2);                   // total particle area
  index_Atot_sum = atom->find_custom("Atot_sum", tmp1, tmp2);           // running sum of contact area minus cap area
  index_ddelta_bar = atom->find_custom("ddelta_bar", tmp1, tmp2);       // change in mean surface displacement
  index_psi = atom->find_custom("psi", tmp1, tmp2);                     // ratio of free surface area to total surface area
  index_sigmaxx = atom->find_custom("sigmaxx", tmp1, tmp2);             // xx-component of the stress tensor, not necessary for force calculation
  index_sigmayy = atom->find_custom("sigmayy", tmp1, tmp2);             // yy-component of the stress tensor, not necessary for force calculation
  index_sigmazz = atom->find_custom("sigmazz", tmp1, tmp2);             // zz-component of the stress tensor, not necessary for force calculation
  index_contacts = atom->find_custom("contacts", tmp1, tmp2);                     // total contacts on particle
  index_adhesive_length = atom->find_custom("adhesive_length", tmp1, tmp2);       // total contacts on particle
}

/* ---------------------------------------------------------------------- */

double GranSubModNormalMDR::calculate_forces()
{
  using namespace Granular_MDR_NS;
  // To understand the structure of the overall code it is important to consider
  // the following:
  //
  // The MDR contact model was developed by imagining individual particles being
  // squished between a number of rigid flats (references below). To allow
  // for many interacting particles, we extend the idea of isolated particles surrounded
  // by rigid flats. In particular, we imagine placing rigid flats at the overlaps
  // between particles. The force is calculated seperately on both sides
  // of the contact assuming interaction with a rigid flat. The two forces are then
  // averaged on either side of the contact to determine the final force. If the
  // contact is between a particle and wall then only one force evaluation is required.
  //
  // Zunker and Kamrin, 2024, Part I: https://doi.org/10.1016/j.jmps.2023.105492
  // Zunker and Kamrin, 2024, Part II: https://doi.org/10.1016/j.jmps.2023.105493
  // Zunker, Dunatunga, Thakur, Tang, and Kamrin, 2025:

  double *Rinitial = atom->dvector[index_Ro];
  double *Vgeo = atom->dvector[index_Vgeo];
  double *Velas = atom->dvector[index_Velas];
  double *Vcaps = atom->dvector[index_Vcaps];
  double *eps_bar = atom->dvector[index_eps_bar];
  double *dRnumerator = atom->dvector[index_dRnumerator];
  double *dRdenominator = atom->dvector[index_dRdenominator];
  double *Acon0 = atom->dvector[index_Acon0];
  double *Acon1 = atom->dvector[index_Acon1];
  double *Atot = atom->dvector[index_Atot];
  double *Atot_sum = atom->dvector[index_Atot_sum];
  double *ddelta_bar = atom->dvector[index_ddelta_bar];
  double *psi = atom->dvector[index_psi];
  double *sigmaxx = atom->dvector[index_sigmaxx];
  double *sigmayy = atom->dvector[index_sigmayy];
  double *sigmazz = atom->dvector[index_sigmazz];
  double *contacts = atom->dvector[index_contacts];
  double *adhesive_length = atom->dvector[index_adhesive_length];

  const int itag_true = atom->tag[gm->i]; // true i particle tag
  const int jtag_true = atom->tag[gm->j]; // true j particle tag
  const int i_true = gm->i;               // true i particle index
  const int j_true = gm->j;               // true j particle index
  const double radi_true = gm->radi;      // true i particle initial radius
  const double radj_true = gm->radj;      // true j particle initial radius

  F = 0.0;                                // average force
  double F0 = 0.0;                        // force on contact side 0
  double F1 = 0.0;                        // force on contact side 1
  double R0 = 0.0;
  double R1 = 0.0;
  double delta = gm->delta;               // apparent overlap

  double * history = & gm->history[history_index]; // load in all history variables

  // Rigid flat placement scheme
  double * deltamax_offset = & history[DELTA_MAX];
  double deltamax = *deltamax_offset;
  double * deltap_offset0 = & history[DELTAP_0];
  double * deltap_offset1 = & history[DELTAP_1];
  double deltap0 = *deltap_offset0;
  double deltap1 = *deltap_offset1;

  if (gm->delta >= *deltamax_offset) *deltamax_offset = gm->delta;
  deltamax = *deltamax_offset;

  double overlap_limit = 0.75;
  for (int contactSide = 0; contactSide < 2; contactSide++) {

    double *delta_offset, *deltao_offset, *delta_MDR_offset, *delta_BULK_offset;
    double *deltamax_MDR_offset, *Yflag_offset, *deltaY_offset, *cA_offset, *aAdh_offset;
    double *Ac_offset, *eps_bar_offset, *penalty_offset, *deltap_offset;

    if (contactSide == 0) {
      if (gm->contact_type == PAIR) {
        if (itag_true > jtag_true) {
          gm->i = i_true;
          gm->j = j_true;
          gm->radi = radi_true;
          gm->radj = radj_true;
        } else {
          gm->i = j_true;
          gm->j = i_true;
          gm->radi = radj_true;
          gm->radj = radi_true;
        }
        R0 = gm->radi;
        R1 = gm->radj;

        double delta_geo, delta_geo_alt;
        double delta_geoOpt1 = deltamax * (deltamax - 2.0 * R1) / (2.0 * (deltamax - R0 - R1));
        double delta_geoOpt2 = deltamax * (deltamax - 2.0 * R0) / (2.0 * (deltamax - R0 - R1));
        (R0 < gm->radj) ? delta_geo = MAX(delta_geoOpt1, delta_geoOpt2) : delta_geo = MIN(delta_geoOpt1, delta_geoOpt2);
        (R0 > gm->radj) ? delta_geo_alt = MAX(delta_geoOpt1, delta_geoOpt2) : delta_geo_alt = MIN(delta_geoOpt1, delta_geoOpt2);

        if (delta_geo / R0 > overlap_limit) {
          delta_geo = R0 * overlap_limit;
        } else if (delta_geo_alt / gm->radj > overlap_limit) {
          delta_geo = deltamax - gm->radj * overlap_limit;
        }

        double deltap = deltap0 + deltap1;
        delta = delta_geo + (deltap0 - delta_geo) / (deltap - deltamax) * (gm->delta - deltamax);
      }
    } else {
      if (gm->contact_type != PAIR) break; // contact with particle-wall requires only one evaluation
      if (itag_true < jtag_true) {
        gm->i = i_true;
        gm->j = j_true;
        gm->radi = radi_true;
        gm->radj = radj_true;
      } else {
        gm->i = j_true;
        gm->j = i_true;
        gm->radi = radj_true;
        gm->radj = radi_true;
      }

      // QUESTION: R0/R1 here are never defined so they default to zero
      //  did you mean to define:
      //    R0 = gm->radi;
      //    R1 = gm->radj;
      //  here to mirror above? After confirming, these two conditions can be easily collapsed to remove duplication and clarify differences (also replace gm->radi/radj with R0/R1)
      // ANSWER: If I am not mistaken R0 & R1 will always be set to non-zero values from the case above where contactSide = 0: 
      //            R0 = gm->radi;
      //            R1 = gm->radj;
      //         Then their values will carry over into this else statement for the case of contactSide = 1. I don't want their values to change after they are set above, 
      //         hence why I didn't redine them.
      //
      //         However, this is written in a bit of a strange way now that I am reviewing it, I will give it some more thought.
       
      double delta_geo, delta_geo_alt;
      double delta_geoOpt1 = deltamax * (deltamax - 2.0 * R1) / (2.0 * (deltamax - R0 - R1));
      double delta_geoOpt2 = deltamax * (deltamax - 2.0 * R0) / (2.0 * (deltamax - R0 - R1));
      (gm->radi < gm->radj) ? delta_geo = MAX(delta_geoOpt1,delta_geoOpt2) : delta_geo = MIN(delta_geoOpt1, delta_geoOpt2);
      (gm->radi > gm->radj) ? delta_geo_alt = MAX(delta_geoOpt1,delta_geoOpt2) : delta_geo_alt = MIN(delta_geoOpt1,delta_geoOpt2);

      if (delta_geo / gm->radi > overlap_limit) {
        delta_geo = gm->radi * overlap_limit;
      } else if (delta_geo_alt / gm->radj > overlap_limit) {
        delta_geo = deltamax - gm->radj * overlap_limit;
      }

      double deltap = deltap0 + deltap1;
      delta = delta_geo + (deltap1 - delta_geo) / (deltap - deltamax) * (gm->delta - deltamax);
    }

    delta_offset = & history[DELTA_0 + contactSide];
    deltao_offset = & history[DELTAO_0 + contactSide];
    delta_MDR_offset = & history[DELTA_MDR_0 + contactSide];
    delta_BULK_offset = & history[DELTA_BULK_0 + contactSide];
    deltamax_MDR_offset = & history[DELTAMAX_MDR_0 + contactSide];
    Yflag_offset = & history[YFLAG_0 + contactSide];
    deltaY_offset = & history[DELTAY_0 + contactSide];
    cA_offset = & history[CA_0 + contactSide];
    aAdh_offset = & history[AADH_0 + contactSide];
    Ac_offset = & history[AC_0 + contactSide];
    eps_bar_offset = & history[EPS_BAR_0 + contactSide];
    deltap_offset = & history[DELTAP_0 + contactSide];

    // temporary i and j indices
    const int i = gm->i;
    const int j = gm->j;

    // geometric property definitions
    const double Ro = Rinitial[i];              // initial radius
    const double R = gm->radi;                  // apparent radius

    // kinematics
    const double ddelta = delta - *delta_offset;
    *delta_offset = delta;

    const double deltao = delta - (R - Ro);
    const double ddeltao = deltao - *deltao_offset;
    *deltao_offset = deltao;

    double ddelta_MDR, ddelta_BULK;
    if (psi[i] < psi_b) { // if true, bulk response has triggered, split displacement increment between the MDR and BULK components
      ddelta_MDR = MIN(ddelta - ddelta_bar[i], delta - *delta_MDR_offset);
      ddelta_BULK = ddelta_bar[i];
    } else { // if false, no bulk response, full displacement increment goes to the MDR component
      ddelta_BULK = 0.0;
      ddelta_MDR = ddelta;
    }
    const double delta_MDR = *delta_MDR_offset + ddelta_MDR; // MDR displacement
    *delta_MDR_offset = delta_MDR; // Update old MDR displacement
    const double delta_BULK = MAX(0.0, *delta_BULK_offset + ddelta_BULK); // bulk displacement
    *delta_BULK_offset = delta_BULK; // update old bulk displacement

    if (delta_MDR > *deltamax_MDR_offset) *deltamax_MDR_offset = delta_MDR;
    const double deltamax_MDR = *deltamax_MDR_offset;

    const double pY = Y * (1.75 * exp(-4.4 * deltamax_MDR / R) + 1.0); // Set value of average pressure along yield surface
    if (*Yflag_offset == 0.0 && delta_MDR >= deltamax_MDR) {
    const double phertz = 4 * Eeff * sqrt(delta_MDR) / (3 * MY_PI * sqrt(R));
      if ( phertz > pY ) {
        *Yflag_offset = 1.0;
        *deltaY_offset = delta_MDR;
        *cA_offset = MY_PI * (pow(*deltaY_offset, 2) - *deltaY_offset * R);
      }
    }

    // MDR force calculation
    double F_MDR;
    double A, Ainv;               // height of elliptical indenter
    double B;                     // width of elliptical indenter
    double deltae1D;              // transformed elastic displacement
    double deltaR;                // displacement correction
    double amax, amaxsq;          // maximum experienced contact radius
    const double cA = *cA_offset; // contact area intercept

    if ( *Yflag_offset == 0.0 ) { // elastic contact
      A = 4.0 * R;
      Ainv = 1.0 / A;
      B = 2.0 * R;
      deltae1D = delta_MDR;
      amax = sqrt(deltamax_MDR * R);
    } else { // plastic contact
      amax = sqrt((2.0 * deltamax_MDR * R - pow(deltamax_MDR, 2)) + cA / MY_PI);
      amaxsq = amax * amax; // TODO later, breaks BFB
      A = 4.0 * pY / Eeff * amax;
      Ainv = 1.0 / A; // TODO later, breaks BFB
      B = 2.0 * amax;
      const double deltae1Dmax = A * 0.5; // maximum transformed elastic displacement
      const double Fmax = Eeff * (A * B * 0.25) * acos(1.0 - 2.0 * deltae1Dmax / A) - (1.0 - 2.0 * deltae1Dmax / A) * sqrt(4.0 * deltae1Dmax / A - 4.0 * pow(deltae1Dmax, 2) / pow(A, 2)); // force caused by full submersion of elliptical indenter to depth of A/2
      const double zR = R - (deltamax_MDR - deltae1Dmax); // depth of particle center
      deltaR = (Fmax * (2 * pow(amax, 2) * (-1 + nu) - (-1 + 2 * nu) * zR * (-zR + sqrt(pow(amax, 2) + pow(zR, 2)))));
      deltaR /= (MY_PI * pow(amax, 2) * 2 * G * sqrt(pow(amax, 2) + pow(zR, 2)));
      deltae1D = (delta_MDR - deltamax_MDR + deltae1Dmax + deltaR) / (1 + deltaR / deltae1Dmax);  // transformed elastic displacement

      // added for rigid flat placement
      *deltap_offset = deltamax_MDR - (deltae1Dmax + deltaR);
    }

    double a_na;
    double a_fac = 0.99;
    (deltae1D >= 0.0) ? a_na = B * sqrt(A - deltae1D) * sqrt(deltae1D) / A : a_na = 0.0;
    double aAdh = *aAdh_offset;
    if (aAdh > a_fac * amax) aAdh = a_fac * amax;

    if ( gamma > 0.0 ) { // adhesive contact
      double g_aAdh;

      if (delta_MDR == deltamax_MDR || a_na >= aAdh ) { // case 1: no tensile springs, purely compressive contact
        if (deltae1D <= 0.0) {
          F_MDR = 0.0;
        } else {
          F_MDR = Eeff * (A * B * 0.25);
          F_MDR *= acos(1.0 - 2.0 * deltae1D / A) - (1.0 - 2.0 * deltae1D / A) * sqrt(4.0 * deltae1D / A - 4.0 * pow(deltae1D, 2) / pow(A, 2));
        }
        if (std::isnan(F_MDR)) {
          error->one(FLERR, "F_MDR is NaN, case 1: no tensile springs");
        }
        *aAdh_offset = a_fac * a_na;
      } else {
        const double lmax = sqrt(2.0 * MY_PI * aAdh * gamma / Eeff);
        g_aAdh = A * 0.5 - A / B * sqrt(pow(B, 2.0) * 0.25 - pow(aAdh, 2));
        const double acrit = (-((pow(B, 2) * gamma * MY_PI) / (pow(A, 2) * Eeff)) + (pow(2, THIRD) * pow(B, 4) * pow(gamma, 2) * pow(MY_PI, FIVETHIRDS)) /
                            (pow(A, 2) * pow(Eeff, 2) * pow((27 * pow(A, 4) * pow(B, 4) * gamma) / Eeff - (2 * pow(B, 6) * pow(gamma, 3) * pow(MY_PI, 2)) / pow(Eeff, 3) + (3 * sqrt(3) * sqrt(27 * pow(A, 8) * pow(B, 8) * pow(Eeff, 2) * pow(gamma, 2) -
                              4 * pow(A, 4) * pow(B, 10) * pow(gamma, 4) * pow(MY_PI, 2))) / pow(Eeff, 2), THIRD)) + (pow(MY_PI * 0.5, THIRD) * pow((27 * pow(A, 4) * pow(B, 4) * gamma) / Eeff -
                            (2 * pow(B, 6) * pow(gamma, 3) * pow(MY_PI, 2)) / pow(Eeff, 3) + (3 * sqrt(3) * sqrt(27 * pow(A, 8) * pow(B, 8) * pow(Eeff, 2) * pow(gamma, 2) - 4 * pow(A, 4) * pow(B, 10) * pow(gamma, 4) * pow(MY_PI, 2))) /
                            pow(Eeff, 2), THIRD)) / pow(A, 2)) / 6;

        if ( (deltae1D + lmax - g_aAdh) >= 0.0) { // case 2: tensile springs, but not exceeding critical length --> deltae + lmax - g(aAdhes) >= 0
          const double deltaeAdh = g_aAdh;
          const double F_na = Eeff * (A * B * 0.25) * (acos(1.0 - 2.0 * deltaeAdh / A) - (1.0 - 2.0 * deltaeAdh / A) * sqrt(4.0 * deltaeAdh / A - 4.0 * pow(deltaeAdh, 2) / pow(A, 2)));
          const double F_Adhes = 2.0 * Eeff * (deltae1D - deltaeAdh) * aAdh;
          F_MDR = F_na + F_Adhes;
          if (std::isnan(F_MDR)) error->one(FLERR, "F_MDR is NaN, case 2: tensile springs, but not exceeding critical length");
        } else { // case 3: tensile springs exceed critical length --> deltae + lmax - g(aAdhes) = 0
          if ( aAdh < acrit ) {
            aAdh = 0.0;
            F_MDR = 0.0;
          } else {
            // newton-raphson to find aAdh
            double aAdh_tmp = aAdh;
            double fa;
            double fa2;
            double dfda;
            for (int lv1 = 0; lv1 < MDR_MAX_IT; ++lv1) {
              fa = deltae1D + sqrt(2.0 * MY_PI * aAdh_tmp * gamma / Eeff) - (A * 0.5 - A / B * sqrt(pow(B, 2) * 0.25 - pow(aAdh_tmp, 2)));
              if (abs(fa) < MDR_EPSILON1) {
                break;
              }
              dfda = -((aAdh_tmp * A) / (B * sqrt(-pow(aAdh_tmp, 2) + pow(B, 2) * 0.25))) + (gamma * sqrt(MY_PI * 0.5)) / (Eeff * sqrt((aAdh_tmp * gamma) / Eeff));
              aAdh_tmp = aAdh_tmp - fa / dfda;
              fa2 = deltae1D + sqrt(2.0 * MY_PI * aAdh_tmp * gamma / Eeff) - (A * 0.5 - A / B * sqrt(pow(B, 2) * 0.25 - pow(aAdh_tmp, 2)));
              if (abs(fa - fa2) < MDR_EPSILON2) {
                break;
              }
              if (lv1 == MDR_MAX_IT - 1){
                aAdh_tmp = 0.0;
              }
            }
            aAdh = aAdh_tmp;

            g_aAdh = A * 0.5 - A / B * sqrt(pow(B, 2) * 0.25 - pow(aAdh, 2));
            const double deltaeAdh = g_aAdh;
            const double F_na = Eeff * (A * B * 0.25) * (acos(1.0 - 2.0 * deltaeAdh / A) - (1.0 - 2.0 * deltaeAdh / A) * sqrt(4.0 * deltaeAdh / A - 4.0 * pow(deltaeAdh, 2) / pow(A, 2)));
            const double F_Adhes = 2.0 * Eeff * (deltae1D - deltaeAdh) * aAdh;
            F_MDR = F_na + F_Adhes;
            if (std::isnan(F_MDR)) error->one(FLERR, "F_MDR is NaN, case 3: tensile springs exceed critical length");
          }
          *aAdh_offset = aAdh;
        }
      }
    } else { // non-adhesive contact
      *aAdh_offset = a_na;
      (deltae1D <= 0.0) ? F_MDR = 0.0 : F_MDR = Eeff * (A * B * 0.25) * (acos(1.0 - 2.0 * deltae1D / A) - (1.0 - 2.0 * deltae1D / A) * sqrt(4.0 * deltae1D / A - 4.0 * pow(deltae1D, 2) / pow(A, 2)));
      if (std::isnan(F_MDR)) {
        error->one(FLERR, "F_MDR is NaN, non-adhesive case");
      }
    }

    contacts[i] += 1;
    adhesive_length[i] += aAdh;

    // contact penalty scheme
    penalty_offset = & history[PENALTY];
    double pij = *penalty_offset;
    const double wij = MAX(1.0 - pij, 0.0);

    // area related calculations
    double Ac;
    (*Yflag_offset == 0.0) ? Ac = MY_PI * delta * R : Ac = MY_PI * ((2.0 * delta * R - pow(delta, 2)) + cA / MY_PI);
    if (Ac < 0.0 ) Ac = 0.0;
    Atot_sum[i] += wij * (Ac - 2.0 * MY_PI * R * (deltamax_MDR + delta_BULK));
    Acon1[i] += wij * Ac;

    // bulk force calculation
    double F_BULK;
    (delta_BULK <= 0.0) ? F_BULK = 0.0 : F_BULK = (1.0 / Vgeo[i]) * Acon0[i] * delta_BULK * kappa * Ac;

    // total force calculation
    (contactSide == 0) ? F0 = F_MDR + F_BULK : F1 = F_MDR + F_BULK;

    // mean surface dipslacement calculation
    *Ac_offset = wij * Ac;

    // radius update scheme quantity calculation
    Vcaps[i] += (MY_PI * THIRD) * pow(delta, 2) * (3.0 * R - delta);

    const double Fntmp = wij * (F_MDR + F_BULK);
    const double fx = Fntmp * gm->nx[0];
    const double fy = Fntmp * gm->nx[1];
    const double fz = Fntmp * gm->nx[2];
    const double bx = -(Ro - deltao) * gm->nx[0];
    const double by = -(Ro - deltao) * gm->nx[1];
    const double bz = -(Ro - deltao) * gm->nx[2];
    const double eps_bar_contact = (1.0 / (3 * kappa * Velas[i])) * (fx * bx + fy * by + fz * bz);
    eps_bar[i] += eps_bar_contact;

    double desp_bar_contact = eps_bar_contact - *eps_bar_offset;
    if(delta_MDR == deltamax_MDR && *Yflag_offset > 0.0 && F_MDR > 0.0){
      const double Vo = (4.0 * THIRD) * MY_PI * pow(Ro, 3);
      dRnumerator[i] += -Vo * (eps_bar_contact - *eps_bar_offset) - wij * MY_PI * ddeltao * (2.0 * deltao * Ro - pow(deltao, 2) + pow(R, 2) - pow(Ro, 2));
      dRdenominator[i] += wij * 2.0 * MY_PI * R * (deltao + R - Ro);
    }
    *eps_bar_offset = eps_bar_contact;

    sigmaxx[i] += (1.0 / Velas[i]) * (fx * bx);
    sigmayy[i] += (1.0 / Velas[i]) * (fy * by);
    sigmazz[i] += (1.0 / Velas[i]) * (fz * bz);
  }

  gm->i = i_true;
  gm->j = j_true;
  gm->radi = radi_true;
  gm->radj = radj_true;

  double * penalty_offset = & history[PENALTY];
  const double pij = *penalty_offset;
  const double wij = MAX(1.0 - pij, 0.0);
  *penalty_offset = 0.0;

  // assign final force
  if (gm->contact_type != PAIR) {
    F = wij * F0;
  } else {
    F = wij * (F0 + F1) * 0.5;
  }

  // calculate damping force
  if (F > 0.0) {
    double Eeff;
    double Reff;
    if (gm->contact_type == PAIR) {
      Eeff = E / (2.0 * (1.0 - pow(nu, 2)));
      Reff = pow((1 / gm->radi + 1 / gm->radj), -1);
    } else {
      Eeff = E / (1.0 - pow(nu, 2));
      Reff = gm->radi;
    }
    const double kn = Eeff * Reff;
    const double beta = -log(CoR) / sqrt(pow(log(CoR), 2) + MY_PI * MY_PI);
    const double damp_prefactor = beta * sqrt(gm->meff * kn);
    const double F_DAMP = -damp_prefactor * gm->vnnr;

    F += wij * F_DAMP;
  }

  return F;
}

/* ---------------------------------------------------------------------- */

void GranSubModNormalMDR::set_fncrit()
{
  Fncrit = fabs(F);
}
