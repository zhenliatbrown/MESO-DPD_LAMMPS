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
   Contributing authors:
   William Zunker (MIT), Sachith Dunatunga (MIT),
   Dan Bolintineanu (SNL), Joel Clemmer (SNL)
----------------------------------------------------------------------- */

#include "fix_mdr_radius_update.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "input.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "variable.h"
#include "granular_model.h"
#include "pair_granular.h"
#include "pair.h"
#include "gran_sub_mod_normal.h"
#include "update.h"
#include "comm.h"
#include <iomanip>
#include <sstream>

using namespace LAMMPS_NS;
using namespace Granular_NS;
using namespace FixConst;
using MathConst::MY_PI;

/* ---------------------------------------------------------------------- */

FixMDRradiusUpdate::FixMDRradiusUpdate(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
 comm_forward = 20; // value needs to match number of values you communicate
}

int FixMDRradiusUpdate::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE | END_OF_STEP;
  return mask;
}

void FixMDRradiusUpdate::setup_pre_force(int /*vflag*/)
{
  int tmp1, tmp2;
  int index_Ro = atom->find_custom("Ro",tmp1,tmp2);
  int index_Vgeo = atom->find_custom("Vgeo",tmp1,tmp2);
  int index_Velas = atom->find_custom("Velas",tmp1,tmp2);
  int index_Vcaps = atom->find_custom("Vcaps",tmp1,tmp2);
  int index_eps_bar = atom->find_custom("eps_bar",tmp1,tmp2);
  int index_dRnumerator = atom->find_custom("dRnumerator",tmp1,tmp2);
  int index_dRdenominator = atom->find_custom("dRdenominator",tmp1,tmp2);
  int index_Acon0 = atom->find_custom("Acon0",tmp1,tmp2);
  int index_Acon1 = atom->find_custom("Acon1",tmp1,tmp2);
  int index_Atot = atom->find_custom("Atot",tmp1,tmp2);
  int index_Atot_sum = atom->find_custom("Atot_sum",tmp1,tmp2);
  int index_ddelta_bar = atom->find_custom("ddelta_bar",tmp1,tmp2);
  int index_psi = atom->find_custom("psi",tmp1,tmp2);
  int index_psi_b = atom->find_custom("psi_b",tmp1,tmp2);
  int index_sigmaxx = atom->find_custom("sigmaxx",tmp1,tmp2);
  int index_sigmayy = atom->find_custom("sigmayy",tmp1,tmp2);
  int index_sigmazz = atom->find_custom("sigmazz",tmp1,tmp2);
  int index_history_setup_flag = atom->find_custom("history_setup_flag",tmp1,tmp2);
  int index_contacts = atom->find_custom("contacts",tmp1,tmp2);
  int index_adhesive_length = atom->find_custom("adhesive_length",tmp1,tmp2);
  Ro = atom->dvector[index_Ro];
  Vgeo = atom->dvector[index_Vgeo];
  Velas = atom->dvector[index_Velas];
  Vcaps = atom->dvector[index_Vcaps];
  eps_bar = atom->dvector[index_eps_bar];
  dRnumerator = atom->dvector[index_dRnumerator];
  dRdenominator = atom->dvector[index_dRdenominator];
  Acon0 = atom->dvector[index_Acon0];
  Acon1 = atom->dvector[index_Acon1];
  Atot = atom->dvector[index_Atot];
  Atot_sum = atom->dvector[index_Atot_sum];
  ddelta_bar = atom->dvector[index_ddelta_bar];
  psi = atom->dvector[index_psi];
  psi_b = atom->dvector[index_psi_b];
  sigmaxx = atom->dvector[index_sigmaxx];
  sigmayy = atom->dvector[index_sigmayy];
  sigmazz = atom->dvector[index_sigmazz];
  history_setup_flag = atom->dvector[index_history_setup_flag];
  contacts = atom->dvector[index_contacts];
  adhesive_length = atom->dvector[index_adhesive_length];

  pre_force(0);
}

void FixMDRradiusUpdate::setup(int /*vflag*/)
{
  int tmp1, tmp2;
  int index_Ro = atom->find_custom("Ro",tmp1,tmp2);
  int index_Vgeo = atom->find_custom("Vgeo",tmp1,tmp2);
  int index_Velas = atom->find_custom("Velas",tmp1,tmp2);
  int index_Vcaps = atom->find_custom("Vcaps",tmp1,tmp2);
  int index_eps_bar = atom->find_custom("eps_bar",tmp1,tmp2);
  int index_dRnumerator = atom->find_custom("dRnumerator",tmp1,tmp2);
  int index_dRdenominator = atom->find_custom("dRdenominator",tmp1,tmp2);
  int index_Acon0 = atom->find_custom("Acon0",tmp1,tmp2);
  int index_Acon1 = atom->find_custom("Acon1",tmp1,tmp2);
  int index_Atot = atom->find_custom("Atot",tmp1,tmp2);
  int index_Atot_sum = atom->find_custom("Atot_sum",tmp1,tmp2);
  int index_ddelta_bar = atom->find_custom("ddelta_bar",tmp1,tmp2);
  int index_psi = atom->find_custom("psi",tmp1,tmp2);
  int index_psi_b = atom->find_custom("psi_b",tmp1,tmp2);
  int index_sigmaxx = atom->find_custom("sigmaxx",tmp1,tmp2);
  int index_sigmayy = atom->find_custom("sigmayy",tmp1,tmp2);
  int index_sigmazz = atom->find_custom("sigmazz",tmp1,tmp2);
  int index_history_setup_flag = atom->find_custom("history_setup_flag",tmp1,tmp2);
  int index_contacts = atom->find_custom("contacts",tmp1,tmp2);
  int index_adhesive_length = atom->find_custom("adhesive_length",tmp1,tmp2);
  Ro = atom->dvector[index_Ro];
  Vgeo = atom->dvector[index_Vgeo];
  Velas = atom->dvector[index_Velas];
  Vcaps = atom->dvector[index_Vcaps];
  eps_bar = atom->dvector[index_eps_bar];
  dRnumerator = atom->dvector[index_dRnumerator];
  dRdenominator = atom->dvector[index_dRdenominator];
  Acon0 = atom->dvector[index_Acon0];
  Acon1 = atom->dvector[index_Acon1];
  Atot = atom->dvector[index_Atot];
  Atot_sum = atom->dvector[index_Atot_sum];
  ddelta_bar = atom->dvector[index_ddelta_bar];
  psi = atom->dvector[index_psi];
  psi_b = atom->dvector[index_psi_b];
  sigmaxx = atom->dvector[index_sigmaxx];
  sigmayy = atom->dvector[index_sigmayy];
  sigmazz = atom->dvector[index_sigmazz];
  history_setup_flag = atom->dvector[index_history_setup_flag];
  contacts = atom->dvector[index_contacts];
  adhesive_length = atom->dvector[index_adhesive_length];

  end_of_step();
}

void FixMDRradiusUpdate::pre_force(int)
{

  PairGranular * pair = dynamic_cast<PairGranular *>(force->pair_match("granular",1));
  class GranularModel* model;
  class GranularModel** models_list = pair->models_list;
  class GranSubModNormalMDR* norm_model = nullptr;
  for (int i = 0; i < pair->nmodels; i++) {
    model = models_list[i];
    if (model->normal_model->name == "mdr") norm_model = dynamic_cast<GranSubModNormalMDR *>(model->normal_model);
  }
  if (norm_model == nullptr) error->all(FLERR, "Did not find mdr model");

  double *radius = atom->radius;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (history_setup_flag[i] < 1e-16) {
      Ro[i] = radius[i];
      Vgeo[i] = 4.0/3.0*MY_PI*pow(Ro[i],3.0);
      Velas[i] = 4.0/3.0*MY_PI*pow(Ro[i],3.0);
      Atot[i] = 4.0*MY_PI*pow(Ro[i],2.0);
      psi[i] = 1.0;
      psi_b[i] = norm_model->psi_b;
      history_setup_flag[i] = 1.0;
    }
    sigmaxx[i] = 0.0;
    sigmayy[i] = 0.0;
    sigmazz[i] = 0.0;
    contacts[i] = 0.0;
    adhesive_length[i] = 0.0;
  }
  comm->forward_comm(this);
}

int FixMDRradiusUpdate::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,int * /*pbc*/)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    buf[m++] = Ro[j];                 // 1
    buf[m++] = Vgeo[j];               // 2
    buf[m++] = Velas[j];              // 3
    buf[m++] = Vcaps[j];              // 4
    buf[m++] = eps_bar[j];            // 5
    buf[m++] = dRnumerator[j];        // 6
    buf[m++] = dRdenominator[j];      // 7
    buf[m++] = Acon0[j];              // 8
    buf[m++] = Acon1[j];              // 9
    buf[m++] = Atot[j];               // 10
    buf[m++] = Atot_sum[j];           // 11
    buf[m++] = ddelta_bar[j];         // 12
    buf[m++] = psi[j];                // 13
    buf[m++] = psi_b[j];              // 14
    buf[m++] = sigmaxx[j];            // 15
    buf[m++] = sigmayy[j];            // 16
    buf[m++] = sigmazz[j];            // 17
    buf[m++] = history_setup_flag[j]; // 18
    buf[m++] = contacts[j];           // 19
    buf[m++] = adhesive_length[j];    // 20
  }
  return m;
}

void FixMDRradiusUpdate::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {
    Ro[i] = buf[m++];                 // 1
    Vgeo[i] = buf[m++];               // 2
    Velas[i] = buf[m++];              // 3
    Vcaps[i] = buf[m++];              // 4
    eps_bar[i] = buf[m++];            // 5
    dRnumerator[i] = buf[m++];        // 6
    dRdenominator[i] = buf[m++];      // 7
    Acon0[i] = buf[m++];              // 8
    Acon1[i] = buf[m++];              // 9
    Atot[i] = buf[m++];               // 10
    Atot_sum[i] = buf[m++];           // 11
    ddelta_bar[i] = buf[m++];         // 12
    psi[i] = buf[m++];                // 13
    psi_b[i] = buf[m++];              // 14
    sigmaxx[i] = buf[m++];            // 15
    sigmayy[i] = buf[m++];            // 16
    sigmazz[i] = buf[m++];            // 17
    history_setup_flag[i] = buf[m++]; // 18
    contacts[i] = buf[m++];           // 19
    adhesive_length[i] = buf[m++];    // 20
  }
}

void FixMDRradiusUpdate::end_of_step()
{
  // update the apparent radius of every particle
  double *radius = atom->radius;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {

    const double R = radius[i];
    Atot[i] = 4.0*MY_PI*pow(R,2.0) + Atot_sum[i];

    const double Vo = 4.0/3.0*MY_PI*pow(Ro[i],3.0);
    const double Vgeoi = 4.0/3.0*MY_PI*pow(R,3.0) - Vcaps[i];
    Vgeo[i] = std::min(Vgeoi,Vo);

    const double Afree = Atot[i] - Acon1[i];
    psi[i] = Afree/Atot[i];

    const double dR = std::max(dRnumerator[i]/(dRdenominator[i] - 4.0*MY_PI*pow(R,2.0)),0.0);
    if (psi_b[i] < psi[i]) {
      if ((radius[i] + dR) < (1.5*Ro[i])) radius[i] += dR;
    }

    Velas[i] = Vo*(1.0 + eps_bar[i]);
    Vcaps[i] = 0.0;
    eps_bar[i] = 0.0;
    dRnumerator[i] = 0.0;
    dRdenominator[i] = 0.0;
    Acon0[i] = Acon1[i];
    Acon1[i] = 0.0;
    Atot_sum[i] = 0.0;
    ddelta_bar[i] = 0.0;
  }
  comm->forward_comm(this);
}
