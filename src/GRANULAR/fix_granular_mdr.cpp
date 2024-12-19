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

#include "fix_granular_mdr.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix_wall_gran_region.h"
#include "fix_neigh_history.h"
#include "force.h"
#include "granular_model.h"
#include "gran_sub_mod_normal.h"
#include "input.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "pair.h"
#include "pair_granular.h"
#include "region.h"
#include "update.h"
#include "variable.h"

#include <algorithm>
#include <iomanip>
#include <sstream>

using namespace LAMMPS_NS;
using namespace Granular_NS;
using namespace FixConst;
using MathConst::MY_PI;

enum {MEAN_SURF_DISP, RADIUS_UPDATE};

/* ---------------------------------------------------------------------- */

FixGranularMDR::FixGranularMDR(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
  comm_forward = 20; // value needs to match number of values you communicate
  id_fix = nullptr;
}

/* ---------------------------------------------------------------------- */

FixGranularMDR::~FixGranularMDR()
{
  if (id_fix) modify->delete_fix(id_fix);
}

/* ---------------------------------------------------------------------- */

int FixGranularMDR::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE | END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGranularMDR::post_constructor()
{
  int tmp1, tmp2;
  id_fix = utils::strdup("MDR_PARTICLE_HISTORY_VARIABLES");
  modify->add_fix(fmt::format("{} all property/atom d_Ro d_Vcaps d_Vgeo d_Velas d_eps_bar d_dRnumerator d_dRdenominator d_Acon0 d_Acon1 d_Atot d_Atot_sum d_ddelta_bar d_psi d_psi_b d_history_setup_flag d_sigmaxx d_sigmayy d_sigmazz d_contacts d_adhesive_length ghost yes", id_fix));

  int index_Ro = atom->find_custom("Ro",tmp1,tmp2);
  int index_Vcaps = atom->find_custom("Vcaps",tmp1,tmp2);
  int index_Vgeo = atom->find_custom("Vgeo",tmp1,tmp2);
  int index_Velas = atom->find_custom("Velas",tmp1,tmp2);
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
  int index_history_setup_flag = atom->find_custom("history_setup_flag",tmp1,tmp2);
  int index_sigmaxx = atom->find_custom("sigmaxx",tmp1,tmp2);
  int index_sigmayy = atom->find_custom("sigmayy",tmp1,tmp2);
  int index_sigmazz = atom->find_custom("sigmazz",tmp1,tmp2);
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
}

/* ---------------------------------------------------------------------- */

void FixGranularMDR::setup_pre_force(int /*vflag*/)
{
  pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixGranularMDR::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixGranularMDR::pre_force(int)
{
  radius_update();

  comm_stage = RADIUS_UPDATE;
  comm->forward_comm(this, 20);

  calculate_contact_penalty();
  mean_surf_disp();

  auto fix_list = modify->get_fix_by_style("wall/gran/region");

  for (int w = 0; w < fix_list.size(); w++) {
    update_fix_gran_wall(fix_list[w]);
  }

  comm_stage = MEAN_SURF_DISP;
  comm->forward_comm(this, 1);
}

/* ---------------------------------------------------------------------- */

int FixGranularMDR::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,int * /*pbc*/)
{
  int m = 0;
  if (comm_stage == RADIUS_UPDATE) {
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
  } else {
    for (int i = 0; i < n; i++) {
      int j = list[i];
      buf[m++] = ddelta_bar[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixGranularMDR::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;

  if (comm_stage == RADIUS_UPDATE) {
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
  } else {
    for (int i = first; i < last; i++) {
      ddelta_bar[i] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGranularMDR::end_of_step()
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
  comm_stage = RADIUS_UPDATE;
  comm->forward_comm(this, 20);
}


/* ----------------------------------------------------------------------
   calculate updated radius for atoms
------------------------------------------------------------------------- */

void FixGranularMDR::radius_update()
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
}

/* ----------------------------------------------------------------------
   ...
------------------------------------------------------------------------- */

void FixGranularMDR::calculate_contact_penalty()
{
  FixNeighHistory * fix_history = dynamic_cast<FixNeighHistory *>(modify->get_fix_by_id("NEIGH_HISTORY_GRANULAR"));
  PairGranular * pair = dynamic_cast<PairGranular *>(force->pair_match("granular",1));
  NeighList * list = pair->list;

  const int size_history = pair->get_size_history();


  int i,j,k,lv1,ii,jj,inum,jnum,itype,jtype,ktype;

  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *history_ij,*history_ik,*history_jk,*history_kj,*allhistory,*allhistory_j,*allhistory_k,**firsthistory;

  bool touchflag = false;

  double **x = atom->x;
  int *type = atom->type;
  double *radius = atom->radius;
  int nlocal = atom->nlocal;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = fix_history->firstflag;
  firsthistory = fix_history->firstvalue;

  // contact penalty calculation
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];
    allhistory = firsthistory[i];
    double radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (int jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      double radj = radius[j];
      const double delx_ij = x[j][0] - xtmp;
      const double dely_ij = x[j][1] - ytmp;
      const double delz_ij = x[j][2] - ztmp;
      const double rsq_ij = delx_ij*delx_ij + dely_ij*dely_ij + delz_ij*delz_ij;
      const double r_ij = sqrt(rsq_ij);
      const double rinv_ij = 1.0/r_ij;
      const double radsum_ij = radi + radj;
      const double deltan_ij = radsum_ij - r_ij;
      if (deltan_ij >= 0.0) {
        for (int kk = jj; kk < jnum; kk++) {
          k = jlist[kk];
          k &= NEIGHMASK;
          ktype = type[k];
          if (kk != jj) {
            const double delx_ik = x[k][0] - xtmp;
            const double dely_ik = x[k][1] - ytmp;
            const double delz_ik = x[k][2] - ztmp;
            const double rsq_ik = delx_ik*delx_ik + dely_ik*dely_ik + delz_ik*delz_ik;
            const double r_ik = sqrt(rsq_ik);
            const double rinv_ik = 1.0/r_ik;
            const double radk = radius[k];
            const double radsum_ik = radi + radk;
            const double deltan_ik = radsum_ik - r_ik;
            const double delx_jk = x[k][0] - x[j][0];
            const double dely_jk = x[k][1] - x[j][1];
            const double delz_jk = x[k][2] - x[j][2];
            const double rsq_jk = delx_jk*delx_jk + dely_jk*dely_jk + delz_jk*delz_jk;
            const double r_jk = sqrt(rsq_jk);
            const double rinv_jk = 1.0/r_jk;
            const double radsum_jk = radj + radk;
            const double deltan_jk = radsum_jk - r_jk;
            if (deltan_ik >= 0.0 && deltan_jk >= 0.0) {

              // pull ij history
              history_ij = &allhistory[size_history * jj];
              double * pij = &history_ij[22]; // penalty for contact i and j

              // pull ik history
              history_ik = &allhistory[size_history * kk];
              double * pik = &history_ik[22]; // penalty for contact i and k

              // we don't know if who owns the contact ahead of time, k might be in j's neigbor list or vice versa,
              // so we need to manually search to figure out the owner check if k is in the neighbor list of j
              double * pjk = NULL;
              int * const jklist = firstneigh[j];
              const int jknum = numneigh[j];
              for (int jk = 0; jk < jknum; jk++) {
                const int kneigh = jklist[jk] & NEIGHMASK;
                if (k == kneigh) {
                  allhistory_j = firsthistory[j];
                  history_jk = &allhistory_j[size_history * jk];
                  pjk = &history_jk[22]; // penalty for contact j and k
                  break;
                }
              }

              // check if j is in the neighbor list of k
              if (pjk == NULL) {
                int * const kjlist = firstneigh[k];
                const int kjnum = numneigh[k];
                for (int kj = 0; kj < kjnum; kj++) {
                  const int jneigh = kjlist[kj] & NEIGHMASK;
                  if (j == jneigh) {
                    allhistory_k = firsthistory[k];
                    history_kj = &allhistory_k[size_history * kj];
                    pjk = &history_kj[22]; // penalty for contact j and k
                    break;
                  }
                }
              }

              std::vector<double> distances = {r_ij,r_ik,r_jk};
              auto maxElement = std::max_element(distances.begin(), distances.end());
              double maxValue = *maxElement;
              int maxIndex = std::distance(distances.begin(), maxElement);
              if (maxIndex == 0) { // the central particle is k
                const double enx_ki = -delx_ik * rinv_ik;
                const double eny_ki = -dely_ik * rinv_ik;
                const double enz_ki = -delz_ik * rinv_ik;
                const double enx_kj = -delx_jk * rinv_jk;
                const double eny_kj = -dely_jk * rinv_jk;
                const double enz_kj = -delz_jk * rinv_jk;
                const double alpha = std::acos(enx_ki*enx_kj + eny_ki*eny_kj + enz_ki*enz_kj);
                pij[0] += 1.0/( 1.0 + std::exp(-50.0*(alpha/MY_PI - 1.0/2.0)));
              } else if (maxIndex == 1) { // the central particle is j
                const double enx_ji = -delx_ij * rinv_ij;
                const double eny_ji = -dely_ij * rinv_ij;
                const double enz_ji = -delz_ij * rinv_ij;
                const double enx_jk = delx_jk * rinv_jk;
                const double eny_jk = dely_jk * rinv_jk;
                const double enz_jk = delz_jk * rinv_jk;
                const double alpha = std::acos(enx_ji*enx_jk + eny_ji*eny_jk + enz_ji*enz_jk);
                pik[0] += 1.0/( 1.0 + std::exp(-50.0*(alpha/MY_PI - 1.0/2.0)));
              } else { // the central particle is i
                if (j < atom->nlocal || k < atom->nlocal) {
                  const double enx_ij = delx_ij * rinv_ij;
                  const double eny_ij = dely_ij * rinv_ij;
                  const double enz_ij = delz_ij * rinv_ij;
                  const double enx_ik = delx_ik * rinv_ik;
                  const double eny_ik = dely_ik * rinv_ik;
                  const double enz_ik = delz_ik * rinv_ik;
                  const double alpha = std::acos(enx_ij*enx_ik + eny_ij*eny_ik + enz_ij*enz_ik);
                  pjk[0] += 1.0/( 1.0 + std::exp(-50.0*(alpha/MY_PI - 1.0/2.0)));
                }
              }
            }
          }
        }
      }
    }
  }
}


/* ----------------------------------------------------------------------
   ...
------------------------------------------------------------------------- */

void FixGranularMDR::mean_surf_disp()
{
  FixNeighHistory * fix_history = dynamic_cast<FixNeighHistory *>(modify->get_fix_by_id("NEIGH_HISTORY_GRANULAR"));
  PairGranular * pair = dynamic_cast<PairGranular *>(force->pair_match("granular",1));
  NeighList * list = pair->list;

  const int size_history = pair->get_size_history();
  int i,j,k,ii,jj,inum,jnum,itype,jtype;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *history,*allhistory,**firsthistory;

  bool touchflag = false;
  class GranularModel* model;
  class GranularModel** models_list = pair->models_list;
  int ** types_indices = pair->types_indices;

  double **x = atom->x;
  int *type = atom->type;
  double *radius = atom->radius;
  int nlocal = atom->nlocal;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = fix_history->firstflag;
  firsthistory = fix_history->firstvalue;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    touch = firsttouch[i];
    allhistory = firsthistory[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      jtype = type[j];
      model = models_list[types_indices[itype][jtype]];

      // Reset model and copy initial geometric data
      model->xi = x[i];
      model->xj = x[j];
      model->radi = radius[i];
      model->radj = radius[j];
      model->i = i;
      model->j = j;
      model->touch = touch[jj];
      touchflag = model->check_contact();

      // is it necessary to clear the history here???
      if (!touchflag) {
        touch[jj] = 0;
        history = &allhistory[size_history * jj];
        for (k = 0; k < size_history; k++) history[k] = 0.0;
        continue;
      }

      touch[jj] = 1;

      history = &allhistory[size_history * jj];
      model->history = history;

      const double delta = model->radsum - sqrt(model->rsq);

      const int delta_offset_0 = 0;           // apparent overlap
      const int delta_offset_1 = 1;
      const int Ac_offset_0 = 18;             // contact area
      const int Ac_offset_1 = 19;
      const int deltamax_offset_ = 23;
      const int deltap_offset_0 = 24;
      const int deltap_offset_1 = 25;

      double deltamax = history[deltamax_offset_];
      double deltap0 = history[deltap_offset_0];
      double deltap1 = history[deltap_offset_1];

      if (delta > deltamax) deltamax = delta;

      double delta0old = history[delta_offset_0];
      double delta1old = history[delta_offset_1];

      int i0;
      int i1;
      if (atom->tag[i] > atom->tag[j]) {
        i0 = i;
        i1 = j;
      } else {
        i0 = j;
        i1 = i;
      }

      double R0 = radius[i0];
      double R1 = radius[i1];

      double delta_geo0;
      double delta_geo1;
      double deltaOpt1 = deltamax*(deltamax - 2.0*R1)/(2.0*(deltamax - R0 - R1));
      double deltaOpt2 = deltamax*(deltamax - 2.0*R0)/(2.0*(deltamax - R0 - R1));
      (R0 < R1) ? delta_geo0 = MAX(deltaOpt1,deltaOpt2) : delta_geo0 = MIN(deltaOpt1,deltaOpt2);
      (R0 < R1) ? delta_geo1 = MIN(deltaOpt1,deltaOpt2) : delta_geo1 = MAX(deltaOpt1,deltaOpt2);

      double overlap_limit = 0.75;

      if (delta_geo0/R0 > overlap_limit) {
        delta_geo0 = R0*overlap_limit;
        delta_geo1 = deltamax - delta_geo0;
      } else if (delta_geo1/R1 > overlap_limit) {
        delta_geo1 = R1*overlap_limit;
        delta_geo0 = deltamax - delta_geo1;
      }

      double deltap = deltap0 + deltap1;

      double delta0 = delta_geo0 + (deltap0 - delta_geo0)/(deltap - deltamax)*(delta-deltamax);
      double delta1 = delta_geo1 + (deltap1 - delta_geo1)/(deltap - deltamax)*(delta-deltamax);

      double ddel0 = delta0 - delta0old;
      double ddel1 = delta1 - delta1old;

      if (Acon0[i0] != 0.0) {
        const double Ac_offset0 = history[Ac_offset_0];
        ddelta_bar[i0] += Ac_offset0/Acon0[i0]*ddel0;
      }

      if (Acon0[i1] != 0.0) {
        const double Ac_offset1 = history[Ac_offset_1];
        ddelta_bar[i1] += Ac_offset1/Acon0[i1]*ddel1;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   Update instance of fix gran/wall
------------------------------------------------------------------------- */

void FixGranularMDR::update_fix_gran_wall(Fix* fix_in)
{
  FixWallGranRegion* fix = dynamic_cast<FixWallGranRegion*>(fix_in);
  GranularModel * model = fix->model;
  Region * region = fix->region;

  const int size_history = model->size_history;
  int i, m, nc, iwall;
  double vwall[3];
  bool touchflag = false;
  int history_update = 1;
  model->history_update = history_update;

  int regiondynamic = region->dynamic_check();
  if (!regiondynamic) vwall[0] = vwall[1] = vwall[2] = 00;

  double **x = atom->x;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (regiondynamic) {
    region->prematch();
    region->set_velocity();
  }

  if (fix->peratom_flag) fix->clear_stored_contacts();

  model->radj = 0.0;

  for (i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if (! region->match(x[i][0], x[i][1], x[i][2])) continue;

    nc = region->surface(x[i][0], x[i][1], x[i][2], radius[i] + model->pulloff_distance(radius[i], 0.0));

    if (nc == 0) {
      fix->ncontact[i] = 0;
      continue;
    }
    if (nc == 1) {
      fix->c2r[0] = 0;
      iwall = region->contact[0].iwall;
      if (fix->ncontact[i] == 0) {
        fix->ncontact[i] = 1;
        fix->walls[i][0] = iwall;
        for (m = 0; m < size_history; m++) fix->history_many[i][0][m] = 0.0;
      } else if (fix->ncontact[i] > 1 || iwall != fix->walls[i][0])
        fix->update_contacts(i, nc);
    } else
      fix->update_contacts(i, nc);


    // process current contacts
    for (int ic = 0; ic < nc; ic++) {

      // Reset model and copy initial geometric data
      model->dx[0] = region->contact[ic].delx;
      model->dx[1] = region->contact[ic].dely;
      model->dx[2] = region->contact[ic].delz;
      model->radi = radius[i];
      model->radj = region->contact[ic].radius;
      model->r = region->contact[ic].r;

      if (model->beyond_contact) model->touch = fix->history_many[i][fix->c2r[ic]][0];

      touchflag = model->check_contact();

      const double wij = 1.0;

      if (Acon0[i] != 0.0) {
        const double delta = model->radsum - model->r;
        const double delta_offset0 = fix->history_many[i][fix->c2r[ic]][0];
        const double ddelta = delta - delta_offset0;
        const double Ac_offset0 = fix->history_many[i][fix->c2r[ic]][18];
        ddelta_bar[i] += wij*Ac_offset0/Acon0[i]*ddelta;
      }
    }
  }
}
