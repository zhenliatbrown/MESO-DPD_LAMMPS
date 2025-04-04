// clang-format off
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
   Contributing authors: Mario Orsi (QMUL), m.orsi@qmul.ac.uk
                         Samuel Genheden (University of Southampton)
------------------------------------------------------------------------- */

#include "pair_lj_sf_dipole_sf.h"

#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

static int warn_single = 0;

/* ---------------------------------------------------------------------- */

PairLJSFDipoleSF::~PairLJSFDipoleSF()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(cut_coul);
    memory->destroy(cut_coulsq);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(scale);
  }
}

/* ---------------------------------------------------------------------- */

void PairLJSFDipoleSF::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fx,fy,fz;
  double rsq,rinv,r2inv,r6inv,r3inv,r5inv;
  double forcecoulx,forcecouly,forcecoulz,crossx,crossy,crossz;
  double tixcoul,tiycoul,tizcoul,tjxcoul,tjycoul,tjzcoul;
  double fq,pdotp,pidotr,pjdotr,pre1,pre2,pre3,pre4;
  double forcelj,factor_coul,factor_lj;
  double presf,afac,bfac,pqfac,qpfac,forceljcut,forceljsf;
  double aforcecoulx,aforcecouly,aforcecoulz;
  double bforcecoulx,bforcecouly,bforcecoulz;
  double rcutlj2inv, rcutcoul2inv,rcutlj6inv;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  double **mu = atom->mu;
  double **torque = atom->torque;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;
        rinv = sqrt(r2inv);

        // atom can have both a charge and dipole
        // i,j = charge-charge, dipole-dipole, dipole-charge, or charge-dipole

        forcecoulx = forcecouly = forcecoulz = 0.0;
        tixcoul = tiycoul = tizcoul = 0.0;
        tjxcoul = tjycoul = tjzcoul = 0.0;

        if (rsq < cut_coulsq[itype][jtype]) {

          rcutcoul2inv=1.0/cut_coulsq[itype][jtype];
          if (qtmp != 0.0 && q[j] != 0.0) {
            pre1 = qtmp*q[j]*rinv*(r2inv-rcutcoul2inv);

            forcecoulx += pre1*delx;
            forcecouly += pre1*dely;
            forcecoulz += pre1*delz;
          }

          if (mu[i][3] > 0.0 && mu[j][3] > 0.0) {
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;

            pdotp = mu[i][0]*mu[j][0] + mu[i][1]*mu[j][1] + mu[i][2]*mu[j][2];
            pidotr = mu[i][0]*delx + mu[i][1]*dely + mu[i][2]*delz;
            pjdotr = mu[j][0]*delx + mu[j][1]*dely + mu[j][2]*delz;

            afac = 1.0 - rsq*rsq * rcutcoul2inv*rcutcoul2inv;
            pre1 = afac * ( pdotp - 3.0 * r2inv * pidotr * pjdotr );
            aforcecoulx = pre1*delx;
            aforcecouly = pre1*dely;
            aforcecoulz = pre1*delz;

            bfac = 1.0 - 4.0*rsq*sqrt(rsq*rcutcoul2inv)*rcutcoul2inv +
              3.0*rsq*rsq*rcutcoul2inv*rcutcoul2inv;
            presf = 2.0 * r2inv * pidotr * pjdotr;
            bforcecoulx = bfac * (pjdotr*mu[i][0]+pidotr*mu[j][0]-presf*delx);
            bforcecouly = bfac * (pjdotr*mu[i][1]+pidotr*mu[j][1]-presf*dely);
            bforcecoulz = bfac * (pjdotr*mu[i][2]+pidotr*mu[j][2]-presf*delz);

            forcecoulx += 3.0 * r5inv * ( aforcecoulx + bforcecoulx );
            forcecouly += 3.0 * r5inv * ( aforcecouly + bforcecouly );
            forcecoulz += 3.0 * r5inv * ( aforcecoulz + bforcecoulz );

            pre2 = 3.0 * bfac * r5inv * pjdotr;
            pre3 = 3.0 * bfac * r5inv * pidotr;
            pre4 = -bfac * r3inv;

            crossx = pre4 * (mu[i][1]*mu[j][2] - mu[i][2]*mu[j][1]);
            crossy = pre4 * (mu[i][2]*mu[j][0] - mu[i][0]*mu[j][2]);
            crossz = pre4 * (mu[i][0]*mu[j][1] - mu[i][1]*mu[j][0]);

            tixcoul += crossx + pre2 * (mu[i][1]*delz - mu[i][2]*dely);
            tiycoul += crossy + pre2 * (mu[i][2]*delx - mu[i][0]*delz);
            tizcoul += crossz + pre2 * (mu[i][0]*dely - mu[i][1]*delx);
            tjxcoul += -crossx + pre3 * (mu[j][1]*delz - mu[j][2]*dely);
            tjycoul += -crossy + pre3 * (mu[j][2]*delx - mu[j][0]*delz);
            tjzcoul += -crossz + pre3 * (mu[j][0]*dely - mu[j][1]*delx);
          }

          if (mu[i][3] > 0.0 && q[j] != 0.0) {
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;
            pidotr = mu[i][0]*delx + mu[i][1]*dely + mu[i][2]*delz;
            pre1 = 3.0 * q[j] * r5inv * pidotr * (1-rsq*rcutcoul2inv);
            pqfac = 1.0 - 3.0*rsq*rcutcoul2inv +
              2.0*rsq*sqrt(rsq*rcutcoul2inv)*rcutcoul2inv;
            pre2 = q[j] * r3inv * pqfac;

            forcecoulx += pre2*mu[i][0] - pre1*delx;
            forcecouly += pre2*mu[i][1] - pre1*dely;
            forcecoulz += pre2*mu[i][2] - pre1*delz;
            tixcoul += pre2 * (mu[i][1]*delz - mu[i][2]*dely);
            tiycoul += pre2 * (mu[i][2]*delx - mu[i][0]*delz);
            tizcoul += pre2 * (mu[i][0]*dely - mu[i][1]*delx);
          }

          if (mu[j][3] > 0.0 && qtmp != 0.0) {
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;
            pjdotr = mu[j][0]*delx + mu[j][1]*dely + mu[j][2]*delz;
            pre1 = 3.0 * qtmp * r5inv * pjdotr * (1-rsq*rcutcoul2inv);
            qpfac = 1.0 - 3.0*rsq*rcutcoul2inv +
              2.0*rsq*sqrt(rsq*rcutcoul2inv)*rcutcoul2inv;
            pre2 = qtmp * r3inv * qpfac;

            forcecoulx += pre1*delx - pre2*mu[j][0];
            forcecouly += pre1*dely - pre2*mu[j][1];
            forcecoulz += pre1*delz - pre2*mu[j][2];
            tjxcoul += -pre2 * (mu[j][1]*delz - mu[j][2]*dely);
            tjycoul += -pre2 * (mu[j][2]*delx - mu[j][0]*delz);
            tjzcoul += -pre2 * (mu[j][0]*dely - mu[j][1]*delx);
          }
        }

        // LJ interaction

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv*r2inv*r2inv;
          forceljcut = r6inv*(lj1[itype][jtype]*r6inv-lj2[itype][jtype])*r2inv;

          rcutlj2inv = 1.0 / cut_ljsq[itype][jtype];
          rcutlj6inv = rcutlj2inv * rcutlj2inv * rcutlj2inv;
          forceljsf = (lj1[itype][jtype]*rcutlj6inv - lj2[itype][jtype]) *
            rcutlj6inv * rcutlj2inv;

          forcelj = factor_lj * (forceljcut - forceljsf);
        } else forcelj = 0.0;

        // total force

        fq = factor_coul*qqrd2e*scale[itype][jtype];
        fx = fq*forcecoulx + delx*forcelj;
        fy = fq*forcecouly + dely*forcelj;
        fz = fq*forcecoulz + delz*forcelj;

        // force & torque accumulation

        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        torque[i][0] += fq*tixcoul;
        torque[i][1] += fq*tiycoul;
        torque[i][2] += fq*tizcoul;

        if (newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
          torque[j][0] += fq*tjxcoul;
          torque[j][1] += fq*tjycoul;
          torque[j][2] += fq*tjzcoul;
        }

        if (eflag) {
          if (rsq < cut_coulsq[itype][jtype]) {
            ecoul = (1.0-sqrt(rsq/cut_coulsq[itype][jtype]));
            ecoul *= ecoul;
            ecoul *= qtmp * q[j] * rinv;
            if (mu[i][3] > 0.0 && mu[j][3] > 0.0)
              ecoul += bfac * (r3inv*pdotp - 3.0*r5inv*pidotr*pjdotr);
            if (mu[i][3] > 0.0 && q[j] != 0.0)
              ecoul += -q[j] * r3inv * pqfac * pidotr;
            if (mu[j][3] > 0.0 && qtmp != 0.0)
              ecoul += qtmp * r3inv * qpfac * pjdotr;
            ecoul *= factor_coul*qqrd2e*scale[itype][jtype];
          } else ecoul = 0.0;

          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype])+
              rcutlj6inv*(6*lj3[itype][jtype]*rcutlj6inv-3*lj4[itype][jtype])*
              rsq*rcutlj2inv+
              rcutlj6inv*(-7*lj3[itype][jtype]*rcutlj6inv+4*lj4[itype][jtype]);
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                                 evdwl,ecoul,fx,fy,fz,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJSFDipoleSF::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(scale,n+1,n+1,"pair:scale");
  for (int i = 1; i <= n; i++) {
    for (int j = i; j <= n; j++) {
      setflag[i][j] = 0;
      scale[i][j] = 1.0;
    }
  }

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut_lj,n+1,n+1,"pair:cut_lj");
  memory->create(cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  memory->create(cut_coul,n+1,n+1,"pair:cut_coul");
  memory->create(cut_coulsq,n+1,n+1,"pair:cut_coulsq");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJSFDipoleSF::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2)
    error->all(FLERR,"Incorrect args in pair_style command");

  if (strcmp(update->unit_style,"electron") == 0)
    error->all(FLERR,"Cannot (yet) use 'electron' units with dipoles");

  cut_lj_global = utils::numeric(FLERR,arg[0],false,lmp);
  if (narg == 1) cut_coul_global = cut_lj_global;
  else cut_coul_global = utils::numeric(FLERR,arg[1],false,lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_lj[i][j] = cut_lj_global;
          cut_coul[i][j] = cut_coul_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJSFDipoleSF::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 8)
    error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double epsilon_one = utils::numeric(FLERR,arg[2],false,lmp);
  double sigma_one = utils::numeric(FLERR,arg[3],false,lmp);

  double cut_lj_one = cut_lj_global;
  double cut_coul_one = cut_coul_global;
  double scale_one = 1.0;
  int iarg = 4;

  if ((narg > iarg) && (strcmp(arg[iarg],"scale") != 0)) {
    cut_coul_one = cut_lj_one = utils::numeric(FLERR,arg[iarg],false,lmp);
    ++iarg;
  }
  if ((narg > iarg) && (strcmp(arg[iarg],"scale") != 0)) {
    cut_coul_one = utils::numeric(FLERR,arg[iarg],false,lmp);
    ++iarg;
  }
  if (narg > iarg) {
    if (strcmp(arg[iarg],"scale") == 0) {
      scale_one = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));
  }
  if (iarg != narg)
    error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut_lj[i][j] = cut_lj_one;
      cut_coul[i][j] = cut_coul_one;
      setflag[i][j] = 1;
      scale[i][j] = scale_one;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJSFDipoleSF::init_style()
{
  if (!atom->q_flag || !atom->mu_flag || !atom->torque_flag)
    error->all(FLERR,"Pair dipole/sf requires atom attributes q, mu, torque");

  neighbor->add_request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJSFDipoleSF::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut_lj[i][j] = mix_distance(cut_lj[i][i],cut_lj[j][j]);
    cut_coul[i][j] = mix_distance(cut_coul[i][i],cut_coul[j][j]);
  }

  double cut = MAX(cut_lj[i][j],cut_coul[i][j]);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];
  cut_coulsq[i][j] = cut_coul[i][j] * cut_coul[i][j];

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);

  cut_ljsq[j][i] = cut_ljsq[i][j];
  cut_coulsq[j][i] = cut_coulsq[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  scale[j][i] = scale[i][j];

  return cut;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJSFDipoleSF::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&cut_lj[i][j],sizeof(double),1,fp);
        fwrite(&cut_coul[i][j],sizeof(double),1,fp);
        fwrite(&scale[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJSFDipoleSF::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&epsilon[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&sigma[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_lj[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&cut_coul[i][j],sizeof(double),1,fp,nullptr,error);
          utils::sfread(FLERR,&scale[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_coul[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&scale[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJSFDipoleSF::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul_global,sizeof(double),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJSFDipoleSF::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_lj_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&cut_coul_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

// PairLJSFDipoleSF: calculation of force is missing (to be implemented)
double PairLJSFDipoleSF::single(int i, int j, int itype, int jtype, double rsq,
                                double factor_coul, double factor_lj,
                                double &fforce)
{
  double r2inv,r6inv;
  double pdotp,pidotr,pjdotr,delx,dely,delz;
  double rinv, r3inv,r5inv, rcutlj2inv, rcutcoul2inv,rcutlj6inv;
  double qtmp,xtmp,ytmp,ztmp,bfac,pqfac,qpfac, ecoul, evdwl;

  double **x = atom->x;
  double *q = atom->q;
  double **mu = atom->mu;

  if (!warn_single) {
    warn_single = 1;
    if (comm->me == 0) {
      error->warning(FLERR,"Single method for lj/sf/dipole/sf does not compute forces");
    }
  }
  qtmp = q[i];
  xtmp = x[i][0];
  ytmp = x[i][1];
  ztmp = x[i][2];

  r2inv = 1.0/rsq;
  rinv = sqrt(r2inv);
  fforce = 0.0;

  if (rsq < cut_coulsq[itype][jtype]) {
    delx = xtmp - x[j][0];
    dely = ytmp - x[j][1];
    delz = ztmp - x[j][2];
    // if (qtmp != 0.0 && q[j] != 0.0) {
    //   pre1 = qtmp*q[j]*rinv*(r2inv-1.0/cut_coulsq[itype][jtype]);
    //   forcecoulx += pre1*delx;
    //   forcecouly += pre1*dely;
    //   forcecoulz += pre1*delz;
    // }
    if (mu[i][3] > 0.0 && mu[j][3] > 0.0) {
      r3inv = r2inv*rinv;
      r5inv = r3inv*r2inv;
      rcutcoul2inv=1.0/cut_coulsq[itype][jtype];
      pdotp = mu[i][0]*mu[j][0] + mu[i][1]*mu[j][1] + mu[i][2]*mu[j][2];
      pidotr = mu[i][0]*delx + mu[i][1]*dely + mu[i][2]*delz;
      pjdotr = mu[j][0]*delx + mu[j][1]*dely + mu[j][2]*delz;
      bfac = 1.0 - 4.0*rsq*sqrt(rsq)*rcutcoul2inv*sqrt(rcutcoul2inv) +
        3.0*rsq*rsq*rcutcoul2inv*rcutcoul2inv;
    }
    if (mu[i][3] > 0.0 && q[j] != 0.0) {
      r3inv = r2inv*rinv;
      r5inv = r3inv*r2inv;
      pidotr = mu[i][0]*delx + mu[i][1]*dely + mu[i][2]*delz;
      rcutcoul2inv=1.0/cut_coulsq[itype][jtype];
      pqfac = 1.0 - 3.0*rsq*rcutcoul2inv +
        2.0*rsq*sqrt(rsq)*rcutcoul2inv*sqrt(rcutcoul2inv);
    }
    if (mu[j][3] > 0.0 && qtmp != 0.0) {
      r3inv = r2inv*rinv;
      r5inv = r3inv*r2inv;
      pjdotr = mu[j][0]*delx + mu[j][1]*dely + mu[j][2]*delz;
      rcutcoul2inv=1.0/cut_coulsq[itype][jtype];
      qpfac = 1.0 - 3.0*rsq*rcutcoul2inv +
        2.0*rsq*sqrt(rsq)*rcutcoul2inv*sqrt(rcutcoul2inv);
    }
  }
  if (rsq < cut_ljsq[itype][jtype]) {
    r6inv = r2inv*r2inv*r2inv;
    rcutlj2inv = 1.0 / cut_ljsq[itype][jtype];
    rcutlj6inv = rcutlj2inv * rcutlj2inv * rcutlj2inv;
  }

  double eng = 0.0;
  if (rsq < cut_coulsq[itype][jtype]) {
    ecoul = (1.0-sqrt(rsq)/sqrt(cut_coulsq[itype][jtype]));
    ecoul *= ecoul;
    ecoul *= qtmp * q[j] * rinv;
    if (mu[i][3] > 0.0 && mu[j][3] > 0.0)
      ecoul += bfac * (r3inv*pdotp - 3.0*r5inv*pidotr*pjdotr);
    if (mu[i][3] > 0.0 && q[j] != 0.0)
      ecoul += -q[j] * r3inv * pqfac * pidotr;
    if (mu[j][3] > 0.0 && qtmp != 0.0)
      ecoul += qtmp * r3inv * qpfac * pjdotr;
    ecoul *= factor_coul*force->qqrd2e*scale[itype][jtype];
    eng += ecoul;
  }
  if (rsq < cut_ljsq[itype][jtype]) {
    evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype])+
      rcutlj6inv*(6*lj3[itype][jtype]*rcutlj6inv-3*lj4[itype][jtype])*
      rsq*rcutlj2inv+
      rcutlj6inv*(-7*lj3[itype][jtype]*rcutlj6inv+4*lj4[itype][jtype]);
    eng += evdwl*factor_lj;
  }

  return eng;
}

/* ---------------------------------------------------------------------- */

void *PairLJSFDipoleSF::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  if (strcmp(str,"scale") == 0) return (void *) scale;
  if (strcmp(str,"cut_coul") == 0) return (void *) cut_coul;
  return nullptr;
}
