// clang-format off
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

/* ----------------------------------------------------------------------
   Contributing authors: Christian Trott (SNL), Stan Moore (SNL)
------------------------------------------------------------------------- */

#include "sna_kokkos.h"
#include "memory_kokkos.h"
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <type_traits>

namespace LAMMPS_NS {

static const double MY_PI  = 3.14159265358979323846; // pi
static const double MY_PI2  = 1.57079632679489661923; // pi/2

template<class DeviceType, typename real_type, int vector_length>
template<class CopyClass>
inline
SNAKokkos<DeviceType, real_type, vector_length>::SNAKokkos(const CopyClass& copy)
  : twojmax(copy.twojmax), d_coeffelem(copy.d_coeffelem), rmin0(copy.rmin0),
    rfac0(copy.rfac0), switch_flag(copy.switchflag), switch_inner_flag(copy.switchinnerflag),
    chem_flag(copy.chemflag), bnorm_flag(copy.bnormflag), wselfall_flag(copy.wselfallflag),
    quadratic_flag(copy.quadraticflag), bzero_flag(copy.bzeroflag)
{
  wself = static_cast<real_type>(1.0);

  if (chem_flag)
    nelements = copy.nelements;
  else
    nelements = 1;

  ncoeff = compute_ncoeff();

  nmax = 0;
  natom = 0;
  natom_pad = 0;

  build_indexlist();

  int jdimpq = twojmax + 2;
  MemKK::realloc_kokkos(rootpqarray,"SNAKokkos::rootpqarray",jdimpq,jdimpq);

  MemKK::realloc_kokkos(cglist,"SNAKokkos::cglist",idxcg_max);

  if (bzero_flag) {
    MemKK::realloc_kokkos(bzero,"sna:bzero",twojmax+1);
    auto h_bzero = Kokkos::create_mirror_view(bzero);

    double www = wself*wself*wself;
    for (int j = 0; j <= twojmax; j++)
      if (bnorm_flag)
        h_bzero[j] = www;
      else
        h_bzero[j] = www*(j+1);
    Kokkos::deep_copy(bzero,h_bzero);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
SNAKokkos<DeviceType, real_type, vector_length>::~SNAKokkos()
{
}

template<class DeviceType, typename real_type, int vector_length>
inline
void SNAKokkos<DeviceType, real_type, vector_length>::build_indexlist()
{
  // index list for cglist

  int jdim = twojmax + 1;
  MemKK::realloc_kokkos(idxcg_block,"SNAKokkos::idxcg_block",jdim,jdim,jdim);
  auto h_idxcg_block = Kokkos::create_mirror_view(idxcg_block);

  int idxcg_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        h_idxcg_block(j1,j2,j) = idxcg_count;
        for (int m1 = 0; m1 <= j1; m1++)
          for (int m2 = 0; m2 <= j2; m2++)
            idxcg_count++;
      }
  idxcg_max = idxcg_count;
  Kokkos::deep_copy(idxcg_block,h_idxcg_block);

  // index list for uarray
  // need to include both halves

  MemKK::realloc_kokkos(idxu_block,"SNAKokkos::idxu_block",jdim);
  auto h_idxu_block = Kokkos::create_mirror_view(idxu_block);

  int idxu_count = 0;

  for (int j = 0; j <= twojmax; j++) {
    h_idxu_block[j] = idxu_count;
    for (int mb = 0; mb <= j; mb++)
      for (int ma = 0; ma <= j; ma++)
        idxu_count++;
  }
  idxu_max = idxu_count;
  Kokkos::deep_copy(idxu_block,h_idxu_block);

  // index list for half uarray
  MemKK::realloc_kokkos(idxu_half_block,"SNAKokkos::idxu_half_block",jdim);
  auto h_idxu_half_block = Kokkos::create_mirror_view(idxu_half_block);

  int idxu_half_count = 0;
  for (int j = 0; j <= twojmax; j++) {
    h_idxu_half_block[j] = idxu_half_count;
    for (int mb = 0; 2*mb <= j; mb++)
      for (int ma = 0; ma <= j; ma++)
        idxu_half_count++;
  }
  idxu_half_max = idxu_half_count;
  Kokkos::deep_copy(idxu_half_block, h_idxu_half_block);

  // mapping between full and half indexing, encoding flipping
  MemKK::realloc_kokkos(idxu_full_half,"SNAKokkos::idxu_full_half",idxu_max);
  auto h_idxu_full_half = Kokkos::create_mirror_view(idxu_full_half);

  idxu_count = 0;
  for (int j = 0; j <= twojmax; j++) {
    int jju_half = h_idxu_half_block[j];
    for (int mb = 0; mb <= j; mb++) {
      for (int ma = 0; ma <= j; ma++) {
        FullHalfMapper mapper;
        if (2*mb <= j) {
          mapper.idxu_half = jju_half + mb * (j + 1) + ma;
          mapper.flip_sign = 0;
        } else {
          mapper.idxu_half = jju_half + (j + 1 - mb) * (j + 1) - (ma + 1);
          mapper.flip_sign = (((ma+mb)%2==0)?1:-1);
        }
        h_idxu_full_half[idxu_count] = mapper;
        idxu_count++;
      }
    }
  }

  Kokkos::deep_copy(idxu_full_half, h_idxu_full_half);

  // index list for "cache" uarray
  // this is the GPU scratch memory requirements
  // applied the CPU structures
  MemKK::realloc_kokkos(idxu_cache_block,"SNAKokkos::idxu_cache_block",jdim);
  auto h_idxu_cache_block = Kokkos::create_mirror_view(idxu_cache_block);

  int idxu_cache_count = 0;
  for (int j = 0; j <= twojmax; j++) {
    h_idxu_cache_block[j] = idxu_cache_count;
    for (int mb = 0; mb < ((j+3)/2); mb++)
      for (int ma = 0; ma <= j; ma++)
        idxu_cache_count++;
  }
  idxu_cache_max = idxu_cache_count;
  Kokkos::deep_copy(idxu_cache_block, h_idxu_cache_block);

  // index list for beta and B

  int idxb_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) idxb_count++;

  idxb_max = idxb_count;
  MemKK::realloc_kokkos(idxb,"SNAKokkos::idxb",idxb_max);
  auto h_idxb = Kokkos::create_mirror_view(idxb);

  idxb_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) {
          h_idxb(idxb_count,0) = j1;
          h_idxb(idxb_count,1) = j2;
          h_idxb(idxb_count,2) = j;
          idxb_count++;
        }
  Kokkos::deep_copy(idxb,h_idxb);

  // reverse index list for beta and b

  MemKK::realloc_kokkos(idxb_block,"SNAKokkos::idxb_block",jdim,jdim,jdim);
  auto h_idxb_block = Kokkos::create_mirror_view(idxb_block);

  idxb_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        if (j >= j1) {
          h_idxb_block(j1,j2,j) = idxb_count;
          idxb_count++;
        }
      }
  Kokkos::deep_copy(idxb_block,h_idxb_block);

  // index list for zlist

  int idxz_count = 0;

  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++)
            idxz_count++;

  idxz_max = idxz_count;
  MemKK::realloc_kokkos(idxz,"SNAKokkos::idxz",idxz_max);
  auto h_idxz = Kokkos::create_mirror_view(idxz);

  MemKK::realloc_kokkos(idxz_block,"SNAKokkos::idxz_block", jdim,jdim,jdim);
  auto h_idxz_block = Kokkos::create_mirror_view(idxz_block);

  // fused table for idxz and idxcg_block
  idxz_count = 0;
  idxcg_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        h_idxz_block(j1,j2,j) = idxz_count;

        // find right beta(ii,jjb) entry
        // multiply and divide by j+1 factors
        // account for multiplicity of 1, 2, or 3

        for (int mb = 0; 2*mb <= j; mb++) {
          for (int ma = 0; ma <= j; ma++) {
            int ma1min = MAX(0, (2 * ma - j - j2 + j1) / 2);
            int ma2max = (2 * ma - j - (2 * ma1min - j1) + j2) / 2;
            int mb1min = MAX(0, (2 * mb - j - j2 + j1) / 2);
            int mb2max = (2 * mb - j - (2 * mb1min - j1) + j2) / 2;
            int na = MIN(j1, (2 * ma - j + j2 + j1) / 2) - ma1min + 1;
            int nb = MIN(j1, (2 * mb - j + j2 + j1) / 2) - mb1min + 1;

            // apply to z(j1,j2,j,ma,mb) to unique element of y(j)
            // ylist is "compressed" via symmetry in its
            // contraction with dulist
            const int jju_half = h_idxu_half_block[j] + (j+1)*mb + ma;

            // idxz_struct's constructor handles all of the data packing
            h_idxz(idxz_count) = idxz_struct(j1, j2, j, ma1min, ma2max, mb1min, mb2max, na, nb, jju_half, idxcg_count);

            idxz_count++;
          }
        }

        // there are different loop bounds for idxcg_count
        for (int m1 = 0; m1 <= j1; m1++)
          for (int m2 = 0; m2 <= j2; m2++)
            idxcg_count++;
      }

  Kokkos::deep_copy(idxz,h_idxz);
  Kokkos::deep_copy(idxz_block,h_idxz_block);

}

/* ---------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
inline
void SNAKokkos<DeviceType, real_type, vector_length>::init()
{
  init_clebsch_gordan();
  init_rootpqarray();
}

template<class DeviceType, typename real_type, int vector_length>
inline
void SNAKokkos<DeviceType, real_type, vector_length>::grow_rij(int newnatom, int newnmax)
{
  if (newnatom <= natom && newnmax <= nmax) return;
  natom = newnatom;
  // Create padded structures
  const int natom_div = (natom + vector_length - 1) / vector_length;
  natom_pad = natom_div * vector_length;
  nmax = newnmax;

  MemKK::realloc_kokkos(rij,"sna:rij",natom_pad,nmax,3);
  MemKK::realloc_kokkos(wj,"sna:wj",natom_pad,nmax);
  MemKK::realloc_kokkos(rcutij,"sna:rcutij",natom_pad,nmax);
  MemKK::realloc_kokkos(sinnerij,"sna:sinnerij",natom_pad,nmax);
  MemKK::realloc_kokkos(dinnerij,"sna:dinnerij",natom_pad,nmax);
  MemKK::realloc_kokkos(inside,"sna:inside",natom_pad,nmax);
  MemKK::realloc_kokkos(element,"sna:element",natom_pad,nmax);
  MemKK::realloc_kokkos(dedr,"sna:dedr",natom_pad,nmax,3);

  MemKK::realloc_kokkos(ulisttot_re,"sna:ulisttot_re", natom_pad, nelements, idxu_half_max);
  MemKK::realloc_kokkos(ulisttot_im,"sna:ulisttot_im", natom_pad, nelements, idxu_half_max);
  MemKK::realloc_kokkos(ulisttot,"sna:ulisttot", natom_pad, nelements, idxu_max);

  MemKK::realloc_kokkos(zlist,"sna:zlist", natom_pad, ndoubles, idxz_max);
  MemKK::realloc_kokkos(blist,"sna:blist", natom_pad, ntriples, idxb_max);

  MemKK::realloc_kokkos(ylist_re,"sna:ylist_re", natom_pad, nelements, idxu_half_max);
  MemKK::realloc_kokkos(ylist_im,"sna:ylist_im", natom_pad, nelements, idxu_half_max);

  if constexpr (!host_flag) {
    MemKK::realloc_kokkos(a_gpu,"sna:a_gpu",natom_pad,nmax);
    MemKK::realloc_kokkos(b_gpu,"sna:b_gpu",natom_pad,nmax);
    MemKK::realloc_kokkos(da_gpu,"sna:da_gpu",natom_pad,nmax,3);
    MemKK::realloc_kokkos(db_gpu,"sna:db_gpu",natom_pad,nmax,3);
    MemKK::realloc_kokkos(sfac_gpu,"sna:sfac_gpu",natom_pad,nmax,4);
    MemKK::realloc_kokkos(ulist_cpu,"sna:ulist_cpu",1,1,1);
    MemKK::realloc_kokkos(dulist_cpu,"sna:dulist_cpu",1,1,1);
  } else {
    MemKK::realloc_kokkos(a_gpu,"sna:a_gpu",1,1);
    MemKK::realloc_kokkos(b_gpu,"sna:b_gpu",1,1);
    MemKK::realloc_kokkos(da_gpu,"sna:da_gpu",1,1,1);
    MemKK::realloc_kokkos(db_gpu,"sna:db_gpu",1,1,1);
    MemKK::realloc_kokkos(sfac_gpu,"sna:sfac_gpu",1,1,1);
    MemKK::realloc_kokkos(ulist_cpu,"sna:ulist_cpu", natom_pad, nmax, idxu_cache_max);
    MemKK::realloc_kokkos(dulist_cpu,"sna:dulist_cpu", natom_pad, nmax, idxu_cache_max);
  }
}

/* ----------------------------------------------------------------------
 * GPU routines
 * ----------------------------------------------------------------------*/


/* ----------------------------------------------------------------------
   Precompute the Cayley-Klein parameters and the derivatives thereof.
   This routine better exploits parallelism than the GPU ComputeUi and
   ComputeFusedDeidrj, which are one warp per atom-neighbor pair.
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType, real_type, vector_length>::compute_cayley_klein(const int& iatom, const int& jnbor) const
{
  const real_type x = rij(iatom,jnbor,0);
  const real_type y = rij(iatom,jnbor,1);
  const real_type z = rij(iatom,jnbor,2);
  const real_type rsq = x * x + y * y + z * z;
  const real_type r = sqrt(rsq);
  const real_type rcut = rcutij(iatom, jnbor);
  const real_type sinner = sinnerij(iatom, jnbor);
  const real_type dinner = dinnerij(iatom, jnbor);
  const real_type rscale0 = rfac0 * static_cast<real_type>(MY_PI) / (rcut - rmin0);
  const real_type theta0 = (r - rmin0) * rscale0;
  const real_type sn = sin(theta0);
  const real_type cs = cos(theta0);
  const real_type z0 = r * cs / sn;
  const real_type dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;

  const real_type wj_local = wj(iatom, jnbor);
  real_type sfac, dsfac;
  compute_s_dsfac(r, rcut, sinner, dinner, sfac, dsfac);
  sfac *= wj_local;
  dsfac *= wj_local;

  const real_type rinv = static_cast<real_type>(1.0) / r;
  const real_type ux = x * rinv;
  const real_type uy = y * rinv;
  const real_type uz = z * rinv;

  const real_type r0inv = static_cast<real_type>(1.0) / sqrt(r * r + z0 * z0);

  const complex a = { z0 * r0inv, -z * r0inv };
  const complex b = { r0inv * y, -r0inv * x };

  const real_type dr0invdr = -r0inv * r0inv * r0inv * (r + z0 * dz0dr);

  const real_type dr0invx = dr0invdr * ux;
  const real_type dr0invy = dr0invdr * uy;
  const real_type dr0invz = dr0invdr * uz;

  const real_type dz0x = dz0dr * ux;
  const real_type dz0y = dz0dr * uy;
  const real_type dz0z = dz0dr * uz;

  const complex dax = { dz0x * r0inv + z0 * dr0invx, -z * dr0invx };
  const complex day = { dz0y * r0inv + z0 * dr0invy, -z * dr0invy };
  const complex daz = { dz0z * r0inv + z0 * dr0invz, -z * dr0invz - r0inv };

  const complex dbx = { y * dr0invx, -x * dr0invx - r0inv };
  const complex dby = { y * dr0invy + r0inv, -x * dr0invy };
  const complex dbz = { y * dr0invz, -x * dr0invz };

  const real_type dsfacux = dsfac * ux;
  const real_type dsfacuy = dsfac * uy;
  const real_type dsfacuz = dsfac * uz;

  a_gpu(iatom,jnbor) = a;
  b_gpu(iatom,jnbor) = b;

  da_gpu(iatom,jnbor,0) = dax;
  db_gpu(iatom,jnbor,0) = dbx;

  da_gpu(iatom,jnbor,1) = day;
  db_gpu(iatom,jnbor,1) = dby;

  da_gpu(iatom,jnbor,2) = daz;
  db_gpu(iatom,jnbor,2) = dbz;

  sfac_gpu(iatom,jnbor,0) = sfac;
  sfac_gpu(iatom,jnbor,1) = dsfacux;
  sfac_gpu(iatom,jnbor,2) = dsfacuy;
  sfac_gpu(iatom,jnbor,3) = dsfacuz;

  // we need to explicitly zero `dedr` somewhere before hitting
  // ComputeFusedDeidrj --- this is just a convenient place to do it.
  dedr(iatom, jnbor, 0) = static_cast<real_type>(0.);
  dedr(iatom, jnbor, 1) = static_cast<real_type>(0.);
  dedr(iatom, jnbor, 2) = static_cast<real_type>(0.);

}

/* ----------------------------------------------------------------------
   Initialize ulisttot with self-energy terms.
   Ulisttot uses a "half" data layout which takes
   advantage of the symmetry of the Wigner U matrices.
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType, real_type, vector_length>::pre_ui(const int& iatom, const int& j, const int& ielem) const
{
  for (int jelem = 0; jelem < nelements; jelem++) {
    int jju_half = idxu_half_block(j);

    // Only diagonal elements get initialized
    // Top half only: gets symmetrized by TransformUi

    for (int mb = 0; 2*mb <= j; mb++) {
      for (int ma = 0; ma <= j; ma++) {

        real_type re_part = static_cast<real_type>(0.);
        if (ma == mb && (!chem_flag || ielem == jelem || wselfall_flag)) { re_part = wself; }

        ulisttot_re(iatom, jelem, jju_half) = re_part;
        ulisttot_im(iatom, jelem, jju_half) = 0;

        jju_half++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute Ui by computing Wigner U-functions for one neighbor and
   accumulating to the total. GPU only.
------------------------------------------------------------------------- */

// Version of the code that exposes additional parallelism by threading over `j_bend` values

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType, real_type, vector_length>::compute_ui_small(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, const int iatom_mod, const int j_bend, const int jnbor, const int iatom_div) const
{
  const int iatom = iatom_mod + vector_length * iatom_div;
  // get shared memory offset
  // scratch size: 32 atoms * (twojmax+1) cached values, no double buffer
  const int tile_size = vector_length * (twojmax + 1);

  const int team_rank = team.team_rank();
  const int scratch_shift = team_rank * tile_size;

  // extract and wrap
  const WignerWrapper<real_type, vector_length> ulist_wrapper((complex*)team.team_shmem().get_shmem(team.team_size() * tile_size * sizeof(complex), 0) + scratch_shift, iatom_mod);

  // load parameters
  const complex a = a_gpu(iatom, jnbor);
  const complex b = b_gpu(iatom, jnbor);
  const real_type sfac = sfac_gpu(iatom, jnbor, 0);

  const int jelem = element(iatom, jnbor);

  // we need to "choose" when to bend
  // this for loop is here for context --- we expose additional
  // parallelism over this loop instead
  //for (int j_bend = 0; j_bend <= twojmax; j_bend++) {
  evaluate_ui_jbend(ulist_wrapper, a, b, sfac, jelem, iatom, j_bend);
}

// Version of the code that loops over all `j_bend` values which reduces integer arithmetic
// and some amount of load imbalance, at the expense of reducing parallelism
template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType, real_type, vector_length>::compute_ui_large(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, const int iatom_mod, const int jnbor, const int iatom_div) const
{
  const int iatom = iatom_mod + vector_length * iatom_div;
  // get shared memory offset
  // scratch size: 32 atoms * (twojmax+1) cached values, no double buffer
  const int tile_size = vector_length * (twojmax + 1);

  const int team_rank = team.team_rank();
  const int scratch_shift = team_rank * tile_size;

  // extract and wrap
  const WignerWrapper<real_type, vector_length> ulist_wrapper((complex*)team.team_shmem().get_shmem(team.team_size() * tile_size * sizeof(complex), 0) + scratch_shift, iatom_mod);

  // load parameters
  const complex a = a_gpu(iatom, jnbor);
  const complex b = b_gpu(iatom, jnbor);
  const real_type sfac = sfac_gpu(iatom, jnbor, 0);

  const int jelem = element(iatom, jnbor);

  // we need to "choose" when to bend
  #ifdef LMP_KK_DEVICE_COMPILE
  #pragma unroll
  #endif
  for (int j_bend = 0; j_bend <= twojmax; j_bend++) {
    evaluate_ui_jbend(ulist_wrapper, a, b, sfac, jelem, iatom, j_bend);
  }
}

// Core "evaluation" kernel that gets reused in `compute_ui_small` and `compute_ui_large`
template<class DeviceType, typename real_type, int vector_length>
KOKKOS_FORCEINLINE_FUNCTION
void SNAKokkos<DeviceType, real_type, vector_length>::evaluate_ui_jbend(const WignerWrapper<real_type, vector_length>& ulist_wrapper,
          const complex& a, const complex& b, const real_type& sfac, const int& jelem,
          const int& iatom, const int& j_bend) const
{
  // utot(j,ma,mb) = 0 for all j,ma,ma
  // utot(j,ma,ma) = 1 for all j,ma
  // for j in neighbors of i:
  //   compute r0 = (x,y,z,z0)
  //   utot(j,ma,mb) += u(r0;j,ma,mb) for all j,ma,mb

  // level 0 is just 1.
  ulist_wrapper.set(0, complex::one());

  // j from before the bend, don't store, mb == 0
  for (int j = 1; j <= j_bend; j++) {

    constexpr int mb = 0; // intentional for readability, compiler should optimize this out

    complex ulist_accum = complex::zero();

    int ma;
    for (ma = 0; ma < j; ma++) {

      // grab the cached value
      const complex ulist_prev = ulist_wrapper.get(ma);

      // ulist_accum += rootpq * a.conj() * ulist_prev;
      real_type rootpq = rootpqarray(j - ma, j - mb);
      ulist_accum.re += rootpq * (a.re * ulist_prev.re + a.im * ulist_prev.im);
      ulist_accum.im += rootpq * (a.re * ulist_prev.im - a.im * ulist_prev.re);

      // store ulist_accum, we atomic accumulate values after the bend, so no atomic add here
      ulist_wrapper.set(ma, ulist_accum);

      // next value
      // ulist_accum = -rootpq * b.conj() * ulist_prev;
      rootpq = rootpqarray(ma + 1, j - mb);
      ulist_accum.re = -rootpq * (b.re * ulist_prev.re + b.im * ulist_prev.im);
      ulist_accum.im = -rootpq * (b.re * ulist_prev.im - b.im * ulist_prev.re);

    }

    ulist_wrapper.set(ma, ulist_accum);
  }

  // now we're after the bend, start storing but only up to the "half way point"
  const int j_half_way = MIN(2 * j_bend, twojmax);

  int mb = 1;
  int j; //= j_bend + 1; // need this value below
  for (j = j_bend + 1; j <= j_half_way; j++) {

    const int jjup = idxu_half_block[j-1] + (mb - 1) * j;

    complex ulist_accum = complex::zero();

    int ma;
    for (ma = 0; ma < j; ma++) {

      // grab the cached value
      const complex ulist_prev = ulist_wrapper.get(ma);

      // atomic add the previous level here
      Kokkos::atomic_add(&(ulisttot_re(iatom, jelem, jjup + ma)), ulist_prev.re * sfac);
      Kokkos::atomic_add(&(ulisttot_im(iatom, jelem, jjup + ma)), ulist_prev.im * sfac);

      // ulist_accum += rootpq * b * ulist_prev;
      real_type rootpq = rootpqarray(j - ma, mb);
      ulist_accum.re += rootpq * (b.re * ulist_prev.re - b.im * ulist_prev.im);
      ulist_accum.im += rootpq * (b.re * ulist_prev.im + b.im * ulist_prev.re);

      // store ulist_accum
      ulist_wrapper.set(ma, ulist_accum);

      // next value
      // ulist_accum = rootpq * a * ulist_prev;
      rootpq = rootpqarray(ma + 1, mb);
      ulist_accum.re = rootpq * (a.re * ulist_prev.re - a.im * ulist_prev.im);
      ulist_accum.im = rootpq * (a.re * ulist_prev.im + a.im * ulist_prev.re);
    }

    ulist_wrapper.set(ma, ulist_accum);
    mb++;
  }

  // atomic add the last level
  const int jjup = idxu_half_block[j-1] + (mb - 1) * j;

  for (int ma = 0; ma < j; ma++) {
    const complex ulist_prev = ulist_wrapper.get(ma);

    // atomic add the previous level here
    Kokkos::atomic_add(&(ulisttot_re(iatom, jelem, jjup + ma)), ulist_prev.re * sfac);
    Kokkos::atomic_add(&(ulisttot_im(iatom, jelem, jjup + ma)), ulist_prev.im * sfac);
  }
}

/* ----------------------------------------------------------------------
   compute Ui by summing over bispectrum components. CPU only.
   This first computes Wigner U-functions for one neighbor.
   `ulisttot` uses a "cached" data layout, matching the amount of
   information stored between layers via scratch memory on the GPU path.
   Next, it adds Wigner U-functions for each neighbor to ulisttot which is
   in a "half" data layout, which is a compressed layout
   which still keeps the recursive calculation simple.
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
template <bool need_atomics> KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType, real_type, vector_length>::compute_ui_cpu(const int& iatom, const int& jnbor) const
{
  // utot(j,ma,mb) = 0 for all j,ma,ma
  // utot(j,ma,ma) = 1 for all j,ma
  // for j in neighbors of i:
  //   compute r0 = (x,y,z,z0)
  //   utot(j,ma,mb) += u(r0;j,ma,mb) for all j,ma,mb

  const real_type x = rij(iatom,jnbor,0);
  const real_type y = rij(iatom,jnbor,1);
  const real_type z = rij(iatom,jnbor,2);
  const real_type rsq = x * x + y * y + z * z;
  const real_type r = sqrt(rsq);

  const real_type theta0 = (r - rmin0) * rfac0 * MY_PI / (rcutij(iatom,jnbor) - rmin0);
  //    theta0 = (r - rmin0) * rscale0;
  const real_type z0 = r / tan(theta0);

  // begin what was "compute_uarray_cpu"

  // compute Cayley-Klein parameters for unit quaternion
  real_type r0inv = static_cast<real_type>(1.0) / sqrt(r * r + z0 * z0);
  complex a = { r0inv * z0, -r0inv * z };
  complex b = { r0inv * y, -r0inv * x };

  // VMK Section 4.8.2
  ulist_cpu(iatom, jnbor, 0) = complex::one();

  for (int j = 1; j <= twojmax; j++) {
    int jju = idxu_cache_block[j]; // removed "const" to work around GCC 7 bug
    int jjup = idxu_cache_block[j-1]; // removed "const" to work around GCC 7 bug

    // fill in left side of matrix layer from previous layer
    for (int mb = 0; 2*mb <= j; mb++) {
      int jju_index = jju + (j + 1) * mb;
      int jjup_index = jjup + j * mb;

      complex ui = complex::zero();

      for (int ma = 0; ma < j; ma++) {
        complex ui_prev = ulist_cpu(iatom, jnbor, jjup_index);

        real_type rootpq = rootpqarray(j - ma, j - mb);
        ui.re += rootpq * (a.re * ui_prev.re + a.im * ui_prev.im);
        ui.im += rootpq * (a.re * ui_prev.im - a.im * ui_prev.re);

        ulist_cpu(iatom, jnbor, jju_index) = ui;

        rootpq = rootpqarray(ma + 1, j - mb);
        ui.re = -rootpq * (b.re * ui_prev.re + b.im * ui_prev.im);
        ui.im = -rootpq * (b.re * ui_prev.im - b.im * ui_prev.re);

        jju_index++;
        jjup_index++;
      }

      ulist_cpu(iatom, jnbor, jju_index) = ui;
    }

    // If j is odd (half-integer in the mathematical convention), we need
    // to add one more row for convenience (for now). This can either be done
    // via symmetry (see the commented code below), or by the equations to fill
    // from the left instead of the right
    if (j % 2 == 1) {
      int mb = j / 2;
      // begin filling in the extra row
      int jju_index = jju + (mb + 1) * (j + 1);
      int jjup_index = jjup + mb * j;

      complex ui = complex::zero();

      for (int ma = 0; ma < j; ma++) {
        complex ui_prev = ulist_cpu(iatom, jnbor, jjup_index);

        real_type rootpq = rootpqarray(j - ma, mb + 1);
        ui.re += rootpq * (b.re * ui_prev.re - b.im * ui_prev.im);
        ui.im += rootpq * (b.re * ui_prev.im + b.im * ui_prev.re);

        ulist_cpu(iatom, jnbor, jju_index) = ui;

        rootpq = rootpqarray(ma + 1, mb + 1);
        ui.re = rootpq * (a.re * ui_prev.re - a.im * ui_prev.im);
        ui.im = rootpq * (a.re * ui_prev.im + a.im * ui_prev.re);

        jju_index++;
        jjup_index++;
      }

      ulist_cpu(iatom, jnbor, jju_index) = ui;
    }
  }

  // begin what was add_uarraytot

  const real_type sfac = compute_sfac(r, rcutij(iatom,jnbor), sinnerij(iatom,jnbor), dinnerij(iatom,jnbor)) * wj(iatom,jnbor);
  const auto jelem = element(iatom, jnbor);

  for (int j = 0; j <= twojmax; j++) {
    int jju_half = idxu_half_block[j]; // index into ulisttot
    int jju_cache = idxu_cache_block[j]; // index into ulist

    int count = 0;
    for (int mb = 0; 2*mb <= j; mb++) {
      for (int ma = 0; ma <= j; ma++) {
        if constexpr (need_atomics) {
          Kokkos::atomic_add(&(ulisttot_re(iatom, jelem, jju_half+count)), sfac * ulist_cpu(iatom, jnbor, jju_cache+count).re);
          Kokkos::atomic_add(&(ulisttot_im(iatom, jelem, jju_half+count)), sfac * ulist_cpu(iatom, jnbor, jju_cache+count).im);
        } else {
          ulisttot_re(iatom, jelem, jju_half+count) += sfac * ulist_cpu(iatom, jnbor, jju_cache+count).re;
          ulisttot_im(iatom, jelem, jju_half+count) += sfac * ulist_cpu(iatom, jnbor, jju_cache+count).im;
        }
        count++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
  De-symmetrize ulisttot_re and _im and pack it into a unified ulisttot
  structure, fused in with zeroing ylist
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType, real_type, vector_length>::transform_ui(const int& iatom, const int& idxu) const
{
  int elem_count = chem_flag ? nelements : 1;

  for (int ielem = 0; ielem < elem_count; ielem++) {

    const FullHalfMapper mapper = idxu_full_half[idxu];

    auto utot_re = ulisttot_re(iatom, ielem, mapper.idxu_half);
    auto utot_im = ulisttot_im(iatom, ielem, mapper.idxu_half);

    if (mapper.flip_sign == 1) {
      utot_im = -utot_im;
    } else if (mapper.flip_sign == -1) {
      utot_re = -utot_re;
    }

    ulisttot(iatom, ielem, idxu) = { utot_re, utot_im };

    if (mapper.flip_sign == 0) {
      ylist_re(iatom, ielem, mapper.idxu_half) = 0.;
      ylist_im(iatom, ielem, mapper.idxu_half) = 0.;
    }
  }
}

/* ----------------------------------------------------------------------
   compute Zi by summing over products of Ui
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
template <bool chemsnap> KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType, real_type, vector_length>::compute_zi(const int& iatom, const int& jjz) const
{
  int j1, j2, j, ma1min, ma2max, mb1min, mb2max, na, nb, idxcg;
  idxz(jjz).get_zi(j1, j2, j, ma1min, ma2max, mb1min, mb2max, na, nb, idxcg);

  const real_type *cgblock = cglist.data() + idxcg;

  if constexpr (chemsnap) {
    int idouble = 0;
    for (int elem1 = 0; elem1 < nelements; elem1++) {
      for (int elem2 = 0; elem2 < nelements; elem2++) {
        zlist(iatom, idouble, jjz) = evaluate_zi(j1, j2, j, ma1min, ma2max, mb1min, mb2max, na, nb, iatom, elem1, elem2, cgblock);
        idouble++;
      } // end loop over elem2
    } // end loop over elem1
  } else {
    zlist(iatom, 0, jjz) = evaluate_zi(j1, j2, j, ma1min, ma2max, mb1min, mb2max, na, nb, iatom, 0, 0, cgblock);
  }
}

/* ----------------------------------------------------------------------
   Core "evaluation" kernel that computes a single zlist value.
   This gets used in both `compute_zi` and `compute_yi`
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_FORCEINLINE_FUNCTION
typename SNAKokkos<DeviceType, real_type, vector_length>::complex SNAKokkos<DeviceType, real_type, vector_length>::evaluate_zi(const int& j1, const int& j2, const int& j,
        const int& ma1min, const int& ma2max, const int& mb1min, const int& mb2max, const int& na, const int& nb,
        const int& iatom, const int& elem1, const int& elem2, const real_type* cgblock) const {
  complex ztmp = complex::zero();

  int jju1 = idxu_block[j1] + (j1+1)*mb1min;
  int jju2 = idxu_block[j2] + (j2+1)*mb2max;
  int icgb = mb1min*(j2+1) + mb2max;
  #ifdef LMP_KK_DEVICE_COMPILE
  #pragma unroll
  #endif
  for (int ib = 0; ib < nb; ib++) {

    int ma1 = ma1min;
    int ma2 = ma2max;
    int icga = ma1min*(j2+1) + ma2max;

    #ifdef LMP_KK_DEVICE_COMPILE
    #pragma unroll
    #endif
    for (int ia = 0; ia < na; ia++) {
      const complex utot1 = ulisttot(iatom, elem1, jju1+ma1);
      const complex utot2 = ulisttot(iatom, elem2, jju2+ma2);
      const real_type cgcoeff_a = cgblock[icga];
      const real_type cgcoeff_b = cgblock[icgb];
      ztmp.re += cgcoeff_a * cgcoeff_b * (utot1.re * utot2.re - utot1.im * utot2.im);
      ztmp.im += cgcoeff_a * cgcoeff_b * (utot1.re * utot2.im + utot1.im * utot2.re);
      ma1++;
      ma2--;
      icga += j2;
    } // end loop over ia

    jju1 += j1 + 1;
    jju2 -= j2 + 1;
    icgb += j2;
  } // end loop over ib

  if (bnorm_flag) {
    const real_type scale = static_cast<real_type>(1) / static_cast<real_type>(j + 1);
    ztmp.re *= scale;
    ztmp.im *= scale;
  }

  return ztmp;
}

/* ----------------------------------------------------------------------
   compute Bi by summing conj(Ui)*Zi
   AoSoA data layout to take advantage of coalescing, avoiding warp
   divergence.
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
template <bool chemsnap> KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType, real_type, vector_length>::compute_bi(const int& iatom, const int& jjb) const
{
  // for j1 = 0,...,twojmax
  //   for j2 = 0,twojmax
  //     for j = |j1-j2|,Min(twojmax,j1+j2),2
  //        b(j1,j2,j) = 0
  //        for mb = 0,...,jmid
  //          for ma = 0,...,j
  //            b(j1,j2,j) +=
  //              2*Conj(u(j,ma,mb))*z(j1,j2,j,ma,mb)

  const int j1 = idxb(jjb,0);
  const int j2 = idxb(jjb,1);
  const int j = idxb(jjb,2);

  const int jjz = idxz_block(j1,j2,j);
  const int jju = idxu_block[j];

  if constexpr (chemsnap) {
    int itriple = 0;
    for (int elem1 = 0; elem1 < nelements; elem1++) {
      for (int elem2 = 0; elem2 < nelements; elem2++) {
        for (int elem3 = 0; elem3 < nelements; elem3++) {
          blist(iatom, itriple, jjb) = evaluate_bi(j, jjz, jju, iatom, elem1, elem2, elem3);
          itriple++;
        } // end loop over elem3
      } // end loop over elem2
    } // end loop over elem1
  } else {
    blist(iatom, 0, jjb) = evaluate_bi(j, jjz, jju, iatom, 0, 0, 0);
  }
}

/* ----------------------------------------------------------------------
   Core "evaluation" kernel that computes a single blist value.
   This gets used in `compute_bi`
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
real_type SNAKokkos<DeviceType, real_type, vector_length>::evaluate_bi(const int& j, const int& jjz, const int& jju, const int& iatom, const int& elem1, const int& elem2, const int& elem3) const
{
  // this computes the:
  //        b(j1,j2,j) = 0
  //        for mb = 0,...,jmid
  //          for ma = 0,...,j
  //            b(j1,j2,j) +=
  //              2*Conj(u(j,ma,mb))*z(j1,j2,j,ma,mb)
  // portion

  const int idouble = elem1 * nelements + elem2;
  real_type sumzu = 0.0;
  real_type sumzu_temp = 0.0;

  for (int mb = 0; 2*mb < j; mb++) {
    for (int ma = 0; ma <= j; ma++) {
      const int jju_index = jju+mb*(j+1)+ma;
      const int jjz_index = jjz+mb*(j+1)+ma;
      if (2*mb == j) return 0; // I think we can remove this?
      const complex utot = ulisttot(iatom, elem3, jju_index);
      const complex zloc = zlist(iatom, idouble, jjz_index);
      sumzu_temp += utot.re * zloc.re + utot.im * zloc.im;
    }
  }
  sumzu += sumzu_temp;

  // For j even, special treatment for middle column
  if (j%2 == 0) {
    sumzu_temp = 0.;

    const int mb = j/2;
    for (int ma = 0; ma < mb; ma++) {
      const int jju_index = jju+(mb-1)*(j+1)+(j+1)+ma;
      const int jjz_index = jjz+(mb-1)*(j+1)+(j+1)+ma;

      const complex utot = ulisttot(iatom, elem3, jju_index);
      const complex zloc = zlist(iatom, idouble, jjz_index);
      sumzu_temp += utot.re * zloc.re + utot.im * zloc.im;

    }
    sumzu += sumzu_temp;

    const int ma = mb;
    const int jju_index = jju+(mb-1)*(j+1)+(j+1)+ma;
    const int jjz_index = jjz+(mb-1)*(j+1)+(j+1)+ma;

    const complex utot = ulisttot(iatom, elem3, jju_index);
    const complex zloc = zlist(iatom, idouble, jjz_index);
    sumzu += static_cast<real_type>(0.5) * (utot.re * zloc.re + utot.im * zloc.im);
  } // end if jeven

  sumzu *= static_cast<real_type>(2.0);
  if (bzero_flag) {
    if (!wselfall_flag) {
      if (elem1 == elem2 && elem1 == elem3) {
        sumzu -= bzero[j];
      }
    } else {
      sumzu -= bzero[j];
    }
  }
  return sumzu;
      //} // end loop over j
    //} // end loop over j1, j2
}


/* ----------------------------------------------------------------------
   compute beta by either appropriately copying it from d_coeffi
   or accumulating the quadratic terms from blist
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType, real_type, vector_length>::compute_beta_linear(const int& iatom, const int& idxb, const int& ielem) const
{
  auto d_coeffi = Kokkos::subview(d_coeffelem, ielem, Kokkos::ALL);

  for (int itriple = 0; itriple < ntriples; itriple++) {
    int icoeff = idxb + itriple * idxb_max;
    d_beta(iatom, icoeff) = d_coeffi[icoeff+1];
  }
}

template<class DeviceType, typename real_type, int vector_length>
template <bool need_atomics>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType, real_type, vector_length>::compute_beta_quadratic(const int& iatom, const int& idxb, const int& ielem) const
{
  auto d_coeffi = Kokkos::subview(d_coeffelem, ielem, Kokkos::ALL);

  // handle quadratic && chemflag as a special case
  if (chem_flag) {
    if (idxb == 0) {
      // no need to use atomics, we're just serializing
      int k = ncoeff+1;
      for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
        const auto idxb = icoeff % idxb_max;
        const auto idx_chem = icoeff / idxb_max;
        real_type bveci = blist(iatom, idx_chem, idxb);
        d_beta(iatom, icoeff) += d_coeffi[k] * bveci;
        k++;
        for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
          const auto jdxb = jcoeff % idxb_max;
          const auto jdx_chem = jcoeff / idxb_max;
          real_type bvecj = blist(iatom, jdx_chem, jdxb);
          d_beta(iatom, icoeff) += d_coeffi[k] * bvecj;
          d_beta(iatom, jcoeff) += d_coeffi[k] * bveci;
          k++;
        }
      }
    }
  } else {
    // Compute triangular partial sum via a closed form to get the starting offset
    int k = (idxb * (1 + 2 * idxb_max - idxb)) / 2 + idxb_max + 1;
    real_type bveci = blist(iatom, 0, idxb);

    // Locally accumulate the contribution to d_beta(iatom, idxb)
    real_type beta_idxb_accum = d_coeffi[k] * bveci;
    k++;

    for (int jdxb = idxb + 1; jdxb < idxb_max; jdxb++) {
      real_type bvecj = blist(iatom, 0, jdxb);
      real_type coeff_k = d_coeffi[k];
      beta_idxb_accum += coeff_k * bvecj;

      // Accumulate "half" contribution into d_beta(iatom, jdxb)
      if constexpr (need_atomics)
        Kokkos::atomic_add(&d_beta(iatom, jdxb), coeff_k * bveci);
      else
        d_beta(iatom, jdxb) += coeff_k * bveci;

      k++;
    }

    if constexpr (need_atomics)
      Kokkos::atomic_add(&d_beta(iatom, idxb), beta_idxb_accum);
    else
      d_beta(iatom, idxb) += beta_idxb_accum;
  }
}


/* ----------------------------------------------------------------------
   Compute Yi from Ui without storing Zi, looping over zlist indices.
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
template <bool chemsnap, bool need_atomics>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType, real_type, vector_length>::compute_yi(const int& iatom, const int& jjz) const
{
  int j1, j2, j, ma1min, ma2max, mb1min, mb2max, na, nb, jju_half, idxcg;
  idxz(jjz).get_yi(j1, j2, j, ma1min, ma2max, mb1min, mb2max, na, nb, jju_half, idxcg);

  const real_type *cgblock = cglist.data() + idxcg;
  //int mb = (2 * (mb1min+mb2max) - j1 - j2 + j) / 2;
  //int ma = (2 * (ma1min+ma2max) - j1 - j2 + j) / 2;

  if constexpr (chemsnap) {
    for (int elem1 = 0; elem1 < nelements; elem1++) {
      for (int elem2 = 0; elem2 < nelements; elem2++) {

        const complex ztmp = evaluate_zi(j1, j2, j, ma1min, ma2max, mb1min, mb2max, na, nb, iatom, elem1, elem2, cgblock);

        // apply to z(j1,j2,j,ma,mb) to unique element of y(j)
        // find right y_list[jju] and beta(iatom,jjb) entries
        // multiply and divide by j+1 factors
        // account for multiplicity of 1, 2, or 3

        // pick out right beta value
        for (int elem3 = 0; elem3 < nelements; elem3++) {

          const real_type betaj = evaluate_beta_scaled<true>(j1, j2, j, iatom, elem1, elem2, elem3);

          if constexpr (need_atomics) {
            Kokkos::atomic_add(&(ylist_re(iatom, elem3, jju_half)), betaj * ztmp.re);
            Kokkos::atomic_add(&(ylist_im(iatom, elem3, jju_half)), betaj * ztmp.im);
          } else {
            ylist_re(iatom, elem3, jju_half) += betaj * ztmp.re;
            ylist_im(iatom, elem3, jju_half) += betaj * ztmp.im;
          }
        } // end loop over elem3
      } // end loop over elem2
    } // end loop over elem1
  } else {
    const complex ztmp = evaluate_zi(j1, j2, j, ma1min, ma2max, mb1min, mb2max, na, nb, iatom, 0, 0, cgblock);
    const real_type betaj = evaluate_beta_scaled<false>(j1, j2, j, iatom, 0, 0, 0);

    if constexpr (need_atomics) {
      Kokkos::atomic_add(&(ylist_re(iatom, 0, jju_half)), betaj * ztmp.re);
      Kokkos::atomic_add(&(ylist_im(iatom, 0, jju_half)), betaj * ztmp.im);
    } else {
      ylist_re(iatom, 0, jju_half) += betaj * ztmp.re;
      ylist_im(iatom, 0, jju_half) += betaj * ztmp.im;
    }
  }
}

/* ----------------------------------------------------------------------
   compute Yi from Ui with the precomputed Zi.
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
template <bool chemsnap, bool need_atomics>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType, real_type, vector_length>::compute_yi_with_zlist(const int& iatom, const int& jjz) const
{
  int j1, j2, j, jju_half;
  idxz(jjz).get_yi_with_zlist(j1, j2, j, jju_half);

  if constexpr (chemsnap) {
    int idouble = 0;
    for (int elem1 = 0; elem1 < nelements; elem1++) {
      for (int elem2 = 0; elem2 < nelements; elem2++) {
        const complex ztmp = zlist(iatom, idouble, jjz);
        // apply to z(j1,j2,j,ma,mb) to unique element of y(j)
        // find right y_list[jju] and beta(iatom,jjb) entries
        // multiply and divide by j+1 factors
        // account for multiplicity of 1, 2, or 3
        // pick out right beta value
        for (int elem3 = 0; elem3 < nelements; elem3++) {

          const real_type betaj = evaluate_beta_scaled<true>(j1, j2, j, iatom, elem1, elem2, elem3);

          if constexpr (need_atomics) {
            Kokkos::atomic_add(&(ylist_re(iatom, elem3, jju_half)), betaj * ztmp.re);
            Kokkos::atomic_add(&(ylist_im(iatom, elem3, jju_half)), betaj * ztmp.im);
          } else {
            ylist_re(iatom, elem3, jju_half) += betaj * ztmp.re;
            ylist_im(iatom, elem3, jju_half) += betaj * ztmp.im;
          }
        } // end loop over elem3
        idouble++;
      } // end loop over elem2
    } // end loop over elem1
  } else {
    const complex ztmp = zlist(iatom, 0, jjz);
    const real_type betaj = evaluate_beta_scaled<false>(j1, j2, j, iatom, 0, 0, 0);

    if constexpr (need_atomics) {
      Kokkos::atomic_add(&(ylist_re(iatom, 0, jju_half)), betaj * ztmp.re);
      Kokkos::atomic_add(&(ylist_im(iatom, 0, jju_half)), betaj * ztmp.im);
    } else {
      ylist_re(iatom, 0, jju_half) += betaj * ztmp.re;
      ylist_im(iatom, 0, jju_half) += betaj * ztmp.im;
    }
  }
}

/* ----------------------------------------------------------------------
   Core "evaluation" kernel that extracts and rescales the appropriate
   `beta` value which gets used in both `compute_yi` and `compute_yi_from_zlist`
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
template <bool chemsnap> KOKKOS_FORCEINLINE_FUNCTION
typename SNAKokkos<DeviceType, real_type, vector_length>::real_type SNAKokkos<DeviceType, real_type, vector_length>::evaluate_beta_scaled(const int& j1, const int& j2, const int& j,
          const int& iatom, const int& elem1, const int& elem2, const int& elem3) const {

  int itriple_jjb = 0;
  real_type factor = 0;

  if constexpr (chemsnap) {
    if (j >= j1) {
      itriple_jjb = ((elem1 * nelements + elem2) * nelements + elem3) * idxb_max + idxb_block(j1, j2, j);
      if (j1 == j) {
        if (j2 == j) factor = 3;
        else factor = 2;
      } else factor = 1;
    } else if (j >= j2) {
      itriple_jjb = ((elem3 * nelements + elem2) * nelements + elem1) * idxb_max + idxb_block(j, j2, j1);
      if (j2 == j) factor = 2;
      else factor = 1;
    } else {
      itriple_jjb = ((elem2 * nelements + elem3) * nelements + elem1) * idxb_max + idxb_block(j2, j, j1);
      factor = 1;
    }
  } else {
    if (j >= j1) {
      itriple_jjb = idxb_block(j1, j2, j);
      if (j1 == j) {
        if (j2 == j) factor = 3;
        else factor = 2;
      } else factor = 1;
    } else if (j >= j2) {
      itriple_jjb = idxb_block(j, j2, j1);
      if (j2 == j) factor = 2;
      else factor = 1;
    } else {
      itriple_jjb = idxb_block(j2, j, j1);
      factor = 1;
    }
  }

  real_type betaj = factor * d_beta(iatom, itriple_jjb);

  if (!bnorm_flag && j1 > j) {
    const real_type scale = static_cast<real_type>(j1 + 1) / static_cast<real_type>(j + 1);
    betaj *= scale;
  }

  return betaj;

}

/* ----------------------------------------------------------------------
   Fused calculation of the derivative of Ui w.r.t. atom j
   and accumulation into dEidRj. GPU only.
------------------------------------------------------------------------- */

// Version of the code that exposes additional parallelism by threading over `j_bend` values
template<class DeviceType, typename real_type, int vector_length>
template<int dir>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType, real_type, vector_length>::compute_fused_deidrj_small(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, const int iatom_mod, const int j_bend, const int jnbor, const int iatom_div) const
{
  const int iatom = iatom_mod + vector_length * iatom_div;
  // get shared memory offset
  // scratch size: 32 atoms * (twojmax+1) cached values, no double buffer
  const int tile_size = vector_length * (twojmax + 1);

  const int team_rank = team.team_rank();
  const int scratch_shift = team_rank * tile_size;

  // extract, wrap shared memory buffer
  WignerWrapper<real_type, vector_length> ulist_wrapper((complex*)team.team_shmem().get_shmem(team.team_size() * tile_size * sizeof(complex), 0) + scratch_shift, iatom_mod);
  WignerWrapper<real_type, vector_length> dulist_wrapper((complex*)team.team_shmem().get_shmem(team.team_size() * tile_size * sizeof(complex), 0) + scratch_shift, iatom_mod);

  // load parameters
  const complex a = a_gpu(iatom, jnbor);
  const complex b = b_gpu(iatom, jnbor);
  const complex da = da_gpu(iatom, jnbor, dir);
  const complex db = db_gpu(iatom, jnbor, dir);
  const real_type sfac = sfac_gpu(iatom, jnbor, 0);
  const real_type dsfacu = sfac_gpu(iatom, jnbor, dir + 1); // dsfac * u

  const int jelem = element(iatom, jnbor);

  // compute the contribution to dedr_full_sum for one "bend" location
  const real_type dedr_full_sum = evaluate_duidrj_jbend(ulist_wrapper, a, b, sfac, dulist_wrapper, da, db, dsfacu,
                                                       jelem, iatom, j_bend);

  // dedr gets zeroed out at the start of each iteration in compute_cayley_klein
  Kokkos::atomic_add(&(dedr(iatom, jnbor, dir)), static_cast<real_type>(2.0) * dedr_full_sum);

}

// Version of the code that loops over all `j_bend` values which reduces integer arithmetic
// and some amount of load imbalance, at the expense of reducing parallelism
template<class DeviceType, typename real_type, int vector_length>
template<int dir>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType, real_type, vector_length>::compute_fused_deidrj_large(const typename Kokkos::TeamPolicy<DeviceType>::member_type& team, const int iatom_mod, const int jnbor, const int iatom_div) const
{
  const int iatom = iatom_mod + vector_length * iatom_div;
  // get shared memory offset
  // scratch size: 32 atoms * (twojmax+1) cached values, no double buffer
  const int tile_size = vector_length * (twojmax + 1);

  const int team_rank = team.team_rank();
  const int scratch_shift = team_rank * tile_size;

  // extract, wrap shared memory buffer
  WignerWrapper<real_type, vector_length> ulist_wrapper((complex*)team.team_shmem().get_shmem(team.team_size() * tile_size * sizeof(complex), 0) + scratch_shift, iatom_mod);
  WignerWrapper<real_type, vector_length> dulist_wrapper((complex*)team.team_shmem().get_shmem(team.team_size() * tile_size * sizeof(complex), 0) + scratch_shift, iatom_mod);

  // load parameters
  const complex a = a_gpu(iatom, jnbor);
  const complex b = b_gpu(iatom, jnbor);
  const complex da = da_gpu(iatom, jnbor, dir);
  const complex db = db_gpu(iatom, jnbor, dir);
  const real_type sfac = sfac_gpu(iatom, jnbor, 0);
  const real_type dsfacu = sfac_gpu(iatom, jnbor, dir + 1); // dsfac * u

  const int jelem = element(iatom, jnbor);

  // compute the contributions to dedr_full_sum for all "bend" locations
  real_type dedr_full_sum = static_cast<real_type>(0);
  #ifdef LMP_KK_DEVICE_COMPILE
  #pragma unroll
  #endif
  for (int j_bend = 0; j_bend <= twojmax; j_bend++) {
    dedr_full_sum += evaluate_duidrj_jbend(ulist_wrapper, a, b, sfac, dulist_wrapper, da, db, dsfacu,
                                          jelem, iatom,j_bend);
  }

  // there's one thread per atom, neighbor pair, so no need to make this atomic
  dedr(iatom, jnbor, dir) = static_cast<real_type>(2.0) * dedr_full_sum;

}

// Core "evaluation" kernel that gets reused in `compute_fused_deidrj_small` and
// `compute_fused_deidrj_large`
template<class DeviceType, typename real_type, int vector_length>
KOKKOS_FORCEINLINE_FUNCTION
typename SNAKokkos<DeviceType, real_type, vector_length>::real_type SNAKokkos<DeviceType, real_type, vector_length>::evaluate_duidrj_jbend(const WignerWrapper<real_type, vector_length>& ulist_wrapper, const complex& a, const complex& b, const real_type& sfac,
                      const WignerWrapper<real_type, vector_length>& dulist_wrapper, const complex& da, const complex& db, const real_type& dsfacu,
                      const int& jelem, const int& iatom, const int& j_bend) const {

  real_type dedr_full_sum = static_cast<real_type>(0);

  // level 0 is just 1, 0
  ulist_wrapper.set(0, complex::one());
  dulist_wrapper.set(0, complex::zero());

  // j from before the bend, don't store, mb == 0
  // this is "creeping up the side"
  for (int j = 1; j <= j_bend; j++) {

    constexpr int mb = 0; // intentional for readability, compiler should optimize this out

    complex ulist_accum = complex::zero();
    complex dulist_accum = complex::zero();

    int ma;
    for (ma = 0; ma < j; ma++) {

      // grab the cached value
      const complex ulist_prev = ulist_wrapper.get(ma);
      const complex dulist_prev = dulist_wrapper.get(ma);

      // ulist_accum += rootpq * a.conj() * ulist_prev;
      real_type rootpq = rootpqarray(j - ma, j - mb);
      ulist_accum.re += rootpq * (a.re * ulist_prev.re + a.im * ulist_prev.im);
      ulist_accum.im += rootpq * (a.re * ulist_prev.im - a.im * ulist_prev.re);

      // product rule of above
      dulist_accum.re += rootpq * (da.re * ulist_prev.re + da.im * ulist_prev.im + a.re * dulist_prev.re + a.im * dulist_prev.im);
      dulist_accum.im += rootpq * (da.re * ulist_prev.im - da.im * ulist_prev.re + a.re * dulist_prev.im - a.im * dulist_prev.re);

      // store ulist_accum, we atomic accumulate values after the bend, so no atomic add here
      ulist_wrapper.set(ma, ulist_accum);
      dulist_wrapper.set(ma, dulist_accum);

      // next value
      // ulist_accum = -rootpq * b.conj() * ulist_prev;
      rootpq = rootpqarray(ma + 1, j - mb);
      ulist_accum.re = -rootpq * (b.re * ulist_prev.re + b.im * ulist_prev.im);
      ulist_accum.im = -rootpq * (b.re * ulist_prev.im - b.im * ulist_prev.re);

      // product rule of above
      dulist_accum.re = -rootpq * (db.re * ulist_prev.re + db.im * ulist_prev.im + b.re * dulist_prev.re + b.im * dulist_prev.im);
      dulist_accum.im = -rootpq * (db.re * ulist_prev.im - db.im * ulist_prev.re + b.re * dulist_prev.im - b.im * dulist_prev.re);

    }

    ulist_wrapper.set(ma, ulist_accum);
    dulist_wrapper.set(ma, dulist_accum);
  }

  // now we're after the bend, start storing but only up to the "half way point"
  const int j_half_way = MIN(2 * j_bend, twojmax);

  int mb = 1;
  int j; //= j_bend + 1; // need this value below
  for (j = j_bend + 1; j <= j_half_way; j++) {

    const int jjup = idxu_half_block[j-1] + (mb - 1) * j;

    complex ulist_accum = complex::zero();
    complex dulist_accum = complex::zero();

    int ma;
    for (ma = 0; ma < j; ma++) {

      // grab y_local early
      // this will never be the last element of a row, no need to rescale.
      complex y_local = complex(ylist_re(iatom, jelem, jjup + ma), ylist_im(iatom, jelem, jjup+ma));

      // grab the cached value
      const complex ulist_prev = ulist_wrapper.get(ma);
      const complex dulist_prev = dulist_wrapper.get(ma);

      // ulist_accum += rootpq * b * ulist_prev;
      real_type rootpq = rootpqarray(j - ma, mb);
      ulist_accum.re += rootpq * (b.re * ulist_prev.re - b.im * ulist_prev.im);
      ulist_accum.im += rootpq * (b.re * ulist_prev.im + b.im * ulist_prev.re);

      // product rule of above
      dulist_accum.re += rootpq * (db.re * ulist_prev.re - db.im * ulist_prev.im + b.re * dulist_prev.re - b.im * dulist_prev.im);
      dulist_accum.im += rootpq * (db.re * ulist_prev.im + db.im * ulist_prev.re + b.re * dulist_prev.im + b.im * dulist_prev.re);

      // store ulist_accum
      ulist_wrapper.set(ma, ulist_accum);
      dulist_wrapper.set(ma, dulist_accum);

      // Directly accumulate deidrj into sum_tmp
      const complex du_prod = (dsfacu * ulist_prev) + (sfac * dulist_prev);
      dedr_full_sum += du_prod.re * y_local.re + du_prod.im * y_local.im;

      // next value
      // ulist_accum = rootpq * a * ulist_prev;
      rootpq = rootpqarray(ma + 1, mb);
      ulist_accum.re = rootpq * (a.re * ulist_prev.re - a.im * ulist_prev.im);
      ulist_accum.im = rootpq * (a.re * ulist_prev.im + a.im * ulist_prev.re);

      // product rule of above
      dulist_accum.re = rootpq * (da.re * ulist_prev.re - da.im * ulist_prev.im + a.re * dulist_prev.re - a.im * dulist_prev.im);
      dulist_accum.im = rootpq * (da.re * ulist_prev.im + da.im * ulist_prev.re + a.re * dulist_prev.im + a.im * dulist_prev.re);

    }

    ulist_wrapper.set(ma, ulist_accum);
    dulist_wrapper.set(ma, dulist_accum);

    mb++;
  }

  // accumulate the last level
  const int jjup = idxu_half_block[j-1] + (mb - 1) * j;

  for (int ma = 0; ma < j; ma++) {
    // grab y_local early
    complex y_local = complex(ylist_re(iatom, jelem, jjup + ma), ylist_im(iatom, jelem, jjup+ma));
    if (j % 2 == 1 && 2*(mb-1) == j-1) { // double check me...
      if (ma == (mb-1)) { y_local = static_cast<real_type>(0.5)*y_local; }
      else if (ma > (mb-1)) { y_local.re = static_cast<real_type>(0.); y_local.im = static_cast<real_type>(0.); } // can probably avoid this outright
      // else the ma < mb gets "double counted", cancelling the 0.5.
    }

    const complex ulist_prev = ulist_wrapper.get(ma);
    const complex dulist_prev = dulist_wrapper.get(ma);

    // Directly accumulate deidrj into sum_tmp
    const complex du_prod = (dsfacu * ulist_prev) + (sfac * dulist_prev);
    dedr_full_sum += du_prod.re * y_local.re + du_prod.im * y_local.im;

  }

  return dedr_full_sum;
}

/* ----------------------------------------------------------------------
   calculate derivative of Ui w.r.t. atom j
   see comments above compute_duarray_cpu for comments on the
   data layout
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType, real_type, vector_length>::compute_duidrj_cpu(const int& iatom, const int& jnbor) const
{
  complex da[3], db[3];
  real_type u[3], dz0[3], dr0inv[3];

  const real_type x = rij(iatom,jnbor,0);
  const real_type y = rij(iatom,jnbor,1);
  const real_type z = rij(iatom,jnbor,2);
  const real_type rsq = x * x + y * y + z * z;
  const real_type r = sqrt(rsq);
  const real_type rscale0 = rfac0 * static_cast<real_type>(MY_PI) / (rcutij(iatom,jnbor) - rmin0);
  const real_type theta0 = (r - rmin0) * rscale0;
  const real_type sn = sin(theta0);
  const real_type cs = cos(theta0);
  const real_type z0 = r * cs / sn;
  const real_type dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;

  // begin what was compute_duarray_cpu

  real_type rinv = 1.0 / r;
  u[0] = x * rinv;
  u[1] = y * rinv;
  u[2] = z * rinv;

  real_type r0inv = 1.0 / sqrt(r * r + z0 * z0);
  complex a = { z0 * r0inv, -z * r0inv };
  complex b = { y * r0inv, -x * r0inv };

  real_type dr0invdr = -r0inv * r0inv * r0inv * (r + z0 * dz0dr);

  dr0inv[0] = dr0invdr * u[0];
  dr0inv[1] = dr0invdr * u[1];
  dr0inv[2] = dr0invdr * u[2];

  dz0[0] = dz0dr * u[0];
  dz0[1] = dz0dr * u[1];
  dz0[2] = dz0dr * u[2];

  for (int k = 0; k < 3; k++) {
    da[k].re = dz0[k] * r0inv + z0 * dr0inv[k];
    da[k].im = -z * dr0inv[k];
  }

  da[2].im += -r0inv;

  for (int k = 0; k < 3; k++) {
    db[k].re = y * dr0inv[k];
    db[k].im = -x * dr0inv[k];
  }

  db[0].im += -r0inv;
  db[1].re += r0inv;

  for (int k = 0; k < 3; k++)
    dulist_cpu(iatom, jnbor, 0, k) = complex::zero();

  for (int j = 1; j <= twojmax; j++) {
    int jju = idxu_cache_block[j];
    int jjup = idxu_cache_block[j-1];

    for (int mb = 0; 2*mb <= j; mb++) {
      int jju_index = jju + mb * (j + 1);
      int jjup_index = jjup + mb * j;

      complex duidrj[3] = { complex::zero(), complex::zero(), complex::zero() };

      for (int ma = 0; ma < j; ma++) {
        complex ui_prev = ulist_cpu(iatom, jnbor, jjup_index);
        complex duidrj_prev[3] = { dulist_cpu(iatom, jnbor, jjup_index, 0),
                                   dulist_cpu(iatom, jnbor, jjup_index, 1),
                                   dulist_cpu(iatom, jnbor, jjup_index, 2) };

        real_type rootpq = rootpqarray(j - ma,j - mb);
        for (int k = 0; k < 3; k++) {
          duidrj[k].re += rootpq * (da[k].re * ui_prev.re + da[k].im * ui_prev.im +
                                    a.re * duidrj_prev[k].re + a.im * duidrj_prev[k].im);

          duidrj[k].im += rootpq * (da[k].re * ui_prev.im - da[k].im * ui_prev.re +
                                    a.re * duidrj_prev[k].im - a.im * duidrj_prev[k].re);

          dulist_cpu(iatom, jnbor, jju_index, k) = duidrj[k];
        }

        rootpq = rootpqarray(ma + 1,j - mb);
        for (int k = 0; k < 3; k++) {
          duidrj[k].re = -rootpq * (db[k].re * ui_prev.re + db[k].im * ui_prev.im +
                                    b.re * duidrj_prev[k].re + b.im * duidrj_prev[k].im);

          duidrj[k].im = -rootpq * (db[k].re * ui_prev.im - db[k].im * ui_prev.re +
                                    b.re * duidrj_prev[k].im - b.im * duidrj_prev[k].re);
        }

        jju_index++;
        jjup_index++;
      }

      for (int k = 0; k < 3; k++) {
        dulist_cpu(iatom, jnbor, jju_index, k) = duidrj[k];
      }
    }

    // Only need to add one symmetrized row for convenience
    // Symmetry gets "unfolded" during the dedr accumulation

    // copy left side to right side with inversion symmetry VMK 4.4(2)
    // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])

    if (j % 2 == 1) {
      int mb = j / 2;
      // begin filling in the extra row
      int jju_index = jju + (mb + 1) * (j + 1);
      int jjup_index = jjup + mb * j;

      complex duidrj[3] = { complex::zero(), complex::zero(), complex::zero() };

      for (int ma = 0; ma < j; ma++) {
        complex ui_prev = ulist_cpu(iatom, jnbor, jjup_index);
        complex duidrj_prev[3] = { dulist_cpu(iatom, jnbor, jjup_index, 0),
                                   dulist_cpu(iatom, jnbor, jjup_index, 1),
                                   dulist_cpu(iatom, jnbor, jjup_index, 2) };

        real_type rootpq = rootpqarray(j - ma, mb + 1);
        for (int k = 0; k < 3; k++) {
          duidrj[k].re += rootpq * (db[k].re * ui_prev.re - db[k].im * ui_prev.im +
                                    b.re * duidrj_prev[k].re - b.im * duidrj_prev[k].im);
          duidrj[k].im += rootpq * (db[k].re * ui_prev.im + db[k].im * ui_prev.re +
                                    b.re * duidrj_prev[k].im + b.im * duidrj_prev[k].re);

          dulist_cpu(iatom, jnbor, jju_index, k) = duidrj[k];
        }

        rootpq = rootpqarray(ma + 1, mb + 1);
        for (int k = 0; k < 3; k++) {
          duidrj[k].re = rootpq * (da[k].re * ui_prev.re - da[k].im * ui_prev.im +
                                   a.re * duidrj_prev[k].re - a.im * duidrj_prev[k].im);
          duidrj[k].im = rootpq * (da[k].re * ui_prev.im + da[k].im * ui_prev.re +
                                   a.re * duidrj_prev[k].im + a.im * duidrj_prev[k].re);
        }

        jju_index++;
        jjup_index++;
      }

      for (int k = 0; k < 3; k++) {
        dulist_cpu(iatom, jnbor, jju_index, k) = duidrj[k];
      }
    }
  }

  real_type sfac, dsfac;
  compute_s_dsfac(r, rcutij(iatom,jnbor), sinnerij(iatom,jnbor), dinnerij(iatom,jnbor), sfac, dsfac);

  sfac *= wj(iatom,jnbor);
  dsfac *= wj(iatom,jnbor);

  // Even though we fill out a full "cached" data layout above,
  // we only need the "half" data for the accumulation into dedr.
  // Thus we skip updating any unnecessary data.
  for (int j = 0; j <= twojmax; j++) {
    int jju = idxu_cache_block[j];
    for (int mb = 0; 2*mb <= j; mb++)
      for (int ma = 0; ma <= j; ma++) {
        for (int k = 0; k < 3; k++) {
          dulist_cpu(iatom, jnbor, jju, k).re = dsfac * ulist_cpu(iatom, jnbor, jju).re * u[k] +
                                    sfac * dulist_cpu(iatom, jnbor, jju, k).re;
          dulist_cpu(iatom, jnbor, jju, k).im = dsfac * ulist_cpu(iatom, jnbor, jju).im * u[k] +
                                    sfac * dulist_cpu(iatom, jnbor, jju, k).im;
        }
        jju++;
      }
  }
}

/* ----------------------------------------------------------------------
   compute dEidRj, CPU path only.
   dulist takes advantage of a `cached` data layout, similar to the
   shared memory layout for the GPU routines, which is efficient for
   compressing the calculation in compute_duarray_cpu. That said,
   dulist only uses the "half" data layout part of that structure.
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType, real_type, vector_length>::compute_deidrj_cpu(const int& iatom, const int& jnbor) const
{
  real_type force_sum[3] = { 0, 0, 0 };
  const int jelem = element(iatom, jnbor);

  for (int j = 0; j <= twojmax; j++) {
    int jju_half = idxu_half_block[j];
    int jju_cache = idxu_cache_block[j];

    for (int mb = 0; 2 * mb < j; mb++) {
      for (int ma = 0; ma <= j; ma++) {
        complex y = { ylist_re(iatom, jelem, jju_half), ylist_im(iatom, jelem, jju_half) };
        for (int k = 0; k < 3; k++)
          force_sum[k] += dulist_cpu(iatom, jnbor, jju_cache, k).re * y.re +
                          dulist_cpu(iatom, jnbor, jju_cache, k).im * y.im;
        jju_half++; jju_cache++;
      }
    } //end loop over ma mb

    // For j even, handle middle column

    if (j % 2 == 0) {
      //int mb = j / 2;
      for (int ma = 0; ma <= j; ma++) {
        complex y = { ylist_re(iatom, jelem, jju_half), ylist_im(iatom, jelem, jju_half) };
        for (int k = 0; k < 3; k++)
          force_sum[k] += static_cast<real_type>(0.5) * (dulist_cpu(iatom, jnbor, jju_cache, k).re * y.re +
                          dulist_cpu(iatom, jnbor, jju_cache, k).im * y.im);
        jju_half++; jju_cache++;
      }
    } // end if jeven

  }

  for (int k = 0; k < 3; k++)
    dedr(iatom, jnbor, k) = 2 * force_sum[k];
}

/* ----------------------------------------------------------------------
   factorial n, wrapper for precomputed table
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
inline
double SNAKokkos<DeviceType, real_type, vector_length>::factorial(int n)
{
  //if (n < 0 || n > nmaxfactorial) {
  //  char str[128];
  //  sprintf(str, "Invalid argument to factorial %d", n);
  //  error->all(FLERR, str);
  //}

  return nfac_table[n];
}

/* ----------------------------------------------------------------------
   factorial n table, size SNA::nmaxfactorial+1
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
const double SNAKokkos<DeviceType, real_type, vector_length>::nfac_table[] = {
  1,
  1,
  2,
  6,
  24,
  120,
  720,
  5040,
  40320,
  362880,
  3628800,
  39916800,
  479001600,
  6227020800,
  87178291200,
  1307674368000,
  20922789888000,
  355687428096000,
  6.402373705728e+15,
  1.21645100408832e+17,
  2.43290200817664e+18,
  5.10909421717094e+19,
  1.12400072777761e+21,
  2.5852016738885e+22,
  6.20448401733239e+23,
  1.5511210043331e+25,
  4.03291461126606e+26,
  1.08888694504184e+28,
  3.04888344611714e+29,
  8.8417619937397e+30,
  2.65252859812191e+32,
  8.22283865417792e+33,
  2.63130836933694e+35,
  8.68331761881189e+36,
  2.95232799039604e+38,
  1.03331479663861e+40,
  3.71993326789901e+41,
  1.37637530912263e+43,
  5.23022617466601e+44,
  2.03978820811974e+46,
  8.15915283247898e+47,
  3.34525266131638e+49,
  1.40500611775288e+51,
  6.04152630633738e+52,
  2.65827157478845e+54,
  1.1962222086548e+56,
  5.50262215981209e+57,
  2.58623241511168e+59,
  1.24139155925361e+61,
  6.08281864034268e+62,
  3.04140932017134e+64,
  1.55111875328738e+66,
  8.06581751709439e+67,
  4.27488328406003e+69,
  2.30843697339241e+71,
  1.26964033536583e+73,
  7.10998587804863e+74,
  4.05269195048772e+76,
  2.35056133128288e+78,
  1.3868311854569e+80,
  8.32098711274139e+81,
  5.07580213877225e+83,
  3.14699732603879e+85,
  1.98260831540444e+87,
  1.26886932185884e+89,
  8.24765059208247e+90,
  5.44344939077443e+92,
  3.64711109181887e+94,
  2.48003554243683e+96,
  1.71122452428141e+98,
  1.19785716699699e+100,
  8.50478588567862e+101,
  6.12344583768861e+103,
  4.47011546151268e+105,
  3.30788544151939e+107,
  2.48091408113954e+109,
  1.88549470166605e+111,
  1.45183092028286e+113,
  1.13242811782063e+115,
  8.94618213078297e+116,
  7.15694570462638e+118,
  5.79712602074737e+120,
  4.75364333701284e+122,
  3.94552396972066e+124,
  3.31424013456535e+126,
  2.81710411438055e+128,
  2.42270953836727e+130,
  2.10775729837953e+132,
  1.85482642257398e+134,
  1.65079551609085e+136,
  1.48571596448176e+138,
  1.3520015276784e+140,
  1.24384140546413e+142,
  1.15677250708164e+144,
  1.08736615665674e+146,
  1.03299784882391e+148,
  9.91677934870949e+149,
  9.61927596824821e+151,
  9.42689044888324e+153,
  9.33262154439441e+155,
  9.33262154439441e+157,
  9.42594775983835e+159,
  9.61446671503512e+161,
  9.90290071648618e+163,
  1.02990167451456e+166,
  1.08139675824029e+168,
  1.14628056373471e+170,
  1.22652020319614e+172,
  1.32464181945183e+174,
  1.44385958320249e+176,
  1.58824554152274e+178,
  1.76295255109024e+180,
  1.97450685722107e+182,
  2.23119274865981e+184,
  2.54355973347219e+186,
  2.92509369349301e+188,
  3.3931086844519e+190,
  3.96993716080872e+192,
  4.68452584975429e+194,
  5.5745857612076e+196,
  6.68950291344912e+198,
  8.09429852527344e+200,
  9.8750442008336e+202,
  1.21463043670253e+205,
  1.50614174151114e+207,
  1.88267717688893e+209,
  2.37217324288005e+211,
  3.01266001845766e+213,
  3.8562048236258e+215,
  4.97450422247729e+217,
  6.46685548922047e+219,
  8.47158069087882e+221,
  1.118248651196e+224,
  1.48727070609069e+226,
  1.99294274616152e+228,
  2.69047270731805e+230,
  3.65904288195255e+232,
  5.01288874827499e+234,
  6.91778647261949e+236,
  9.61572319694109e+238,
  1.34620124757175e+241,
  1.89814375907617e+243,
  2.69536413788816e+245,
  3.85437071718007e+247,
  5.5502938327393e+249,
  8.04792605747199e+251,
  1.17499720439091e+254,
  1.72724589045464e+256,
  2.55632391787286e+258,
  3.80892263763057e+260,
  5.71338395644585e+262,
  8.62720977423323e+264,
  1.31133588568345e+267,
  2.00634390509568e+269,
  3.08976961384735e+271,
  4.78914290146339e+273,
  7.47106292628289e+275,
  1.17295687942641e+278,
  1.85327186949373e+280,
  2.94670227249504e+282,
  4.71472363599206e+284,
  7.59070505394721e+286,
  1.22969421873945e+289,
  2.0044015765453e+291,
  3.28721858553429e+293,
  5.42391066613159e+295,
  9.00369170577843e+297,
  1.503616514865e+300, // nmaxfactorial = 167
};

/* ----------------------------------------------------------------------
   the function delta given by VMK Eq. 8.2(1)
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
inline
double SNAKokkos<DeviceType, real_type, vector_length>::deltacg(int j1, int j2, int j)
{
  double sfaccg = factorial((j1 + j2 + j) / 2 + 1);
  return sqrt(factorial((j1 + j2 - j) / 2) *
              factorial((j1 - j2 + j) / 2) *
              factorial((-j1 + j2 + j) / 2) / sfaccg);
}

/* ----------------------------------------------------------------------
   assign Clebsch-Gordan coefficients using
   the quasi-binomial formula VMK 8.2.1(3)
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
inline
void SNAKokkos<DeviceType, real_type, vector_length>::init_clebsch_gordan()
{
  auto h_cglist = Kokkos::create_mirror_view(cglist);

  double sum,dcg,sfaccg;
  int m, aa2, bb2, cc2;
  int ifac;

  int idxcg_count = 0;
  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        for (int m1 = 0; m1 <= j1; m1++) {
          aa2 = 2 * m1 - j1;

          for (int m2 = 0; m2 <= j2; m2++) {

            // -c <= cc <= c

            bb2 = 2 * m2 - j2;
            m = (aa2 + bb2 + j) / 2;

            if (m < 0 || m > j) {
              h_cglist[idxcg_count] = 0.0;
              idxcg_count++;
              continue;
            }

            sum = 0.0;

            for (int z = MAX(0, MAX(-(j - j2 + aa2)
                                    / 2, -(j - j1 - bb2) / 2));
                 z <= MIN((j1 + j2 - j) / 2,
                          MIN((j1 - aa2) / 2, (j2 + bb2) / 2));
                 z++) {
              ifac = z % 2 ? -1 : 1;
              sum += ifac /
                (factorial(z) *
                 factorial((j1 + j2 - j) / 2 - z) *
                 factorial((j1 - aa2) / 2 - z) *
                 factorial((j2 + bb2) / 2 - z) *
                 factorial((j - j2 + aa2) / 2 + z) *
                 factorial((j - j1 - bb2) / 2 + z));
            }

            cc2 = 2 * m - j;
            dcg = deltacg(j1, j2, j);
            sfaccg = sqrt(factorial((j1 + aa2) / 2) *
                          factorial((j1 - aa2) / 2) *
                          factorial((j2 + bb2) / 2) *
                          factorial((j2 - bb2) / 2) *
                          factorial((j  + cc2) / 2) *
                          factorial((j  - cc2) / 2) *
                          (j + 1));

            h_cglist[idxcg_count] = sum * dcg * sfaccg;
            idxcg_count++;
          }
        }
      }
  Kokkos::deep_copy(cglist,h_cglist);
}

/* ----------------------------------------------------------------------
   pre-compute table of sqrt[p/m2], p, q = 1,twojmax
   the p = 0, q = 0 entries are allocated and skipped for convenience.
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
inline
void SNAKokkos<DeviceType, real_type, vector_length>::init_rootpqarray()
{
  auto h_rootpqarray = Kokkos::create_mirror_view(rootpqarray);
  for (int p = 1; p <= twojmax; p++)
    for (int q = 1; q <= twojmax; q++)
      h_rootpqarray(p,q) = static_cast<real_type>(sqrt(static_cast<double>(p)/q));
  Kokkos::deep_copy(rootpqarray,h_rootpqarray);
}


/* ---------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
inline
int SNAKokkos<DeviceType, real_type, vector_length>::compute_ncoeff()
{
  int ncount;

  ncount = 0;

  for (int j1 = 0; j1 <= twojmax; j1++)
    for (int j2 = 0; j2 <= j1; j2++)
      for (int j = j1 - j2;
           j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) ncount++;

  ndoubles = nelements*nelements;
  ntriples = nelements*nelements*nelements;
  if (chem_flag) ncount *= ntriples;

  return ncount;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
real_type SNAKokkos<DeviceType, real_type, vector_length>::compute_sfac(real_type r, real_type rcut, real_type sinner, real_type  dinner) const
{
  real_type sfac_outer;
  constexpr real_type one = static_cast<real_type>(1.0);
  constexpr real_type zero = static_cast<real_type>(0.0);
  constexpr real_type onehalf = static_cast<real_type>(0.5);
  if (switch_flag == 0) sfac_outer = one;
  if (switch_flag == 1) {
    if (r <= rmin0) sfac_outer = one;
    else if (r > rcut) return zero;
    else {
      real_type rcutfac = static_cast<real_type>(MY_PI) / (rcut - rmin0);
      sfac_outer = onehalf * (cos((r - rmin0) * rcutfac) + one);
    }
  }

  if (switch_inner_flag == 0) return sfac_outer;
  if (switch_inner_flag == 1) {
    if (r >= sinner + dinner)
        return sfac_outer;
    else if (r > sinner - dinner) {
      real_type rcutfac = static_cast<real_type>(MY_PI2) / dinner;
      return sfac_outer *
        onehalf * (one - cos(static_cast<real_type>(MY_PI2) + (r - sinner) * rcutfac));
    } else return zero;
  }
  return zero; // dummy return
}

/* ---------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
real_type SNAKokkos<DeviceType, real_type, vector_length>::compute_dsfac(real_type r, real_type rcut, real_type sinner, real_type dinner) const
{
  real_type sfac_outer, dsfac_outer, sfac_inner, dsfac_inner;
  constexpr real_type one = static_cast<real_type>(1.0);
  constexpr real_type zero = static_cast<real_type>(0.0);
  constexpr real_type onehalf = static_cast<real_type>(0.5);
  if (switch_flag == 0) dsfac_outer = zero;
  if (switch_flag == 1) {
    if (r <= rmin0) dsfac_outer = zero;
    else if (r > rcut) return zero;
    else {
      real_type rcutfac = static_cast<real_type>(MY_PI) / (rcut - rmin0);
      dsfac_outer = -onehalf * sin((r - rmin0) * rcutfac) * rcutfac;
    }
  }

  if (switch_inner_flag == 0) return dsfac_outer;
  if (switch_inner_flag == 1) {
    if (r >= sinner + dinner)
      return dsfac_outer;
    else if (r > sinner - dinner) {

      // calculate sfac_outer

      if (switch_flag == 0) sfac_outer = one;
      if (switch_flag == 1) {
        if (r <= rmin0) sfac_outer = one;
        else if (r > rcut) sfac_outer = zero;
        else {
          real_type rcutfac = static_cast<real_type>(MY_PI) / (rcut - rmin0);
          sfac_outer = onehalf * (cos((r - rmin0) * rcutfac) + one);
        }
      }

      // calculate sfac_inner

      real_type rcutfac = static_cast<real_type>(MY_PI2) / dinner;
      sfac_inner = onehalf * (one - cos(static_cast<real_type>(MY_PI2) + (r - sinner) * rcutfac));
      dsfac_inner = onehalf * rcutfac * sin(static_cast<real_type>(MY_PI2) + (r - sinner) * rcutfac);
      return dsfac_outer * sfac_inner + sfac_outer * dsfac_inner;

    } else return zero;
  }
  return zero; // dummy return
}

template<class DeviceType, typename real_type, int vector_length>
KOKKOS_INLINE_FUNCTION
void SNAKokkos<DeviceType, real_type, vector_length>::compute_s_dsfac(const real_type r, const real_type rcut, const real_type sinner, const real_type dinner, real_type& sfac, real_type& dsfac) const {

  real_type sfac_outer, dsfac_outer, sfac_inner, dsfac_inner;
  constexpr real_type one = static_cast<real_type>(1.0);
  constexpr real_type zero = static_cast<real_type>(0.0);
  constexpr real_type onehalf = static_cast<real_type>(0.5);

  if (switch_flag == 0) { sfac_outer = one; dsfac_outer = zero; }
  else if (switch_flag == 1) {
    if (r <= rmin0) { sfac_outer = one; dsfac_outer = zero; }
    else if (r > rcut) { sfac = zero; dsfac = zero; return; }
    else {
      const real_type rcutfac = static_cast<real_type>(MY_PI) / (rcut - rmin0);
      const real_type theta0 = (r - rmin0) * rcutfac;
      const real_type sn = sin(theta0);
      const real_type cs = cos(theta0);
      sfac_outer = onehalf * (cs + one);
      dsfac_outer = -onehalf * sn * rcutfac;
    }
  } else { sfac = zero; dsfac = zero; return; } // dummy return

  if (switch_inner_flag == 0) { sfac = sfac_outer; dsfac = dsfac_outer; return; }
  else if (switch_inner_flag == 1) {
    if (r >= sinner + dinner) { sfac = sfac_outer; dsfac = dsfac_outer; return; }
    else if (r > sinner - dinner) {
      real_type rcutfac = static_cast<real_type>(MY_PI2) / dinner;
      sfac_inner = onehalf * (one - cos(static_cast<real_type>(MY_PI2) + (r - sinner) * rcutfac));
      dsfac_inner = onehalf * rcutfac * sin(static_cast<real_type>(MY_PI2) + (r - sinner) * rcutfac);
      sfac = sfac_outer * sfac_inner;
      dsfac = dsfac_outer * sfac_inner + sfac_outer * dsfac_inner;
      return;
    } else { sfac = zero; dsfac = zero; return; }
  } else { sfac = zero; dsfac = zero; return; } // dummy return

}

/* ----------------------------------------------------------------------
   memory usage of arrays
------------------------------------------------------------------------- */

template<class DeviceType, typename real_type, int vector_length>
double SNAKokkos<DeviceType, real_type, vector_length>::memory_usage()
{
  double bytes = 0;

  bytes += MemKK::memory_usage(rootpqarray);
  bytes += MemKK::memory_usage(cglist);

  bytes += MemKK::memory_usage(ulisttot_re);
  bytes += MemKK::memory_usage(ulisttot_im);
  bytes += MemKK::memory_usage(ulisttot);

  bytes += MemKK::memory_usage(zlist);
  bytes += MemKK::memory_usage(blist);

  bytes += MemKK::memory_usage(ylist_re);
  bytes += MemKK::memory_usage(ylist_im);

  if constexpr (!host_flag) {
    bytes += MemKK::memory_usage(a_gpu);
    bytes += MemKK::memory_usage(b_gpu);
    bytes += MemKK::memory_usage(da_gpu);
    bytes += MemKK::memory_usage(db_gpu);
    bytes += MemKK::memory_usage(sfac_gpu);
  } else {
    bytes += MemKK::memory_usage(ulist_cpu);
    bytes += MemKK::memory_usage(dulist_cpu);
  }

  bytes += MemKK::memory_usage(dedr);

  bytes += MemKK::memory_usage(idxcg_block);
  bytes += MemKK::memory_usage(idxu_block);
  bytes += MemKK::memory_usage(idxu_half_block);
  bytes += MemKK::memory_usage(idxu_full_half);
  bytes += MemKK::memory_usage(idxu_cache_block);
  bytes += MemKK::memory_usage(idxz_block);
  bytes += MemKK::memory_usage(idxb_block);

  bytes += MemKK::memory_usage(idxz);
  bytes += MemKK::memory_usage(idxb);

  bytes += MemKK::memory_usage(bzero);

  bytes += MemKK::memory_usage(rij);
  bytes += MemKK::memory_usage(inside);
  bytes += MemKK::memory_usage(wj);
  bytes += MemKK::memory_usage(rcutij);
  bytes += MemKK::memory_usage(sinnerij);
  bytes += MemKK::memory_usage(dinnerij);

  return bytes;
}

} // namespace LAMMPS_NS
