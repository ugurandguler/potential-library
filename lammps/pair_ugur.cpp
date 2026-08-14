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
   Contributing authors: Gokay Ugur, Sule Ugur (Gazi University),
   Melek Guler, Emre Guler (Ankara Haci Bayram Veli University)
   - Akgun-Ugur pair styles.

   The file-reading and parameter-indexing structure (read_file, setup_params,
   elem3param, the one-line-per-ordered-triple format) follows LAMMPS's own
   Tersoff style, pair_tersoff.cpp, by Aidan Thompson (SNL).  That is deliberate
   and not incidental: phi3 depends on its two legs through r1 + r2 and does not
   factorise into per-leg pieces, so its parameters belong to a triple rather
   than to a pair - the same shape of problem Tersoff has, so it gets the same
   file format and the same machinery for reading it.

   The physics is in ugurpot.h and is original work, validated against
   standalone/latdyn.py to 1e-11 on energy and pressure.
------------------------------------------------------------------------- */

/* See pair_ugur.h.  The physics is in ugurpot.c; this walks neighbours. */
#include "pair_ugur.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>
#include <string>

static constexpr int NPARAMS_PER_LINE = 15;
static constexpr int DELTA = 4;

using namespace LAMMPS_NS;

PairUgur::PairUgur(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;      /* the three-body term has no pair decomposition */
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  cutmax = 0.0;
  params = nullptr;
  nparams = maxparam = 0;
  elem3param = nullptr;
  angular = 0;
}

PairUgur::~PairUgur()
{
  if (copymode) return;
  memory->destroy(params);
  memory->destroy(elem3param);
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(map);
  }
}

void PairUgur::allocate()
{
  allocated = 1;
  const int n = atom->ntypes + 1;
  memory->create(setflag, n, n, "pair:setflag");
  memory->create(cutsq, n, n, "pair:cutsq");
  /* map is what map_element2type writes the type-to-element table into, and
     it has to exist first.  Leaving it out does not fail loudly - LAMMPS
     stops after printing its banner, with no error line anywhere. */
  memory->create(map, n, "pair:map");
  for (int i = 1; i < n; i++)
    for (int j = i; j < n; j++) setflag[i][j] = 0;
}

void PairUgur::settings(int narg, char **)
{
  if (narg != 0) error->all(FLERR, "pair_style ugur takes no arguments");
}

/* Parameter file: one line of ten numbers, in the order the fit stores them.
 *
 *   m D alpha r0 gamma C alpha3 rcut2 rcut3 taper
 *
 * `taper` is not optional and a negative value means hard truncation.  It is
 * written by the exporter rather than defaulted here on purpose: a potential
 * fitted with the switch on and run with it off is a different potential, and
 * that class of mistake has already cost this project once, with the
 * three-body cutoff.
 */
void PairUgur::read_file(char *file)
{
  memory->sfree(params);
  params = nullptr;
  nparams = maxparam = 0;

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, file, "ugur", unit_convert_flag);
    char *line;
    while ((line = reader.next_line(NPARAMS_PER_LINE))) {
      try {
        ValueTokenizer v(line);
        const std::string e1 = v.next_string();
        const std::string e2 = v.next_string();
        const std::string e3 = v.next_string();
        /* keep only the triples whose three labels are all in the pair_coeff
           list; a file may legitimately carry more elements than a run uses */
        int i, j, k;
        for (i = 0; i < nelements; i++) if (e1 == elements[i]) break;
        if (i == nelements) continue;
        for (j = 0; j < nelements; j++) if (e2 == elements[j]) break;
        if (j == nelements) continue;
        for (k = 0; k < nelements; k++) if (e3 == elements[k]) break;
        if (k == nelements) continue;

        if (nparams == maxparam) {
          maxparam += DELTA;
          params = (Param *) memory->srealloc(params, maxparam * sizeof(Param),
                                              "pair:params");
          std::memset(params + nparams, 0, DELTA * sizeof(Param));
        }
        Param &pm = params[nparams];
        pm.ielement = i;
        pm.jelement = j;
        pm.kelement = k;
        pm.p.m = v.next_double();
        pm.p.D = v.next_double();
        pm.p.alpha = v.next_double();
        pm.p.r0 = v.next_double();
        pm.p.gamma = v.next_double();
        pm.p.C = v.next_double();
        pm.p.alpha3 = v.next_double();
        pm.p.rcut2 = v.next_double();
        pm.p.rcut3 = v.next_double();
        pm.p.taper = v.next_double();
        /* The UG Legendre weights are in the format so that a second,
           incompatible format is not needed later.  This style implements the
           published angle-free phi3, so it refuses a file that carries them
           rather than silently ignoring a term the user asked for. */
        pm.p.lam2 = v.has_next() ? v.next_double() : 0.0;
        pm.p.lam4 = v.has_next() ? v.next_double() : 0.0;
        if ((pm.p.lam2 != 0.0 || pm.p.lam4 != 0.0) && !angular)
          error->one(FLERR, "pair_style ugur implements the angle-free phi3; "
                            "use pair_style ugur/ang for a file with nonzero "
                            "lam2/lam4");
        nparams++;
      } catch (TokenizerException &e) {
        error->one(FLERR, "Invalid ugur potential file: {}", e.what());
      }
    }
  }
  MPI_Bcast(&nparams, 1, MPI_INT, 0, world);
  MPI_Bcast(&maxparam, 1, MPI_INT, 0, world);
  if (comm->me != 0) {
    params = (Param *) memory->srealloc(params, maxparam * sizeof(Param),
                                        "pair:params");
  }
  MPI_Bcast(params, maxparam * sizeof(Param), MPI_BYTE, 0, world);
  setup_params();
}

/* Index the triples and check the file is complete.  A missing entry is a
   silent wrong answer otherwise: the run would use whatever happened to be in
   elem3param, so it is an error here instead. */
void PairUgur::setup_params()
{
  memory->destroy(elem3param);
  memory->create(elem3param, nelements, nelements, nelements, "pair:elem3param");

  for (int i = 0; i < nelements; i++)
    for (int j = 0; j < nelements; j++)
      for (int k = 0; k < nelements; k++) {
        int n = -1;
        for (int m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement &&
              k == params[m].kelement) {
            if (n >= 0) error->all(FLERR, "Duplicate entry for {} {} {} in "
                                          "ugur potential file",
                                   elements[i], elements[j], elements[k]);
            n = m;
          }
        }
        if (n < 0) error->all(FLERR, "No entry for {} {} {} in ugur "
                                     "potential file",
                              elements[i], elements[j], elements[k]);
        elem3param[i][j][k] = n;
      }

  /* One window for the whole file.  Mixing taper values would be a different
     model per bond, not a potential, and it is exactly the kind of thing that
     runs happily and means nothing. */
  const double t0 = params[0].p.taper;
  for (int m = 1; m < nparams; m++)
    if (params[m].p.taper != t0)
      error->all(FLERR, "ugur potential file mixes taper values ({} and {}); "
                        "the switch window is part of the model, not a "
                        "per-bond parameter", t0, params[m].p.taper);

  cutmax = 0.0;
  for (int m = 0; m < nparams; m++) {
    if (params[m].p.rcut2 > cutmax) cutmax = params[m].p.rcut2;
    if (params[m].p.rcut3 > cutmax) cutmax = params[m].p.rcut3;
  }
}

void PairUgur::coeff(int narg, char **arg)
{
  if (!allocated) allocate();
  map_element2type(narg - 3, arg + 3);
  read_file(arg[2]);
}

void PairUgur::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR, "pair_style ugur requires newton pair on");
  /* full list: the three-body term needs every neighbour of the centre atom,
     not half of them */
  neighbor->add_request(this, NeighConst::REQ_FULL);
}

double PairUgur::init_one(int, int) { return cutmax; }

void PairUgur::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  const int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  int *type = atom->type;

  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    const int ei = map[type[i]];
    const double xi = x[i][0], yi = x[i][1], zi = x[i][2];
    int *jlist = firstneigh[i];
    const int jnum = numneigh[i];

    /* ---- pair term ------------------------------------------------------
       A full list visits every bond twice, so each visit carries half.     */
    for (int jj = 0; jj < jnum; jj++) {
      const int j = jlist[jj] & NEIGHMASK;
      /* the pair term for an unlike bond lives on the (i,j,j) entry, which is
         Tersoff's rule and is symmetric in i,j only if the file is - the
         generator writes it so and setup_params does not assume it */
      const ugur_param *pp = pair_par(ei, map[type[j]]);
      const double dx = xi - x[j][0], dy = yi - x[j][1], dz = zi - x[j][2];
      const double rsq = dx * dx + dy * dy + dz * dz;
      if (rsq >= pp->rcut2 * pp->rcut2) continue;
      const double r = std::sqrt(rsq);
      double dphi;
      const double e = ugur_phi2(pp, r, &dphi);
      const double fpair = -dphi / r;
      f[i][0] += 0.5 * fpair * dx;
      f[i][1] += 0.5 * fpair * dy;
      f[i][2] += 0.5 * fpair * dz;
      f[j][0] -= 0.5 * fpair * dx;
      f[j][1] -= 0.5 * fpair * dy;
      f[j][2] -= 0.5 * fpair * dz;
      if (evflag) ev_tally(i, j, atom->nlocal, 1, 0.5 * e, 0.0,
                           0.5 * fpair, dx, dy, dz);
    }

    /* ---- three-body term -------------------------------------------------
       Unordered leg pairs about the centre i, each pair once, matching
       latdyn.triplets.  phi3 depends on the legs only through r1 + r2, so
       there is no angular derivative: each leg is pushed along its own unit
       vector and the centre takes the reaction.                            */
    {
      for (int jj = 0; jj < jnum; jj++) {
        const int j = jlist[jj] & NEIGHMASK;
        const int ej = map[type[j]];
        /* the switch acts per leg, so its cutoff belongs to the bond i-j and
           not to the triplet */
        const double rc3a = pair_par(ei, ej)->rcut3;
        const double ax = x[j][0] - xi, ay = x[j][1] - yi, az = x[j][2] - zi;
        const double rasq = ax * ax + ay * ay + az * az;
        if (rasq >= rc3a * rc3a) continue;
        const double ra = std::sqrt(rasq);
        for (int kk = jj + 1; kk < jnum; kk++) {
          const int k = jlist[kk] & NEIGHMASK;
          const int ek = map[type[k]];
          /* C and alpha3 belong to the triple; the radial shape - D, m, r0,
             gamma - is taken from the CENTRE's own entry.  Both halves of that
             are forced.  A triple must be symmetric under swapping its legs,
             and the (i,j,k) and (i,k,j) lines carry different pair columns by
             construction, so reading the shape off the triple line would give
             one triplet two different energies depending on which neighbour
             the loop reached first.  The centre is the only symmetric source
             available without adding columns, and it is the atom whose bonds
             are being bent.  For one element it is the same entry, so nothing
             about the single-species validation changes. */
          ugur_param tri = *pair_par(ei, ei);
          {
            const ugur_param *t3 = &params[elem3param[ei][ej][ek]].p;
            tri.C = t3->C;
            tri.alpha3 = t3->alpha3;
          }
          const ugur_param *tp = &tri;
          if (tp->C == 0.0) continue;
          const double rc3b = pair_par(ei, ek)->rcut3;
          const double bx = x[k][0] - xi, by = x[k][1] - yi,
                       bz = x[k][2] - zi;
          const double rbsq = bx * bx + by * by + bz * bz;
          if (rbsq >= rc3b * rc3b) continue;
          const double rb = std::sqrt(rbsq);
          double d1, d2;
          /* dE/dR_j = d1 * n_a, so the force on j is -d1 * n_a.  Written
             into fj and fk in the layout ev_tally3 expects - the same
             convention Stillinger-Weber uses: delr1 runs i to j, delr2 runs
             i to k, the two arrays hold the leg forces and the centre takes
             the reaction. */
          double dcos = 0.0, e;
          if (angular) {
            const double cth = (ax * bx + ay * by + az * bz) / (ra * rb);
            e = ugur_phi3_ang(tp, ra, rb, cth, rc3a, rc3b, &d1, &d2, &dcos);
            /* The angular part is where the force stops being along the legs.
                 c = (a.b)/(ra rb)
                 dc/da = b/(ra rb) - c a/ra^2,   and the mirror for dc/db
               Both terms are needed: dropping the second leaves a force with
               a radial component that does not belong to it, and the run
               still looks plausible while conserving nothing. */
            const double ia = 1.0 / ra, ib = 1.0 / rb, iab = ia * ib;
            const double fa = d1 * ia, fb = d2 * ib;
            double fj[3], fk[3], delr1[3], delr2[3];
            delr1[0] = ax;  delr1[1] = ay;  delr1[2] = az;
            delr2[0] = bx;  delr2[1] = by;  delr2[2] = bz;
            fj[0] = -(fa * ax + dcos * (bx * iab - cth * ax * ia * ia));
            fj[1] = -(fa * ay + dcos * (by * iab - cth * ay * ia * ia));
            fj[2] = -(fa * az + dcos * (bz * iab - cth * az * ia * ia));
            fk[0] = -(fb * bx + dcos * (ax * iab - cth * bx * ib * ib));
            fk[1] = -(fb * by + dcos * (ay * iab - cth * by * ib * ib));
            fk[2] = -(fb * bz + dcos * (az * iab - cth * bz * ib * ib));
            f[j][0] += fj[0];  f[j][1] += fj[1];  f[j][2] += fj[2];
            f[k][0] += fk[0];  f[k][1] += fk[1];  f[k][2] += fk[2];
            f[i][0] -= fj[0] + fk[0];
            f[i][1] -= fj[1] + fk[1];
            f[i][2] -= fj[2] + fk[2];
            if (evflag) ev_tally3(i, j, k, e, 0.0, fj, fk, delr1, delr2);
            continue;
          }
          e = ugur_phi3_ab(tp, ra, rb, rc3a, rc3b, &d1, &d2);
          const double fa = d1 / ra, fb = d2 / rb;
          double fj[3], fk[3], delr1[3], delr2[3];
          delr1[0] = ax;  delr1[1] = ay;  delr1[2] = az;
          delr2[0] = bx;  delr2[1] = by;  delr2[2] = bz;
          fj[0] = -fa * ax;  fj[1] = -fa * ay;  fj[2] = -fa * az;
          fk[0] = -fb * bx;  fk[1] = -fb * by;  fk[2] = -fb * bz;
          f[j][0] += fj[0];  f[j][1] += fj[1];  f[j][2] += fj[2];
          f[k][0] += fk[0];  f[k][1] += fk[1];  f[k][2] += fk[2];
          f[i][0] -= fj[0] + fk[0];
          f[i][1] -= fj[1] + fk[1];
          f[i][2] -= fj[2] + fk[2];
          if (evflag) ev_tally3(i, j, k, e, 0.0, fj, fk, delr1, delr2);
        }
      }
    }
  }
  if (vflag_fdotr) virial_fdotr_compute();
}
