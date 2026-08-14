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

/* LAMMPS pair styles for the Akgun-Ugur family.
 *
 * Thin by design: every number comes from ugurpot.h, which is validated
 * against standalone/latdyn.py.  Nothing here computes physics; it walks
 * neighbours and hands pairs and triplets to the kernel.
 *
 *     pair_style  ugur          the published angle-free phi3 (AU / MAU)
 *     pair_style  ugur/ang      the same plus h(cos) = 1 + lam2 P2 + lam4 P4
 *     pair_coeff  * * CuNi.ugur.alloy Cu Ni
 *
 * Both are multi-species and read the same file, in LAMMPS's Tersoff format:
 * one line per ordered triple (centre, leg, leg), two-body parameters taken
 * only from the lines whose two legs agree.  See ALLOYS.md for why the problem
 * has that shape.
 *
 * They are two styles rather than one flag because they are two models.  With
 * lam2 = lam4 = 0 they compute the same thing, and `ugur` refuses a file that
 * carries nonzero weights rather than ignoring a term the user asked for -
 * silently dropping physics that is written in the file is the failure mode
 * worth engineering against.
 */
#ifdef PAIR_CLASS
// clang-format off
PairStyle(ugur, PairUgur);
PairStyle(ugur/ang, PairUgurAng);
// clang-format on
#else

#ifndef LMP_PAIR_UGUR_H
#define LMP_PAIR_UGUR_H

#include "pair.h"

#include "ugurpot.h"

namespace LAMMPS_NS {

class PairUgur : public Pair {
 public:
  PairUgur(class LAMMPS *);
  ~PairUgur() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

 protected:
  /* Storage follows Tersoff's, because the problem has Tersoff's shape:
     phi3 depends on the legs only through r1 + r2 and so does not factorise
     into per-leg pieces, which means C and alpha3 belong to an ordered triple
     (centre; leg, leg) and cannot be built from pair quantities.  params[]
     holds one entry per triple as read; elem3param maps (i,j,k) element
     indices into it.  Two-body parameters are taken only from the entries
     with j == k, exactly as Tersoff does, so a file written for one style
     reads like a file written for the other. */
  struct Param {
    ugur_param p;
    int ielement, jelement, kelement;
  };
  Param *params;
  int nparams, maxparam;
  int angular;                 /* 0: published phi3; 1: with h(cos theta) */
  int ***elem3param;
  double cutmax;

  void allocate();
  void read_file(char *);
  void setup_params();
  /* the (i,j,j) entry, which is where the pair term lives */
  inline const ugur_param *pair_par(int i, int j) const
  {
    return &params[elem3param[i][j][j]].p;
  }
};

/* The UG form.  Everything is inherited; only the flag differs, so the two
   styles cannot drift apart in the parts they share. */
class PairUgurAng : public PairUgur {
 public:
  PairUgurAng(class LAMMPS *lmp) : PairUgur(lmp) { angular = 1; }
};

}    // namespace LAMMPS_NS
#endif
#endif
