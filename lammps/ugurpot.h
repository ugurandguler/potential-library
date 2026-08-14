/* Copyright (C) 2026  Gokay Ugur, Sule Ugur, Melek Guler and Emre Guler
 * Gazi University and Ankara Haci Bayram Veli University
 *
 * This file is part of the Ugur interatomic potential library and is free
 * software: you may redistribute it and/or modify it under the terms of the
 * GNU General Public License, version 2, as published by the Free Software
 * Foundation.  See COPYING at the root of the repository.
 *
 * It has no LAMMPS dependency - header-only C99 over libm - and is the one
 * place the physics lives.  The GPL applies because the repository is
 * distributed as one work with a LAMMPS pair style in it.
 */

/* ugurpot - the Akgun-Ugur / modified Akgun-Ugur potential, energy and forces.
 *
 * Header-only C99, no dependencies beyond libm.  The physics lives here and
 * nowhere else; the LAMMPS pair style, a Python binding and anything else are
 * thin wrappers over this file, so a paper can point at one implementation
 * rather than at a patch against a particular release.
 *
 * Header-only on purpose.  A separate .c would be a second translation unit to
 * wire into LAMMPS's build, and a contributed style is easier to accept when
 * it adds sources rather than build rules.  It also lets the compiler inline
 * the kernel into the neighbour loop, which is where all the time goes.
 *
 * The reference implementation is standalone/latdyn.py.  Any change here must
 * still reproduce its energies and forces - validate_kernel.py does exactly
 * that comparison and is the only thing that says this file is correct.
 *
 *     phi2(r)      = D (r0/r)^g [ e^{m a (r0-r)} - m e^{a (r0-r)} ] / (m-1)
 *     phi3(r1,r2)  = C D (r0/x)^g [ e^{m a3 (r0-x)} - m e^{a3 (r0-x)} ] / (m-1)
 *                    S(r1) S(r2),    x = r1 + r2
 *
 * phi3 depends on the two legs only through their sum, so it has no angular
 * derivatives - three lines where Stillinger-Weber and Tersoff need pages.
 * The same property is what puts a floor under the C44/C' ratio the form can
 * reach, so the simplicity here and the limitation there are one fact.
 *
 * S is the cutoff switch.  Without it phi2 does not vanish at its cutoff - for
 * ruthenium it is -0.346 eV there - and a neighbour crossing the sphere changes
 * the energy in one step, which is invisible at a fixed geometry and fatal in
 * dynamics.  taper <= 0 disables it and reproduces hard truncation exactly.
 */
#ifndef UGURPOT_H
#define UGURPOT_H

#include <math.h>

typedef struct {
    double m, D, alpha, r0, gamma;   /* pair term                       */
    double C, alpha3;                /* three-body weight and range     */
    double rcut2, rcut3;             /* cutoffs (angstrom)              */
    double taper;                    /* switch start, as a fraction of  */
                                     /* each cutoff; <= 0 turns it off  */
    double lam2, lam4;               /* UG Legendre weights; both zero  */
                                     /* is the published angle-free form*/
} ugur_param;

/* The shape function and its radial derivative, shared by both terms.
 *
 *     f(r) = (r0/r)^g A(r),   A = [e^{m u} - m e^{u}]/(m-1),  u = al (r0 - r)
 *
 * written exactly as latdyn.Potential._f so the two can be compared term by
 * term when they disagree. */
static inline double ugur_shape(const ugur_param *p, double r, double al,
                                double *df)
{
    const double m = p->m, r0 = p->r0, g = p->gamma;
    const double u = al * (r0 - r);
    const double e1 = exp(u), em = exp(m * u);
    const double A = (em - m * e1) / (m - 1.0);
    const double P = pow(r0 / r, g);
    if (df) {
        const double dA = m * al * (e1 - em) / (m - 1.0);
        const double dP = -g * P / r;
        *df = P * dA + dP * A;
    }
    return P * A;
}

/* The switch itself, exposed for testing against latdyn.Potential.switch. */
static inline double ugur_switch(const ugur_param *p, double r, double rc,
                                 double *ds)
{
    if (p->taper <= 0.0 || rc <= 0.0) { if (ds) *ds = 0.0; return 1.0; }
    const double r_on = p->taper * rc;
    if (r <= r_on) { if (ds) *ds = 0.0; return 1.0; }
    if (r >= rc)   { if (ds) *ds = 0.0; return 0.0; }
    const double w = rc - r_on;
    const double t = (r - r_on) / w;
    /* quintic: S(0)=1, S(1)=0, S' and S'' both zero at each end, so phi'' -
     * and therefore the dynamical matrix - stays continuous across the join */
    if (ds) *ds = -30.0 * t * t * (1.0 - t) * (1.0 - t) / w;
    return 1.0 - t * t * t * (10.0 - 15.0 * t + 6.0 * t * t);
}

/* Pair term.  Returns phi2(r); *dphi receives dphi2/dr when non-NULL. */
static inline double ugur_phi2(const ugur_param *p, double r, double *dphi)
{
    if (p->rcut2 > 0.0 && r >= p->rcut2) { if (dphi) *dphi = 0.0; return 0.0; }
    double df = 0.0, dS = 0.0;
    const double f = ugur_shape(p, r, p->alpha, dphi ? &df : 0);
    const double S = ugur_switch(p, r, p->rcut2, dphi ? &dS : 0);
    if (dphi) *dphi = p->D * (df * S + f * dS);
    return p->D * f * S;
}

/* Three-body term for one triplet, switches included.  Returns the energy;
 * *d1 and *d2 receive dE/dr1 and dE/dr2 when non-NULL.  These differ once the
 * switch is on, because the term stops being a function of r1 + r2 alone. */
/* Two cutoffs, because the switch acts per leg and in an alloy the two legs
 * are different bonds: S_AB(r1) uses rcut3 of the A-B pair and S_AC(r2) that
 * of A-C, while everything else in the triplet - C, alpha3, m, r0, gamma -
 * comes from the (A;B,C) entry.  With rc1 == rc2 this is the single-species
 * function, which is what ugur_phi3 below calls it as, so the whole existing
 * validation chain still measures the same code. */
static inline double ugur_phi3_ab(const ugur_param *p, double r1, double r2,
                                  double rc1, double rc2,
                                  double *d1, double *d2)
{
    if (p->C == 0.0) { if (d1) *d1 = 0.0; if (d2) *d2 = 0.0; return 0.0; }
    if ((rc1 > 0.0 && r1 >= rc1) || (rc2 > 0.0 && r2 >= rc2)) {
        if (d1) *d1 = 0.0;
        if (d2) *d2 = 0.0;
        return 0.0;
    }
    const int want = (d1 != 0) || (d2 != 0);
    double dg = 0.0;
    const double g = p->C * p->D *
                     ugur_shape(p, r1 + r2, p->alpha3, want ? &dg : 0);
    if (want) dg *= p->C * p->D;
    double dS1 = 0.0, dS2 = 0.0;
    const double S1 = ugur_switch(p, r1, rc1, want ? &dS1 : 0);
    const double S2 = ugur_switch(p, r2, rc2, want ? &dS2 : 0);
    /* E = g(r1 + r2) S(r1) S(r2); dg/dr1 = dg/dr2 = g', but the switches make
     * the two derivatives differ, which is the whole reason they are separate
     * arguments here */
    if (d1) *d1 = dg * S1 * S2 + g * dS1 * S2;
    if (d2) *d2 = dg * S1 * S2 + g * S1 * dS2;
    return g * S1 * S2;
}

/* The UG angular factor and its derivative.
 *
 *   h(c) = 1 + lam2 P2(c) + lam4 P4(c),  c = cos(theta)
 *   P2 = (3c^2 - 1)/2,   P4 = (35c^4 - 30c^2 + 3)/8
 *
 * With both weights zero this is 1 and h' is 0, which is what makes the
 * published form a strict special case rather than a separate branch. */
static inline double ugur_h(const ugur_param *p, double c, double *dh)
{
    if (p->lam2 == 0.0 && p->lam4 == 0.0) {
        if (dh) *dh = 0.0;
        return 1.0;
    }
    const double c2 = c * c;
    double h = 1.0 + p->lam2 * 0.5 * (3.0 * c2 - 1.0);
    double d = p->lam2 * 3.0 * c;
    if (p->lam4 != 0.0) {
        h += p->lam4 * 0.125 * (35.0 * c2 * c2 - 30.0 * c2 + 3.0);
        d += p->lam4 * 0.125 * (140.0 * c2 * c - 60.0 * c);
    }
    if (dh) *dh = d;
    return h;
}

/* Full UG triplet: E = g(r1+r2) h(cos) S(r1) S(r2).
 *
 * The angular factor changes the SHAPE of the force, not just its size.  With
 * h = 1 the two legs are pushed along their own unit vectors and nothing acts
 * perpendicular to them - that is why the published form needs no angular
 * derivative anywhere.  As soon as h depends on the angle there are components
 * across the legs, so this returns the radial derivatives and dE/dcos
 * separately and leaves the caller to build the vectors, which is the
 * Stillinger-Weber arrangement and what ev_tally3 expects. */
static inline double ugur_phi3_ang(const ugur_param *p, double r1, double r2,
                                   double c, double rc1, double rc2,
                                   double *d1, double *d2, double *dcos)
{
    if (p->C == 0.0) {
        if (d1) *d1 = 0.0;
        if (d2) *d2 = 0.0;
        if (dcos) *dcos = 0.0;
        return 0.0;
    }
    if ((rc1 > 0.0 && r1 >= rc1) || (rc2 > 0.0 && r2 >= rc2)) {
        if (d1) *d1 = 0.0;
        if (d2) *d2 = 0.0;
        if (dcos) *dcos = 0.0;
        return 0.0;
    }
    const int want = (d1 != 0) || (d2 != 0) || (dcos != 0);
    double dg = 0.0;
    const double g = p->C * p->D *
                     ugur_shape(p, r1 + r2, p->alpha3, want ? &dg : 0);
    if (want) dg *= p->C * p->D;
    double dh = 0.0;
    const double h = ugur_h(p, c, want ? &dh : 0);
    double dS1 = 0.0, dS2 = 0.0;
    const double S1 = ugur_switch(p, r1, rc1, want ? &dS1 : 0);
    const double S2 = ugur_switch(p, r2, rc2, want ? &dS2 : 0);
    if (d1) *d1 = (dg * S1 + g * dS1) * h * S2;
    if (d2) *d2 = (dg * S2 + g * dS2) * h * S1;
    if (dcos) *dcos = g * dh * S1 * S2;
    return g * h * S1 * S2;
}

static inline double ugur_phi3(const ugur_param *p, double r1, double r2,
                               double *d1, double *d2)
{
    return ugur_phi3_ab(p, r1, r2, p->rcut3, p->rcut3, d1, d2);
}

#endif
