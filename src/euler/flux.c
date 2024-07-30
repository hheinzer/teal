#include "flux.h"

#include <math.h>

#include "core/utils.h"
#include "euler.h"
#include "riemann.h"
#include "rotate.h"
#include "update.h"

static void physical_flux(const double *u, double *f);

void godunov(const Equations *eqns, const double *n, const double *ul, const double *ur, double *f)
{
    // rotate global to local
    double lul[N_VARS], lur[N_VARS];
    global_to_local(n, ul, lul);
    global_to_local(n, ur, lur);
    prim_to_cons(eqns, lul);
    prim_to_cons(eqns, lur);

    // Toro 2009, sec. 6.3
    const double gamma = eqns->scalar.value[GAMMA];
    double u[N_VARS];
    riemann(&u[D], &u[U], &u[P], lul[D], lul[U], lul[P], lur[D], lur[U], lur[P], 0, gamma);
    if (u[U] >= 0) {
        u[V] = lul[V];
        u[W] = lul[W];
    }
    else {
        u[V] = lur[V];
        u[W] = lur[W];
    }
    prim_to_cons(eqns, u);

    // compute flux
    double lf[N_VARS];
    physical_flux(u, lf);
    local_to_global(n, lf, f);
}

void roe(const Equations *eqns, const double *n, const double *ul, const double *ur, double *f)
{
    // rotate global to local
    double lul[N_VARS], lur[N_VARS];
    global_to_local(n, ul, lul);
    global_to_local(n, ur, lur);
    prim_to_cons(eqns, lul);
    prim_to_cons(eqns, lur);

    // Toro 2009, sec. 11.2.2
    const double gamma = eqns->scalar.value[GAMMA];
    const double Hl = (lul[DE] + lul[P]) / lul[D];
    const double Hr = (lur[DE] + lur[P]) / lur[D];

    // compute Roe averages
    const double sqrt_dl = sqrt(lul[D]);
    const double sqrt_dr = sqrt(lur[D]);
    const double denom = 1 / (sqrt_dl + sqrt_dr);
    const double us = (sqrt_dl * lul[U] + sqrt_dr * lur[U]) * denom;
    const double vs = (sqrt_dl * lul[V] + sqrt_dr * lur[V]) * denom;
    const double ws = (sqrt_dl * lul[W] + sqrt_dr * lur[W]) * denom;
    const double Hs = (sqrt_dl * Hl + sqrt_dr * Hr) * denom;
    const double V2s = sq(us) + sq(vs) + sq(ws);
    const double as = sqrt((gamma - 1) * (Hs - 0.5 * V2s));

    // compute averaged eigenvalues
    double es[3] = {us - as, us, us + as};

    // Harten and Hyman entropy fix, Pelanti et al. 2018, sec. 4.1 and 4.3
    const double al = sqrt(gamma * lul[P] / lul[D]);
    const double ar = sqrt(gamma * lur[P] / lur[D]);
    const double el[3] = {lul[U] - al, lul[U], lul[U] + al};
    const double er[3] = {lur[U] - ar, lur[U], lur[U] + ar};
    for (long i = 0; i < 3; ++i) {
        const double delta = max(EPS, max(es[i] - el[i], er[i] - es[i]));
        if (fabs(es[i]) < delta) {
            es[i] = delta;  // first
            // es[i] = 0.5 * (SQ(es[i]) / delta + delta);  // second
        }
        else {
            es[i] = fabs(es[i]);
        }
    }

    // compute averaged right eigenvectors
    const double Ks[5][5] = {
        {1, us - as, vs, ws, Hs - us * as},
        {1, us, vs, ws, 0.5 * V2s},
        {0, 0, 1, 0, vs},
        {0, 0, 0, 1, ws},
        {1, us + as, vs, ws, Hs + us * as},
    };

    // compute wave strengths
    const double du1 = lur[D] - lul[D];
    const double du2 = lur[DU] - lul[DU];
    const double du3 = lur[DV] - lul[DV];
    const double du4 = lur[DW] - lul[DW];
    const double du5 = lur[DE] - lul[DE] - (du3 - vs * du1) * vs - (du4 - ws * du1) * ws;
    const double alpha3 = du3 - vs * du1;
    const double alpha4 = du4 - ws * du1;
    const double alpha2 = (gamma - 1) / sq(as) * (du1 * (Hs - sq(us)) + us * du2 - du5);
    const double alpha1 = 0.5 / as * (du1 * (us + as) - du2 - as * alpha2);
    const double alpha5 = du1 - (alpha1 + alpha2);
    const double alpha[5] = {alpha1, alpha2, alpha3, alpha4, alpha5};

    // compute flux correction
    const double lambda[5] = {es[0], es[1], es[1], es[1], es[2]};
    double df[5] = {0};
    for (long i = 0; i < 5; ++i)
        for (long j = 0; j < 5; ++j) df[j] += alpha[i] * lambda[i] * Ks[i][j];

    // compute flux
    double fl[N_VARS], fr[N_VARS], lf[N_VARS];
    physical_flux(lul, fl);
    physical_flux(lur, fr);
    lf[D] = 0.5 * (fl[D] + fr[D] - df[0]);
    lf[DU] = 0.5 * (fl[DU] + fr[DU] - df[1]);
    lf[DV] = 0.5 * (fl[DV] + fr[DV] - df[2]);
    lf[DW] = 0.5 * (fl[DW] + fr[DW] - df[3]);
    lf[DE] = 0.5 * (fl[DE] + fr[DE] - df[4]);
    local_to_global(n, lf, f);
}

void hll(const Equations *eqns, const double *n, const double *ul, const double *ur, double *f)
{
    // rotate global to local
    double lul[N_VARS], lur[N_VARS];
    global_to_local(n, ul, lul);
    global_to_local(n, ur, lur);
    prim_to_cons(eqns, lul);
    prim_to_cons(eqns, lur);

    // Toro 2009, sec. 10.3
    const double gamma = eqns->scalar.value[GAMMA];
    const double Hl = (lul[DE] + lul[P]) / lul[D];
    const double Hr = (lur[DE] + lur[P]) / lur[D];

    // compute Roe averages
    const double sqrt_dl = sqrt(lul[D]);
    const double sqrt_dr = sqrt(lur[D]);
    const double denom = 1 / (sqrt_dl + sqrt_dr);
    const double us = (sqrt_dl * lul[U] + sqrt_dr * lur[U]) * denom;
    const double vs = (sqrt_dl * lul[V] + sqrt_dr * lur[V]) * denom;
    const double ws = (sqrt_dl * lul[W] + sqrt_dr * lur[W]) * denom;
    const double Hs = (sqrt_dl * Hl + sqrt_dr * Hr) * denom;
    const double V2s = sq(us) + sq(vs) + sq(ws);
    const double as = sqrt((gamma - 1) * (Hs - 0.5 * V2s));

    // compute signal speeds
    const double al = sqrt(gamma * lul[P] / lul[D]);
    const double ar = sqrt(gamma * lur[P] / lur[D]);
    const double sl = min(lul[U] - al, us - as);
    const double sr = max(lur[U] + ar, us + as);

    // compute flux
    double fl[N_VARS], fr[N_VARS], lf[N_VARS];
    if (0 <= sl)
        physical_flux(lul, lf);
    else if (0 >= sr)
        physical_flux(lur, lf);
    else {
        physical_flux(lul, fl);
        physical_flux(lur, fr);
        lf[D] = (sr * fl[D] - sl * fr[D] + sl * sr * (lur[D] - lul[D])) / (sr - sl);
        lf[DU] = (sr * fl[DU] - sl * fr[DU] + sl * sr * (lur[DU] - lul[DU])) / (sr - sl);
        lf[DV] = (sr * fl[DV] - sl * fr[DV] + sl * sr * (lur[DV] - lul[DV])) / (sr - sl);
        lf[DW] = (sr * fl[DW] - sl * fr[DW] + sl * sr * (lur[DW] - lul[DW])) / (sr - sl);
        lf[DE] = (sr * fl[DE] - sl * fr[DE] + sl * sr * (lur[DE] - lul[DE])) / (sr - sl);
    }
    local_to_global(n, lf, f);
}

void hllc(const Equations *eqns, const double *n, const double *ul, const double *ur, double *f)
{
    // rotate global to local
    double lul[N_VARS], lur[N_VARS];
    global_to_local(n, ul, lul);
    global_to_local(n, ur, lur);
    prim_to_cons(eqns, lul);
    prim_to_cons(eqns, lur);

    // Toro 2009, sec. 10.4.2
    const double gamma = eqns->scalar.value[GAMMA];
    const double Hl = (lul[DE] + lul[P]) / lul[D];
    const double Hr = (lur[DE] + lur[P]) / lur[D];

    // compute Roe averages
    const double sqrt_dl = sqrt(lul[D]);
    const double sqrt_dr = sqrt(lur[D]);
    const double denom = 1 / (sqrt_dl + sqrt_dr);
    const double us = (sqrt_dl * lul[U] + sqrt_dr * lur[U]) * denom;
    const double vs = (sqrt_dl * lul[V] + sqrt_dr * lur[V]) * denom;
    const double ws = (sqrt_dl * lul[W] + sqrt_dr * lur[W]) * denom;
    const double Hs = (sqrt_dl * Hl + sqrt_dr * Hr) * denom;
    const double V2s = sq(us) + sq(vs) + sq(ws);
    const double as = sqrt((gamma - 1) * (Hs - 0.5 * V2s));

    // compute signal speeds
    const double al = sqrt(gamma * lul[P] / lul[D]);
    const double ar = sqrt(gamma * lur[P] / lur[D]);
    const double sl = min(lul[U] - al, us - as);
    const double sr = max(lur[U] + ar, us + as);

    // compute flux
    double lf[N_VARS];
    if (0 <= sl)
        physical_flux(lul, lf);
    else if (0 >= sr)
        physical_flux(lur, lf);
    else {
        const double ss = (lur[P] - lul[P] + lul[DU] * (sl - lul[U]) - lur[DU] * (sr - lur[U])) /
                          (lul[D] * (sl - lul[U]) - lur[D] * (sr - lur[U]));
        if (sl <= 0 && 0 <= ss) {
            const double fac = lul[D] * (sl - lul[U]) / (sl - ss);
            const double usl[5] = {
                fac,
                fac * ss,
                fac * lul[V],
                fac * lul[W],
                fac * (lul[DE] / lul[D] + (ss - lul[U]) * (ss + lul[P] / (lul[D] * (sl - lul[U])))),
            };
            physical_flux(lul, lf);
            lf[D] += sl * (usl[0] - lul[D]);
            lf[DU] += sl * (usl[1] - lul[DU]);
            lf[DV] += sl * (usl[2] - lul[DV]);
            lf[DW] += sl * (usl[3] - lul[DW]);
            lf[DE] += sl * (usl[4] - lul[DE]);
        }
        else {
            const double fac = lur[D] * (sr - lur[U]) / (sr - ss);
            const double usr[5] = {
                fac,
                fac * ss,
                fac * lur[V],
                fac * lur[W],
                fac * (lur[DE] / lur[D] + (ss - lur[U]) * (ss + lur[P] / (lur[D] * (sr - lur[U])))),
            };
            physical_flux(lur, lf);
            lf[D] += sr * (usr[0] - lur[D]);
            lf[DU] += sr * (usr[1] - lur[DU]);
            lf[DV] += sr * (usr[2] - lur[DV]);
            lf[DW] += sr * (usr[3] - lur[DW]);
            lf[DE] += sr * (usr[4] - lur[DE]);
        }
    }
    local_to_global(n, lf, f);
}

void hlle(const Equations *eqns, const double *n, const double *ul, const double *ur, double *f)
{
    // rotate global to local
    double lul[N_VARS], lur[N_VARS];
    global_to_local(n, ul, lul);
    global_to_local(n, ur, lur);
    prim_to_cons(eqns, lul);
    prim_to_cons(eqns, lur);

    // Einfeldt 1988
    const double gamma = eqns->scalar.value[GAMMA];

    // compute Roe averages
    const double sqrt_dl = sqrt(lul[D]);
    const double sqrt_dr = sqrt(lur[D]);
    const double denom = 1 / (sqrt_dl + sqrt_dr);
    const double us = (sqrt_dl * lul[U] + sqrt_dr * lur[U]) * denom;

    // compute signal speeds
    const double al = sqrt(gamma * lul[P] / lul[D]);
    const double ar = sqrt(gamma * lur[P] / lur[D]);
    const double d = sqrt((sqrt(lul[D]) * al * al + sqrt(lur[D]) * ar * ar) * denom +
                          0.5 * (sqrt(lul[D]) * sqrt(lur[D])) * sq(denom) * (lur[U] - lul[U]) *
                              (lur[U] - lul[U]));
    const double sl = min(lul[U] - al, us - d);
    const double sr = max(lur[U] + ar, us + d);

    // compute flux
    double fl[N_VARS], fr[N_VARS], lf[N_VARS];
    if (0 <= sl)
        physical_flux(lul, lf);
    else if (0 >= sr)
        physical_flux(lur, lf);
    else {
        physical_flux(lul, fl);
        physical_flux(lur, fr);
        lf[D] = (sr * fl[D] - sl * fr[D] + sl * sr * (lur[D] - lul[D])) / (sr - sl);
        lf[DU] = (sr * fl[DU] - sl * fr[DU] + sl * sr * (lur[DU] - lul[DU])) / (sr - sl);
        lf[DV] = (sr * fl[DV] - sl * fr[DV] + sl * sr * (lur[DV] - lul[DV])) / (sr - sl);
        lf[DW] = (sr * fl[DW] - sl * fr[DW] + sl * sr * (lur[DW] - lul[DW])) / (sr - sl);
        lf[DE] = (sr * fl[DE] - sl * fr[DE] + sl * sr * (lur[DE] - lul[DE])) / (sr - sl);
    }
    local_to_global(n, lf, f);
}

void lxf(const Equations *eqns, const double *n, const double *ul, const double *ur, double *f)
{
    // rotate global to local
    double lul[N_VARS], lur[N_VARS];
    global_to_local(n, ul, lul);
    global_to_local(n, ur, lur);
    prim_to_cons(eqns, lul);
    prim_to_cons(eqns, lur);

    // Lax-Friedrichs flux
    const double gamma = eqns->scalar.value[GAMMA];
    const double cl = sqrt(gamma * lul[P] / lul[D]);
    const double cr = sqrt(gamma * lur[P] / lur[D]);
    const double s = max(fabs(lul[U]) + cl, fabs(lur[U]) + cr);

    // compute flux
    double fl[N_VARS], fr[N_VARS], lf[N_VARS];
    physical_flux(lul, fl);
    physical_flux(lur, fr);
    lf[D] = 0.5 * (fl[D] + fr[D] - s * (lur[D] - lul[D]));
    lf[DU] = 0.5 * (fl[DU] + fr[DU] - s * (lur[DU] - lul[DU]));
    lf[DV] = 0.5 * (fl[DV] + fr[DV] - s * (lur[DV] - lul[DV]));
    lf[DW] = 0.5 * (fl[DW] + fr[DW] - s * (lur[DW] - lul[DW]));
    lf[DE] = 0.5 * (fl[DE] + fr[DE] - s * (lur[DE] - lul[DE]));
    local_to_global(n, lf, f);
}

static void physical_flux(const double *u, double *f)
{
    f[D] = u[DU];
    f[DU] = u[DU] * u[U] + u[P];
    f[DV] = u[DU] * u[V];
    f[DW] = u[DU] * u[W];
    f[DE] = (u[DE] + u[P]) * u[U];
}
