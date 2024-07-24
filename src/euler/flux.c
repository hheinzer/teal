#include "flux.h"

#include <math.h>

#include "core/utils.h"
#include "euler.h"
#include "riemann.h"
#include "rotate.h"
#include "update.h"

static void physical_flux(const double *u, double *f);

void godunov(const Equations *eqns, const double *n, const double *g_ul, const double *g_ur,
             double *g_f)
{
    // rotate global to local
    double ul[N_VARS], ur[N_VARS];
    global_to_local(n, g_ul, ul);
    global_to_local(n, g_ur, ur);
    prim_to_cons(eqns, ul);
    prim_to_cons(eqns, ur);

    // Toro 2009, sec. 6.3
    const double gamma = eqns->scalar.value[GAMMA];
    double u[N_VARS];
    riemann(&u[D], &u[U], &u[P], ul[D], ul[U], ul[P], ur[D], ur[U], ur[P], 0, gamma);
    if (u[U] >= 0) {
        u[V] = ul[V];
        u[W] = ul[W];
    }
    else {
        u[V] = ur[V];
        u[W] = ur[W];
    }
    prim_to_cons(eqns, u);

    // compute flux
    double f[N_VARS];
    physical_flux(u, f);
    local_to_global(n, f, g_f);
}

void roe(const Equations *eqns, const double *n, const double *g_ul, const double *g_ur,
         double *g_f)
{
    // rotate global to local
    double ul[N_VARS], ur[N_VARS];
    global_to_local(n, g_ul, ul);
    global_to_local(n, g_ur, ur);
    prim_to_cons(eqns, ul);
    prim_to_cons(eqns, ur);

    // Toro 2009, sec. 11.2.2
    const double gamma = eqns->scalar.value[GAMMA];
    const double Hl = (ul[DE] + ul[P]) / ul[D];
    const double Hr = (ur[DE] + ur[P]) / ur[D];

    // compute Roe averages
    const double sqrt_dl = sqrt(ul[D]);
    const double sqrt_dr = sqrt(ur[D]);
    const double denom = 1 / (sqrt_dl + sqrt_dr);
    const double us = (sqrt_dl * ul[U] + sqrt_dr * ur[U]) * denom;
    const double vs = (sqrt_dl * ul[V] + sqrt_dr * ur[V]) * denom;
    const double ws = (sqrt_dl * ul[W] + sqrt_dr * ur[W]) * denom;
    const double Hs = (sqrt_dl * Hl + sqrt_dr * Hr) * denom;
    const double V2s = sq(us) + sq(vs) + sq(ws);
    const double as = sqrt((gamma - 1) * (Hs - 0.5 * V2s));

    // compute averaged eigenvalues
    double es[3] = {us - as, us, us + as};

    // Harten and Hyman entropy fix, Pelanti et al. 2018, sec. 4.1 and 4.3
    const double al = sqrt(gamma * ul[P] / ul[D]);
    const double ar = sqrt(gamma * ur[P] / ur[D]);
    const double el[3] = {ul[U] - al, ul[U], ul[U] + al};
    const double er[3] = {ur[U] - ar, ur[U], ur[U] + ar};
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
    const double du1 = ur[D] - ul[D];
    const double du2 = ur[DU] - ul[DU];
    const double du3 = ur[DV] - ul[DV];
    const double du4 = ur[DW] - ul[DW];
    const double du5 = ur[DE] - ul[DE] - (du3 - vs * du1) * vs - (du4 - ws * du1) * ws;
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
    double fl[N_VARS], fr[N_VARS], f[N_VARS];
    physical_flux(ul, fl);
    physical_flux(ur, fr);
    f[D] = 0.5 * (fl[D] + fr[D] - df[0]);
    f[DU] = 0.5 * (fl[DU] + fr[DU] - df[1]);
    f[DV] = 0.5 * (fl[DV] + fr[DV] - df[2]);
    f[DW] = 0.5 * (fl[DW] + fr[DW] - df[3]);
    f[DE] = 0.5 * (fl[DE] + fr[DE] - df[4]);
    local_to_global(n, f, g_f);
}

void hll(const Equations *eqns, const double *n, const double *g_ul, const double *g_ur,
         double *g_f)
{
    // rotate global to local
    double ul[N_VARS], ur[N_VARS];
    global_to_local(n, g_ul, ul);
    global_to_local(n, g_ur, ur);
    prim_to_cons(eqns, ul);
    prim_to_cons(eqns, ur);

    // Toro 2009, sec. 10.3
    const double gamma = eqns->scalar.value[GAMMA];
    const double Hl = (ul[DE] + ul[P]) / ul[D];
    const double Hr = (ur[DE] + ur[P]) / ur[D];

    // compute Roe averages
    const double sqrt_dl = sqrt(ul[D]);
    const double sqrt_dr = sqrt(ur[D]);
    const double denom = 1 / (sqrt_dl + sqrt_dr);
    const double us = (sqrt_dl * ul[U] + sqrt_dr * ur[U]) * denom;
    const double vs = (sqrt_dl * ul[V] + sqrt_dr * ur[V]) * denom;
    const double ws = (sqrt_dl * ul[W] + sqrt_dr * ur[W]) * denom;
    const double Hs = (sqrt_dl * Hl + sqrt_dr * Hr) * denom;
    const double V2s = sq(us) + sq(vs) + sq(ws);
    const double as = sqrt((gamma - 1) * (Hs - 0.5 * V2s));

    // compute signal speeds
    const double al = sqrt(gamma * ul[P] / ul[D]);
    const double ar = sqrt(gamma * ur[P] / ur[D]);
    const double sl = min(ul[U] - al, us - as);
    const double sr = max(ur[U] + ar, us + as);

    // compute flux
    double fl[N_VARS], fr[N_VARS], f[N_VARS];
    if (0 <= sl)
        physical_flux(ul, f);
    else if (0 >= sr)
        physical_flux(ur, f);
    else {
        physical_flux(ul, fl);
        physical_flux(ur, fr);
        f[D] = (sr * fl[D] - sl * fr[D] + sl * sr * (ur[D] - ul[D])) / (sr - sl);
        f[DU] = (sr * fl[DU] - sl * fr[DU] + sl * sr * (ur[DU] - ul[DU])) / (sr - sl);
        f[DV] = (sr * fl[DV] - sl * fr[DV] + sl * sr * (ur[DV] - ul[DV])) / (sr - sl);
        f[DW] = (sr * fl[DW] - sl * fr[DW] + sl * sr * (ur[DW] - ul[DW])) / (sr - sl);
        f[DE] = (sr * fl[DE] - sl * fr[DE] + sl * sr * (ur[DE] - ul[DE])) / (sr - sl);
    }
    local_to_global(n, f, g_f);
}

void hllc(const Equations *eqns, const double *n, const double *g_ul, const double *g_ur,
          double *g_f)
{
    // rotate global to local
    double ul[N_VARS], ur[N_VARS];
    global_to_local(n, g_ul, ul);
    global_to_local(n, g_ur, ur);
    prim_to_cons(eqns, ul);
    prim_to_cons(eqns, ur);

    // Toro 2009, sec. 10.4.2
    const double gamma = eqns->scalar.value[GAMMA];
    const double Hl = (ul[DE] + ul[P]) / ul[D];
    const double Hr = (ur[DE] + ur[P]) / ur[D];

    // compute Roe averages
    const double sqrt_dl = sqrt(ul[D]);
    const double sqrt_dr = sqrt(ur[D]);
    const double denom = 1 / (sqrt_dl + sqrt_dr);
    const double us = (sqrt_dl * ul[U] + sqrt_dr * ur[U]) * denom;
    const double vs = (sqrt_dl * ul[V] + sqrt_dr * ur[V]) * denom;
    const double ws = (sqrt_dl * ul[W] + sqrt_dr * ur[W]) * denom;
    const double Hs = (sqrt_dl * Hl + sqrt_dr * Hr) * denom;
    const double V2s = sq(us) + sq(vs) + sq(ws);
    const double as = sqrt((gamma - 1) * (Hs - 0.5 * V2s));

    // compute signal speeds
    const double al = sqrt(gamma * ul[P] / ul[D]);
    const double ar = sqrt(gamma * ur[P] / ur[D]);
    const double sl = min(ul[U] - al, us - as);
    const double sr = max(ur[U] + ar, us + as);

    // compute flux
    double f[N_VARS];
    if (0 <= sl)
        physical_flux(ul, f);
    else if (0 >= sr)
        physical_flux(ur, f);
    else {
        const double ss = (ur[P] - ul[P] + ul[DU] * (sl - ul[U]) - ur[DU] * (sr - ur[U])) /
                          (ul[D] * (sl - ul[U]) - ur[D] * (sr - ur[U]));
        if (sl <= 0 && 0 <= ss) {
            const double fac = ul[D] * (sl - ul[U]) / (sl - ss);
            const double usl[5] = {
                fac,
                fac * ss,
                fac * ul[V],
                fac * ul[W],
                fac * (ul[DE] / ul[D] + (ss - ul[U]) * (ss + ul[P] / (ul[D] * (sl - ul[U])))),
            };
            physical_flux(ul, f);
            f[D] += sl * (usl[0] - ul[D]);
            f[DU] += sl * (usl[1] - ul[DU]);
            f[DV] += sl * (usl[2] - ul[DV]);
            f[DW] += sl * (usl[3] - ul[DW]);
            f[DE] += sl * (usl[4] - ul[DE]);
        }
        else {
            const double fac = ur[D] * (sr - ur[U]) / (sr - ss);
            const double usr[5] = {
                fac,
                fac * ss,
                fac * ur[V],
                fac * ur[W],
                fac * (ur[DE] / ur[D] + (ss - ur[U]) * (ss + ur[P] / (ur[D] * (sr - ur[U])))),
            };
            physical_flux(ur, f);
            f[D] += sr * (usr[0] - ur[D]);
            f[DU] += sr * (usr[1] - ur[DU]);
            f[DV] += sr * (usr[2] - ur[DV]);
            f[DW] += sr * (usr[3] - ur[DW]);
            f[DE] += sr * (usr[4] - ur[DE]);
        }
    }
    local_to_global(n, f, g_f);
}

void hlle(const Equations *eqns, const double *n, const double *g_ul, const double *g_ur,
          double *g_f)
{
    // rotate global to local
    double ul[N_VARS], ur[N_VARS];
    global_to_local(n, g_ul, ul);
    global_to_local(n, g_ur, ur);
    prim_to_cons(eqns, ul);
    prim_to_cons(eqns, ur);

    // Einfeldt 1988
    const double gamma = eqns->scalar.value[GAMMA];

    // compute Roe averages
    const double sqrt_dl = sqrt(ul[D]);
    const double sqrt_dr = sqrt(ur[D]);
    const double denom = 1 / (sqrt_dl + sqrt_dr);
    const double us = (sqrt_dl * ul[U] + sqrt_dr * ur[U]) * denom;

    // compute signal speeds
    const double al = sqrt(gamma * ul[P] / ul[D]);
    const double ar = sqrt(gamma * ur[P] / ur[D]);
    const double d =
        sqrt((sqrt(ul[D]) * al * al + sqrt(ur[D]) * ar * ar) * denom +
             0.5 * (sqrt(ul[D]) * sqrt(ur[D])) * sq(denom) * (ur[U] - ul[U]) * (ur[U] - ul[U]));
    const double sl = min(ul[U] - al, us - d);
    const double sr = max(ur[U] + ar, us + d);

    // compute flux
    double fl[N_VARS], fr[N_VARS], f[N_VARS];
    if (0 <= sl)
        physical_flux(ul, f);
    else if (0 >= sr)
        physical_flux(ur, f);
    else {
        physical_flux(ul, fl);
        physical_flux(ur, fr);
        f[D] = (sr * fl[D] - sl * fr[D] + sl * sr * (ur[D] - ul[D])) / (sr - sl);
        f[DU] = (sr * fl[DU] - sl * fr[DU] + sl * sr * (ur[DU] - ul[DU])) / (sr - sl);
        f[DV] = (sr * fl[DV] - sl * fr[DV] + sl * sr * (ur[DV] - ul[DV])) / (sr - sl);
        f[DW] = (sr * fl[DW] - sl * fr[DW] + sl * sr * (ur[DW] - ul[DW])) / (sr - sl);
        f[DE] = (sr * fl[DE] - sl * fr[DE] + sl * sr * (ur[DE] - ul[DE])) / (sr - sl);
    }
    local_to_global(n, f, g_f);
}

void lxf(const Equations *eqns, const double *n, const double *g_ul, const double *g_ur,
         double *g_f)
{
    // rotate global to local
    double ul[N_VARS], ur[N_VARS];
    global_to_local(n, g_ul, ul);
    global_to_local(n, g_ur, ur);
    prim_to_cons(eqns, ul);
    prim_to_cons(eqns, ur);

    // Lax-Friedrichs flux
    const double gamma = eqns->scalar.value[GAMMA];
    const double cl = sqrt(gamma * ul[P] / ul[D]);
    const double cr = sqrt(gamma * ur[P] / ur[D]);
    const double s = max(fabs(ul[U]) + cl, fabs(ur[U]) + cr);

    // compute flux
    double fl[N_VARS], fr[N_VARS], f[N_VARS];
    physical_flux(ul, fl);
    physical_flux(ur, fr);
    f[D] = 0.5 * (fl[D] + fr[D] - s * (ur[D] - ul[D]));
    f[DU] = 0.5 * (fl[DU] + fr[DU] - s * (ur[DU] - ul[DU]));
    f[DV] = 0.5 * (fl[DV] + fr[DV] - s * (ur[DV] - ul[DV]));
    f[DW] = 0.5 * (fl[DW] + fr[DW] - s * (ur[DW] - ul[DW]));
    f[DE] = 0.5 * (fl[DE] + fr[DE] - s * (ur[DE] - ul[DE]));
    local_to_global(n, f, g_f);
}

static void physical_flux(const double *u, double *f)
{
    f[D] = u[DU];
    f[DU] = u[DU] * u[U] + u[P];
    f[DV] = u[DU] * u[V];
    f[DW] = u[DU] * u[W];
    f[DE] = (u[DE] + u[P]) * u[U];
}
