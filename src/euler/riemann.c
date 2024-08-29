#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "euler.h"
#include "teal/utils.h"

#define TOLPRE 1e-6
#define NRITER 20
#define QUSER 2

static void starpu(double *pm, double *um, double dl, double ul, double pl, double dr, double ur,
                   double pr, double cl, double cr, const double *g);

static double guessp(double dl, double ul, double pl, double dr, double ur, double pr, double cl,
                     double cr, const double *g);

static void prefun(double *f, double *fd, double p, double dk, double pk, double ck,
                   const double *g);

static void sample(double *d, double *u, double *p, double dl, double ul, double pl, double dr,
                   double ur, double pr, double cl, double cr, double pm, double um, double s,
                   const double *g);

void euler_riemann(double *d, double *u, double *p, double dl, double ul, double pl, double dr,
                   double ur, double pr, double s, double gamma)
{
    // Toro 1999, sec. 4.9
    const double g[9] = {
        gamma,
        (gamma - 1) / (2 * gamma),
        (gamma + 1) / (2 * gamma),
        2 * gamma / (gamma - 1),
        2 / (gamma - 1),
        2 / (gamma + 1),
        (gamma - 1) / (gamma + 1),
        (gamma - 1) / 2,
        gamma - 1,
    };
    const double cl = sqrt(gamma * pl / dl);
    const double cr = sqrt(gamma * pr / dr);
    assert(g[4] * (cl + cr) > (ur - ul));

    double pm, um;
    starpu(&pm, &um, dl, ul, pl, dr, ur, pr, cl, cr, g);
    sample(d, u, p, dl, ul, pl, dr, ur, pr, cl, cr, pm, um, s, g);
}

static void starpu(double *pm, double *um, double dl, double ul, double pl, double dr, double ur,
                   double pr, double cl, double cr, const double *g)
{
    const double udiff = ur - ul;
    double pold = guessp(dl, ul, pl, dr, ur, pr, cl, cr, g);
    double fl, fld, fr, frd;
    for (long iter = 0; iter < NRITER; ++iter) {
        prefun(&fl, &fld, pold, dl, pl, cl, g);
        prefun(&fr, &frd, pold, dr, pr, cr, g);
        *pm = pold - (fl + fr + udiff) / (fld + frd);
        if (2 * fabs((*pm - pold) / (*pm + pold)) <= TOLPRE) {
            *um = 0.5 * (ul + ur + fr - fl);
            return;
        }
        pold = *pm = max(TOLPRE, *pm);
    }
    abort();
}

static double guessp(double dl, double ul, double pl, double dr, double ur, double pr, double cl,
                     double cr, const double *g)
{
    const double cup = 0.25 * (dl + dr) * (cl + cr);
    const double ppv = max(0.5 * (pl + pr) + 0.5 * (ul - ur) * cup, 0);
    const double pmin = min(pl, pr);
    const double pmax = max(pl, pr);
    const double qmax = pmax / pmin;
    if (qmax <= QUSER && pmin <= ppv && ppv <= pmax) return ppv;
    if (ppv < pmin) {
        const double pq = pow(pl / pr, g[1]);
        const double um = (pq * ul / cl + ur / cr + g[4] * (pq - 1)) / (pq / cl + 1 / cr);
        const double ptl = 1 + g[7] * (ul - um) / cl;
        const double ptr = 1 + g[7] * (um - ur) / cr;
        return 0.5 * (pl * pow(ptl, g[3]) + pr * pow(ptr, g[3]));
    }
    const double gel = sqrt((g[5] / dl) / (g[6] * pl + ppv));
    const double ger = sqrt((g[5] / dr) / (g[6] * pr + ppv));
    return (gel * pl + ger * pr - (ur - ul)) / (gel + ger);
}

static void prefun(double *f, double *fd, double p, double dk, double pk, double ck,
                   const double *g)
{
    if (p <= pk) {
        const double prat = p / pk;
        *f = g[4] * ck * (pow(prat, g[1]) - 1);
        *fd = (1 / (dk * ck)) * pow(prat, -g[2]);
    }
    else {
        const double ak = g[5] / dk;
        const double bk = g[6] * pk;
        const double qrt = sqrt(ak / (bk + p));
        *f = (p - pk) * qrt;
        *fd = (1 - 0.5 * (p - pk) / (bk + p)) * qrt;
    }
}

static void sample(double *d, double *u, double *p, double dl, double ul, double pl, double dr,
                   double ur, double pr, double cl, double cr, double pm, double um, double s,
                   const double *g)
{
    if (s <= um) {
        if (pm <= pl) {
            const double shl = ul - cl;
            if (s <= shl) {
                *d = dl;
                *u = ul;
                *p = pl;
            }
            else {
                const double cml = cl * pow(pm / pl, g[1]);
                const double stl = um - cml;
                if (s > stl) {
                    *d = dl * pow(pm / pl, 1 / g[0]);
                    *u = um;
                    *p = pm;
                }
                else {
                    const double c = g[5] * (cl + g[7] * (ul - s));
                    *d = dl * pow(c / cl, g[4]);
                    *u = g[5] * (cl + g[7] * ul + s);
                    *p = pl * pow(c / cl, g[3]);
                }
            }
        }
        else {
            const double pml = pm / pl;
            const double sl = ul - cl * sqrt(g[2] * pml + g[1]);
            if (s <= sl) {
                *d = dl;
                *u = ul;
                *p = pl;
            }
            else {
                *d = dl * (pml + g[6]) / (pml * g[6] + 1);
                *u = um;
                *p = pm;
            }
        }
    }
    else {
        if (pm > pr) {
            const double pmr = pm / pr;
            const double sr = ur + cr * sqrt(g[2] * pmr + g[1]);
            if (s >= sr) {
                *d = dr;
                *u = ur;
                *p = pr;
            }
            else {
                *d = dr * (pmr + g[6]) / (pmr * g[6] + 1);
                *u = um;
                *p = pm;
            }
        }
        else {
            const double shr = ur + cr;
            if (s >= shr) {
                *d = dr;
                *u = ur;
                *p = pr;
            }
            else {
                const double cmr = cr * pow(pm / pr, g[1]);
                const double str = um + cmr;
                if (s <= str) {
                    *d = dr * pow(pm / pr, 1 / g[0]);
                    *u = um;
                    *p = pm;
                }
                else {
                    const double c = g[5] * (cr - g[7] * (ur - s));
                    *d = dr * pow(c / cr, g[4]);
                    *u = g[5] * (-cr + g[7] * ur + s);
                    *p = pr * pow(c / cr, g[3]);
                }
            }
        }
    }
}
