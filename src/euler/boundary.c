#include "boundary.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "euler.h"
#include "rotate.h"
#include "teal/utils.h"
#include "update.h"

static ApplyBC symmetry, supersonic_inflow, supersonic_outflow, subsonic_inflow, subsonic_outflow,
    farfield;

ApplyBC *euler_bc(const char *name)
{
    if (!strcmp(name, "symmetry") || !strcmp(name, "slipwall")) return symmetry;
    if (!strcmp(name, "supersonic inflow")) return supersonic_inflow;
    if (!strcmp(name, "supersonic outflow")) return supersonic_outflow;
    if (!strcmp(name, "subsonic inflow")) return subsonic_inflow;
    if (!strcmp(name, "subsonic outflow")) return subsonic_outflow;
    if (!strcmp(name, "farfield")) return farfield;
    abort();
}

static void symmetry(const Equations *, const Matrix3d b, const double *, const double *ui,
                     double *ug)
{
    // Blazek 2015, sec. 8.6
    const double *n = b[X];
    const double vn = ui[U] * n[X] + ui[V] * n[Y] + ui[W] * n[Z];
    ug[D] = ui[D];
    ug[U] = ui[U] - 2 * vn * n[X];
    ug[V] = ui[V] - 2 * vn * n[Y];
    ug[W] = ui[W] - 2 * vn * n[Z];
    ug[P] = ui[P];
}

static void supersonic_inflow(const Equations *, const Matrix3d, const double *state,
                              const double *, double *ug)
{
    // Blazek 2015, sec. 8.3.1
    ug[D] = state[D];
    ug[U] = state[U];
    ug[V] = state[V];
    ug[W] = state[W];
    ug[P] = state[P];
}

static void supersonic_outflow(const Equations *, const Matrix3d, const double *, const double *ui,
                               double *ug)
{
    // Blazek 2015, sec. 8.3.1
    ug[D] = ui[D];
    ug[U] = ui[U];
    ug[V] = ui[V];
    ug[W] = ui[W];
    ug[P] = ui[P];
}

static void subsonic_inflow(const Equations *eqns, const Matrix3d b, const double *state,
                            const double *ui, double *ug)
{
    // Blazek 2015, sec. 8.3.1
    const double *n = b[X];
    const double gamma = eqns->scalar.value[GAMMA];
    const double d0 = ui[D];
    const double c0 = sqrt(gamma * ui[P] / ui[D]);
    ug[P] =
        0.5 *
        (state[P] + ui[P] -
         d0 * c0 *
             (n[X] * (state[U] - ui[U]) + n[Y] * (state[V] - ui[V]) + n[Z] * (state[W] - ui[W])));
    ug[D] = state[D] + (ug[P] - state[P]) / sq(c0);
    ug[U] = state[U] - n[X] * (state[P] - ug[P]) / (d0 * c0);
    ug[V] = state[V] - n[Y] * (state[P] - ug[P]) / (d0 * c0);
    ug[W] = state[W] - n[Z] * (state[P] - ug[P]) / (d0 * c0);
}

static void subsonic_outflow(const Equations *eqns, const Matrix3d b, const double *state,
                             const double *ui, double *ug)
{
    // Blazek 2015, sec. 8.3.1
    const double *n = b[X];
    const double gamma = eqns->scalar.value[GAMMA];
    const double d0 = ui[D];
    const double c0 = sqrt(gamma * ui[P] / ui[D]);
    ug[P] = state[P];
    ug[D] = ui[D] + (ug[P] - ui[P]) / sq(c0);
    ug[U] = ui[U] + n[X] * (ui[P] - ug[P]) / (d0 * c0);
    ug[V] = ui[V] + n[Y] * (ui[P] - ug[P]) / (d0 * c0);
    ug[W] = ui[W] + n[Z] * (ui[P] - ug[P]) / (d0 * c0);
}

static void farfield(const Equations *eqns, const Matrix3d b, const double *state, const double *ui,
                     double *ug)
{
    // compute initial ghost cell state
    ug[D] = state[D];
    ug[U] = state[U];
    ug[V] = state[V];
    ug[W] = state[W];
    ug[P] = state[P];

    // rotate variables
    double lui[N_VARS], lug[N_VARS];
    euler_global_to_local(b, ui, lui);
    euler_global_to_local(b, ug, lug);
    euler_prim_to_cons(eqns, lui);
    euler_prim_to_cons(eqns, lug);

    // compute characteristic variables
    const double gamma = eqns->scalar.value[GAMMA];
    const double g = gamma - 1;
    const double a = sqrt(gamma * lui[P] / lui[D]);
    const double u = lui[U];
    const double v = lui[V];
    const double w = lui[W];
    const double V2 = sq(u) + sq(v) + sq(w);
    const double H = 0.5 * V2 + sq(a) / g;
    const double fac = g / (2 * sq(a));
    const double invK[5][5] = {
        {H + a / g * (u - a), -(u + a / g), -v, -w, 1},
        {-2 * H + 4 / g * sq(a), 2 * u, 2 * v, 2 * w, -2},
        {-2 * v * sq(a) / g, 0, 2 * sq(a) / g, 0, 0},
        {-2 * w * sq(a) / g, 0, 0, 2 * sq(a) / g, 0},
        {H - a / g * (u + a), -u + a / g, -v, -w, 1},
    };
    const long cons[5] = {D, DU, DV, DW, DE};
    double ei[5] = {0}, eg[5] = {0};
    for (long i = 0; i < 5; ++i) {
        for (long j = 0; j < 5; ++j) {
            ei[i] += invK[i][j] * lui[cons[j]];
            eg[i] += invK[i][j] * lug[cons[j]];
        }
        ei[i] *= fac;
        eg[i] *= fac;
    }

    // compute characteristic variables of ghost cell
    const double ag = sqrt(gamma * lug[P] / lug[D]);
    const double vg = lug[U];
    if (vg - ag > 0) eg[0] = ei[0];
    if (vg > 0) {
        eg[1] = ei[1];
        eg[2] = ei[2];
        eg[3] = ei[3];
    }
    if (vg + ag > 0) eg[4] = ei[4];

    // compute variables of ghost cell
    const double K[5][5] = {
        {1, 1, 0, 0, 1},
        {u - a, u, 0, 0, u + a},
        {v, v, 1, 0, v},
        {w, w, 0, 1, w},
        {H - u * a, 0.5 * V2, v, w, H + u * a},
    };
    for (long i = 0; i < 5; ++i) {
        lug[cons[i]] = 0;
        for (long j = 0; j < 5; ++j) lug[cons[i]] += K[i][j] * eg[j];
    }
    euler_local_to_global(b, lug, ug);
    euler_cons_to_prim(eqns, ug);
}
