#include "boundary.h"

#include <math.h>

#include "core/utils.h"
#include "euler.h"
#include "rotate.h"
#include "update.h"

void symmetry(const Equations *, const double *n, const double *, const double *ui, double *ug)
{
    // Blazek 2015, sec. 8.6
    const double vn = ui[U] * n[X] + ui[V] * n[Y] + ui[W] * n[Z];
    ug[D] = ui[D];
    ug[U] = ui[U] - 2 * vn * n[X];
    ug[V] = ui[V] - 2 * vn * n[Y];
    ug[W] = ui[W] - 2 * vn * n[Z];
    ug[P] = ui[P];
}

void supersonic_inflow(const Equations *, const double *, const double *state, const double *,
                       double *ug)
{
    // Blazek 2015, sec. 8.3.1
    ug[D] = state[D];
    ug[U] = state[U];
    ug[V] = state[V];
    ug[W] = state[W];
    ug[P] = state[P];
}

void supersonic_outflow(const Equations *, const double *, const double *, const double *ui,
                        double *ug)
{
    // Blazek 2015, sec. 8.3.1
    ug[D] = ui[D];
    ug[U] = ui[U];
    ug[V] = ui[V];
    ug[W] = ui[W];
    ug[P] = ui[P];
}

void subsonic_inflow(const Equations *eqns, const double *n, const double *state, const double *ui,
                     double *ug)
{
    // Blazek 2015, sec. 8.3.1
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

void subsonic_outflow(const Equations *eqns, const double *n, const double *state, const double *ui,
                      double *ug)
{
    // Blazek 2015, sec. 8.3.1
    const double gamma = eqns->scalar.value[GAMMA];
    const double d0 = ui[D];
    const double c0 = sqrt(gamma * ui[P] / ui[D]);
    ug[P] = state[P];
    ug[D] = ui[D] + (ug[P] - ui[P]) / sq(c0);
    ug[U] = ui[U] + n[X] * (ui[P] - ug[P]) / (d0 * c0);
    ug[V] = ui[V] + n[Y] * (ui[P] - ug[P]) / (d0 * c0);
    ug[W] = ui[W] + n[Z] * (ui[P] - ug[P]) / (d0 * c0);
}

void farfield(const Equations *eqns, const double *n, const double *state, const double *ui,
              double *ug)
{
    // compute initial ghost cell state
    ug[D] = state[D];
    ug[U] = state[U];
    ug[V] = state[V];
    ug[W] = state[W];
    ug[P] = state[P];

    // rotate variables
    double l_ui[N_VARS], l_ug[N_VARS];
    global_to_local(n, ui, l_ui);
    global_to_local(n, ug, l_ug);
    prim_to_cons(eqns, l_ui);
    prim_to_cons(eqns, l_ug);

    // compute characteristic variables
    const double gamma = eqns->scalar.value[GAMMA];
    const double g = gamma - 1;
    const double a = sqrt(gamma * l_ui[P] / l_ui[D]);
    const double u = l_ui[U];
    const double v = l_ui[V];
    const double w = l_ui[W];
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
            ei[i] += invK[i][j] * l_ui[cons[j]];
            eg[i] += invK[i][j] * l_ug[cons[j]];
        }
        ei[i] *= fac;
        eg[i] *= fac;
    }

    // compute characteristic variables of ghost cell
    const double ag = sqrt(gamma * l_ug[P] / l_ug[D]);
    const double vg = l_ug[U];
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
        l_ug[cons[i]] = 0;
        for (long j = 0; j < 5; ++j) l_ug[cons[i]] += K[i][j] * eg[j];
    }
    local_to_global(n, l_ug, ug);
    cons_to_prim(eqns, ug);
}
