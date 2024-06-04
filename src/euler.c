#include "euler.h"

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "equations.h"
#include "global.h"
#include "memory.h"
#include "mesh.h"
#include "riemann.h"
#include "utils.h"

static double m_gamma = 1.4;  // heat capacity ratio
static char m_flux[128];      // flux function

static void compute_conserved(double *u);
static void compute_primitive(double *u);
static TimeStep time_step;
static Flux godunov, roe, hll, hllc, hlle, lxf;
static ApplyBC symmetry, supersonic_inflow, supersonic_outflow, subsonic_inflow, subsonic_outflow,
    characteristic, custom;
static void rotate_primitives_g2l(const double *n, const double *u, double *lu);
static void rotate_fluxes_l2g(const double *n, const double *lu, double *u);
static void physical_flux(const double *u, double *flux);

Equations euler_create(const Mesh *mesh, const Fields *user)
{
    Fields vars = {.n_fields = N_VARS_EULER};
    vars.name = memory_calloc(vars.n_fields, sizeof(*vars.name));
    vars.name[D] = utils_strdup("density");
    vars.name[U] = utils_strdup("velocity-x");
    vars.name[V] = utils_strdup("velocity-y");
    vars.name[P] = utils_strdup("pressure");
    vars.name[DU] = utils_strdup("momentum-x");
    vars.name[DV] = utils_strdup("momentum-y");
    vars.name[DE] = utils_strdup("energy");

    Equations eqns = equations_create(mesh, &vars, user);
    eqns.time_step = time_step;
    eqns.update = compute_primitive;
    eqns.flux = euler_flux("roe");

    const ALIAS(bc, eqns.mesh->entity.bc.name);
    ALIAS(apply, eqns.mesh->entity.bc.apply);
    for (long e = 0; e < eqns.mesh->n_entities; ++e)
        if (bc[e]) apply[e] = euler_boundary_condition(bc[e]);

    return eqns;
}

void euler_set_gamma(const double gamma) { m_gamma = gamma; }

void euler_compute_conserved(Equations *eqns)
{
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    FIELDS(u, eqns->vars);

    for (long i = 0; i < n_inner_cells; ++i) compute_conserved(u[i]);
}

void euler_compute_primitive(Equations *eqns)
{
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    FIELDS(u, eqns->vars);

    for (long i = 0; i < n_inner_cells; ++i) compute_primitive(u[i]);
}

void euler_print(const Equations *eqns)
{
    equations_print(eqns, "Euler");
    if (eqns->mesh->rank == 0) {
        printf(" | " FMT_KEY ": %s\n", "flux function", m_flux);
        printf(" | " FMT_KEY ": %g\n", "heat capacity ratio", m_gamma);
    }
}

Flux *euler_flux(const char *flux)
{
    strncpy(m_flux, flux, sizeof(m_flux) - 1);

    if (!strcmp(flux, "godunov")) {
        return godunov;
    }
    else if (!strcmp(flux, "roe")) {
        return roe;
    }
    else if (!strcmp(flux, "hll")) {
        return hll;
    }
    else if (!strcmp(flux, "hllc")) {
        return hllc;
    }
    else if (!strcmp(flux, "hlle")) {
        return hlle;
    }
    else if (!strcmp(flux, "lxf")) {
        return lxf;
    }
    else {
        assert(!"unsupported flux function");
    }
}

ApplyBC *euler_boundary_condition(const char *bc)
{
    if (!strcmp(bc, "symmetry") || !strcmp(bc, "wall")) {
        return symmetry;
    }
    else if (!strcmp(bc, "supersonic_inflow") || !strcmp(bc, "inflow")) {
        return supersonic_inflow;
    }
    else if (!strcmp(bc, "supersonic_outflow") || !strcmp(bc, "outflow")) {
        return supersonic_outflow;
    }
    else if (!strcmp(bc, "subsonic_inflow") || !strcmp(bc, "inlet")) {
        return subsonic_inflow;
    }
    else if (!strcmp(bc, "subsonic_outflow") || !strcmp(bc, "outlet")) {
        return subsonic_outflow;
    }
    else if (!strcmp(bc, "characteristic")) {
        return characteristic;
    }
    else if (!strcmp(bc, "custom")) {
        return custom;
    }
    else {
        assert(!"unsupported boundary condition");
    }
}

static void compute_conserved(double *u)
{
    u[DU] = u[D] * u[U];
    u[DV] = u[D] * u[V];
    u[DE] = u[P] / (m_gamma - 1) + 0.5 * (u[DU] * u[U] + u[DV] * u[V]);
}

static void compute_primitive(double *u)
{
    u[D] = MAX(EPS, u[D]);
    u[U] = u[DU] / u[D];
    u[V] = u[DV] / u[D];
    u[P] = MAX(EPS, (m_gamma - 1) * (u[DE] - 0.5 * (u[DU] * u[U] + u[DV] * u[V])));
}

static double time_step(const Equations *eqns)
{
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const FIELDS(u, eqns->vars);
    const ALIAS(v, eqns->mesh->cell.volume);
    const ALIAS(p, eqns->mesh->cell.projection);

    // Blazek 2015, eqs. (6.22)
    double dt = DBL_MAX;
    for (long i = 0; i < n_inner_cells; ++i) {
        const double c = sqrt(m_gamma * u[i][P] / u[i][D]);
        dt = MIN(dt, v[i] / ((fabs(u[i][U]) + c) * p[i][0] + (fabs(u[i][V]) + c) * p[i][1]));
    }
    return dt;
}

static void godunov(const double *n, const double *ul_g, const double *ur_g, double *f_g)
{
    double ul[N_VARS_EULER], ur[N_VARS_EULER];
    rotate_primitives_g2l(n, ul_g, ul);
    rotate_primitives_g2l(n, ur_g, ur);

    // Toro 2009, sec. 6.3
    double u[N_VARS_EULER];
    riemann(m_gamma, ul[D], ul[U], ul[P], ur[D], ur[U], ur[P], 0, &u[D], &u[U], &u[P]);
    u[V] = (u[U] >= 0 ? ul[V] : ur[V]);
    compute_conserved(u);

    double f[N_VARS_EULER];
    physical_flux(u, f);
    rotate_fluxes_l2g(n, f, f_g);
}

static void roe(const double *n, const double *ul_g, const double *ur_g, double *f_g)
{
    double ul[N_VARS_EULER], ur[N_VARS_EULER];
    rotate_primitives_g2l(n, ul_g, ul);
    rotate_primitives_g2l(n, ur_g, ur);
    compute_conserved(ul);
    compute_conserved(ur);

    // Toro 2009, sec. 11.2.2
    const double hl = (ul[DE] + ul[P]) / ul[D];
    const double hr = (ur[DE] + ur[P]) / ur[D];
    const double ds = sqrt(ul[D]) + sqrt(ur[D]);
    const double us = (sqrt(ul[D]) * ul[U] + sqrt(ur[D]) * ur[U]) / ds;
    const double vs = (sqrt(ul[D]) * ul[V] + sqrt(ur[D]) * ur[V]) / ds;
    const double hs = (sqrt(ul[D]) * hl + sqrt(ur[D]) * hr) / ds;
    const double q2 = us * us + vs * vs;
    const double cs = sqrt((m_gamma - 1) * (hs - 0.5 * q2));
    double e[4] = {us - cs, us, us, us + cs};

    // second entropy fix of Harten and Hyman, Pelanti 2018
    const double cl = sqrt(m_gamma * ul[P] / ul[D]);
    const double cr = sqrt(m_gamma * ur[P] / ur[D]);
    const double el[4] = {ul[U] - cl, ul[U], ul[U], ul[U] + cl};
    const double er[4] = {ur[U] - cr, ur[U], ur[U], ur[U] + cr};
    for (long i = 0; i < 4; ++i) {
        const double delta = MAX(EPS, MAX(e[i] - el[i], er[i] - e[i]));
        if (fabs(e[i]) < delta) {
            e[i] = 0.5 * (e[i] * e[i] / delta + delta);
        }
        else {
            e[i] = fabs(e[i]);
        }
    }

    const double v[4][4] = {
        {1, 1, 0, 1},
        {us - cs, us, 0, us + cs},
        {vs, vs, 1, vs},
        {hs - us * cs, 0.5 * q2, vs, hs + us * cs},
    };
    const double dd = ur[D] - ul[D];
    const double ddu = ur[DU] - ul[DU];
    const double ddv = ur[DV] - ul[DV];
    const double dde = ur[DE] - ul[DE] - (ddv - vs * dd) * vs;
    double a[4];
    a[2] = ddv - vs * dd;
    a[1] = (m_gamma - 1) / (cs * cs) * (dd * (hs - us * us) + us * ddu - dde);
    a[0] = 0.5 / cs * (dd * (us + cs) - ddu - cs * a[1]);
    a[3] = dd - (a[0] + a[1]);
    double df[4] = {};
    for (long i = 0; i < 4; ++i)
        for (long j = 0; j < 4; ++j) df[i] += a[j] * fabs(e[j]) * v[i][j];

    double fl[N_VARS_EULER], fr[N_VARS_EULER], f[N_VARS_EULER];
    physical_flux(ul, fl);
    physical_flux(ur, fr);
    f[D] = 0.5 * (fl[D] + fr[D] - df[0]);
    f[DU] = 0.5 * (fl[DU] + fr[DU] - df[1]);
    f[DV] = 0.5 * (fl[DV] + fr[DV] - df[2]);
    f[DE] = 0.5 * (fl[DE] + fr[DE] - df[3]);
    rotate_fluxes_l2g(n, f, f_g);
}

static void hll(const double *n, const double *ul_g, const double *ur_g, double *f_g)
{
    double ul[N_VARS_EULER], ur[N_VARS_EULER];
    rotate_primitives_g2l(n, ul_g, ul);
    rotate_primitives_g2l(n, ur_g, ur);
    compute_conserved(ul);
    compute_conserved(ur);

    // Toro 2009, sec. 10.3
    const double hl = (ul[DE] + ul[P]) / ul[D];
    const double hr = (ur[DE] + ur[P]) / ur[D];
    const double ds = sqrt(ul[D]) + sqrt(ur[D]);
    const double us = (sqrt(ul[D]) * ul[U] + sqrt(ur[D]) * ur[U]) / ds;
    const double vs = (sqrt(ul[D]) * ul[V] + sqrt(ur[D]) * ur[V]) / ds;
    const double hs = (sqrt(ul[D]) * hl + sqrt(ur[D]) * hr) / ds;
    const double q2 = us * us + vs * vs;
    const double cs = sqrt((m_gamma - 1) * (hs - 0.5 * q2));
    const double cl = sqrt(m_gamma * ul[P] / ul[D]);
    const double cr = sqrt(m_gamma * ur[P] / ur[D]);
    const double sl = MIN(ul[U] - cl, us - cs);
    const double sr = MAX(ur[U] + cr, us + cs);

    double fl[N_VARS_EULER], fr[N_VARS_EULER], f[N_VARS_EULER];
    if (0 <= sl) {
        physical_flux(ul, f);
    }
    else if (0 >= sr) {
        physical_flux(ur, f);
    }
    else {
        const double du[4] = {
            ur[D] - ul[D],
            ur[DU] - ul[DU],
            ur[DV] - ul[DV],
            ur[DE] - ul[DE],
        };
        physical_flux(ul, fl);
        physical_flux(ur, fr);
        f[D] = (sr * fl[D] - sl * fr[D] + sl * sr * du[0]) / (sr - sl);
        f[DU] = (sr * fl[DU] - sl * fr[DU] + sl * sr * du[1]) / (sr - sl);
        f[DV] = (sr * fl[DV] - sl * fr[DV] + sl * sr * du[2]) / (sr - sl);
        f[DE] = (sr * fl[DE] - sl * fr[DE] + sl * sr * du[3]) / (sr - sl);
    }
    rotate_fluxes_l2g(n, f, f_g);
}

static void hllc(const double *n, const double *ul_g, const double *ur_g, double *f_g)
{
    double ul[N_VARS_EULER], ur[N_VARS_EULER];
    rotate_primitives_g2l(n, ul_g, ul);
    rotate_primitives_g2l(n, ur_g, ur);
    compute_conserved(ul);
    compute_conserved(ur);

    // Toro 2009, sec. 10.4.2
    const double hl = (ul[DE] + ul[P]) / ul[D];
    const double hr = (ur[DE] + ur[P]) / ur[D];
    const double ds = sqrt(ul[D]) + sqrt(ur[D]);
    const double us = (sqrt(ul[D]) * ul[U] + sqrt(ur[D]) * ur[U]) / ds;
    const double vs = (sqrt(ul[D]) * ul[V] + sqrt(ur[D]) * ur[V]) / ds;
    const double hs = (sqrt(ul[D]) * hl + sqrt(ur[D]) * hr) / ds;
    const double q2 = us * us + vs * vs;
    const double cs = sqrt((m_gamma - 1) * (hs - 0.5 * q2));
    const double cl = sqrt(m_gamma * ul[P] / ul[D]);
    const double cr = sqrt(m_gamma * ur[P] / ur[D]);
    const double sl = MIN(ul[U] - cl, us - cs);
    const double sr = MAX(ur[U] + cr, us + cs);

    double f[N_VARS_EULER];
    if (0 <= sl) {
        physical_flux(ul, f);
    }
    else if (0 >= sr) {
        physical_flux(ur, f);
    }
    else {
        const double ss = (ur[P] - ul[P] + ul[DU] * (sl - ul[U]) - ur[DU] * (sr - ur[U])) /
                          (ul[D] * (sl - ul[U]) - ur[D] * (sr - ur[U]));
        if (sl <= 0 && 0 <= ss) {
            const double fac = ul[D] * (sl - ul[U]) / (sl - ss);
            const double u_s[4] = {
                1,
                ss,
                ul[V],
                ul[DE] / ul[D] + (ss - ul[U]) * (ss + ul[P] / (ul[D] * (sl - ul[U]))),
            };
            physical_flux(ul, f);
            f[D] += sl * (fac * u_s[0] - ul[D]);
            f[DU] += sl * (fac * u_s[1] - ul[DU]);
            f[DV] += sl * (fac * u_s[2] - ul[DV]);
            f[DE] += sl * (fac * u_s[3] - ul[DE]);
        }
        else {
            const double fac = ur[D] * (sr - ur[U]) / (sr - ss);
            const double u_s[4] = {
                1,
                ss,
                ur[V],
                ur[DE] / ur[D] + (ss - ur[U]) * (ss + ur[P] / (ur[D] * (sr - ur[U]))),
            };
            physical_flux(ur, f);
            f[D] += sr * (fac * u_s[0] - ur[D]);
            f[DU] += sr * (fac * u_s[1] - ur[DU]);
            f[DV] += sr * (fac * u_s[2] - ur[DV]);
            f[DE] += sr * (fac * u_s[3] - ur[DE]);
        }
    }
    rotate_fluxes_l2g(n, f, f_g);
}

static void hlle(const double *n, const double *ul_g, const double *ur_g, double *f_g)
{
    double ul[N_VARS_EULER], ur[N_VARS_EULER];
    rotate_primitives_g2l(n, ul_g, ul);
    rotate_primitives_g2l(n, ur_g, ur);
    compute_conserved(ul);
    compute_conserved(ur);

    // Einfeldt 1988
    const double ds = sqrt(ul[D]) + sqrt(ur[D]);
    const double us = (sqrt(ul[D]) * ul[U] + sqrt(ur[D]) * ur[U]) / ds;
    const double cl = sqrt(m_gamma * ul[P] / ul[D]);
    const double cr = sqrt(m_gamma * ur[P] / ur[D]);
    const double d =
        sqrt((sqrt(ul[D]) * cl * cl + sqrt(ur[D]) * cr * cr) / ds +
             0.5 * (sqrt(ul[D]) * sqrt(ur[D])) / (ds * ds) * (ur[U] - ul[U]) * (ur[U] - ul[U]));
    const double sl = MIN(ul[U] - cl, us - d);
    const double sr = MAX(ur[U] + cr, us + d);

    double fl[N_VARS_EULER], fr[N_VARS_EULER], f[N_VARS_EULER];
    if (0 <= sl) {
        physical_flux(ul, f);
    }
    else if (0 >= sr) {
        physical_flux(ur, f);
    }
    else {
        const double du[4] = {
            ur[D] - ul[D],
            ur[DU] - ul[DU],
            ur[DV] - ul[DV],
            ur[DE] - ul[DE],
        };
        physical_flux(ul, fl);
        physical_flux(ur, fr);
        f[D] = (sr * fl[D] - sl * fr[D] + sl * sr * du[0]) / (sr - sl);
        f[DU] = (sr * fl[DU] - sl * fr[DU] + sl * sr * du[1]) / (sr - sl);
        f[DV] = (sr * fl[DV] - sl * fr[DV] + sl * sr * du[2]) / (sr - sl);
        f[DE] = (sr * fl[DE] - sl * fr[DE] + sl * sr * du[3]) / (sr - sl);
    }
    rotate_fluxes_l2g(n, f, f_g);
}

static void lxf(const double *n, const double *ul_g, const double *ur_g, double *f_g)
{
    double ul[N_VARS_EULER], ur[N_VARS_EULER];
    rotate_primitives_g2l(n, ul_g, ul);
    rotate_primitives_g2l(n, ur_g, ur);
    compute_conserved(ul);
    compute_conserved(ur);

    // Lax-Friedrichs flux
    const double cl = sqrt(m_gamma * ul[P] / ul[D]);
    const double cr = sqrt(m_gamma * ur[P] / ur[D]);
    const double s = MAX(fabs(ul[U]) + cl, fabs(ur[U]) + cr);

    double fl[N_VARS_EULER], fr[N_VARS_EULER], f[N_VARS_EULER];
    physical_flux(ul, fl);
    physical_flux(ur, fr);
    f[D] = 0.5 * (fl[D] + fr[D] - s * (ur[D] - ul[D]));
    f[DU] = 0.5 * (fl[DU] + fr[DU] - s * (ur[DU] - ul[DU]));
    f[DV] = 0.5 * (fl[DV] + fr[DV] - s * (ur[DV] - ul[DV]));
    f[DE] = 0.5 * (fl[DE] + fr[DE] - s * (ur[DE] - ul[DE]));
    rotate_fluxes_l2g(n, f, f_g);
}

static void symmetry(const BCContext, const double *n, const double *ui, double *ug)
{
    // Blazek 2015, sec. 8.6
    const double vn = n[0] * ui[U] + n[1] * ui[V];
    ug[D] = ui[D];
    ug[U] = ui[U] - 2 * vn * n[0];
    ug[V] = ui[V] - 2 * vn * n[1];
    ug[P] = ui[P];
}

static void supersonic_inflow(const BCContext context, const double *, const double *, double *ug)
{
    // Blazek 2015, sec. 8.3.1
    const double us[] = {
        [D] = context.state[0],
        [U] = context.state[1],
        [V] = context.state[2],
        [P] = context.state[3],
    };
    ug[D] = us[D];
    ug[U] = us[U];
    ug[V] = us[V];
    ug[P] = us[P];
}

static void supersonic_outflow(const BCContext, const double *, const double *ui, double *ug)
{
    // Blazek 2015, sec. 8.3.1
    ug[D] = ui[D];
    ug[U] = ui[U];
    ug[V] = ui[V];
    ug[P] = ui[P];
}

static void subsonic_inflow(const BCContext context, const double *n, const double *ui, double *ug)
{
    // Blazek 2015, sec. 8.3.1
    const double us[] = {
        [D] = context.state[0],
        [U] = context.state[1],
        [V] = context.state[2],
        [P] = context.state[3],
    };
    const double rho0 = ui[D];
    const double c0 = sqrt(m_gamma * ui[P] / ui[D]);
    const double dvn = n[0] * (us[U] - ui[U]) + n[1] * (us[V] - ui[V]);
    ug[P] = 0.5 * (us[P] + ui[P] - rho0 * c0 * dvn);
    ug[D] = us[D] + (ug[P] - us[P]) / (c0 * c0);
    ug[U] = us[U] - n[0] * (us[P] - ug[P]) / (rho0 * c0);
    ug[V] = us[V] - n[1] * (us[P] - ug[P]) / (rho0 * c0);
}

static void subsonic_outflow(const BCContext context, const double *n, const double *ui, double *ug)
{
    // Blazek 2015, sec. 8.3.1
    const double us[] = {
        [P] = context.state[0],
    };
    const double rho0 = ui[D];
    const double c0 = sqrt(m_gamma * ui[P] / ui[D]);
    ug[P] = us[P];
    ug[D] = ui[D] + (ug[P] - ui[P]) / (c0 * c0);
    ug[U] = ui[U] - n[0] * (ui[P] - ug[P]) / (rho0 * c0);
    ug[V] = ui[V] - n[1] * (ui[P] - ug[P]) / (rho0 * c0);
}

static void characteristic(const BCContext context, const double *n, const double *ui_g,
                           double *ug_g)
{
    // compute initial ghost cell state
    ug_g[D] = context.state[0];
    ug_g[U] = context.state[1];
    ug_g[V] = context.state[2];
    ug_g[P] = context.state[3];

    // compute eigenvalues of ghost cell
    const double cg = sqrt(m_gamma * ug_g[P] / ug_g[D]);
    const double vg = n[0] * ug_g[U] + n[1] * ug_g[V];

    // rotate variables
    double ui[N_VARS_EULER], ug[N_VARS_EULER];
    rotate_primitives_g2l(n, ui_g, ui);
    rotate_primitives_g2l(n, ug_g, ug);
    compute_conserved(ui);
    compute_conserved(ug);

    // compute characteristic variables
    const double c = sqrt(m_gamma * ui_g[P] / ui_g[D]);
    const double u = ui_g[U];
    const double v = ui_g[V];
    const double q2 = u * u + v * v;
    const double h = c * c / (m_gamma - 1) + 0.5 * q2;
    const double a = 2 * h - q2;
    const double Ki[4][4] = {
        {u / (2 * c) + q2 / (2 * a), -1 / (2 * c) - u / a, -v / a, 1 / a},
        {(a - q2) / a, 2 * u / a, 2 * v / a, -2 / a},
        {-v, 0, 1, 0},
        {-u / (2 * c) + q2 / (2 * a), 1 / (2 * c) - u / a, -v / a, 1 / a},
    };
    double ei[4], eg[4];
    for (long i = 0; i < 4; ++i)
        ei[i] = Ki[i][0] * ui[D] + Ki[i][1] * ui[DU] + Ki[i][2] * ui[DV] + Ki[i][3] * ui[DE];
    for (long i = 0; i < 4; ++i)
        eg[i] = Ki[i][0] * ug[D] + Ki[i][1] * ug[DU] + Ki[i][2] * ug[DV] + Ki[i][3] * ug[DE];

    // compute characteristic variables of ghost cell
    if (vg - cg > 0) eg[0] = ei[0];
    if (vg > 0) {
        eg[1] = ei[1];
        eg[2] = ei[2];
    }
    if (vg + cg > 0) eg[3] = ei[3];

    // compute variables of ghost cell
    const double K[4][4] = {
        {1, 1, 0, 1},
        {-c + u, u, 0, c + u},
        {v, v, 1, v},
        {-c * u + h, q2 / 2, v, c * u + h},
    };
    ug[D] = K[0][0] * eg[0] + K[0][1] * eg[1] + K[0][2] * eg[2] + K[0][3] * eg[3];
    ug[DU] = K[1][0] * eg[0] + K[1][1] * eg[1] + K[1][2] * eg[2] + K[1][3] * eg[3];
    ug[DV] = K[2][0] * eg[0] + K[2][1] * eg[1] + K[2][2] * eg[2] + K[2][3] * eg[3];
    ug[DE] = K[3][0] * eg[0] + K[3][1] * eg[1] + K[3][2] * eg[2] + K[3][3] * eg[3];
    compute_primitive(ug);

    // rotate ghost cell variables
    ug_g[D] = ug[D];
    ug_g[U] = n[0] * ug[U] - n[1] * ug[V];
    ug_g[V] = n[0] * ug[V] + n[1] * ug[U];
    ug_g[P] = ug[P];
}

static void custom(const BCContext context, const double *n, const double *ui, double *ug)
{
    const double *x = n;
    const double time = context.state[0];
    context.custom(x, time, ui, ug);
}

static void rotate_primitives_g2l(const double *n, const double *u_g, double *u)
{
    u[D] = u_g[D];
    u[U] = n[0] * u_g[U] + n[1] * u_g[V];
    u[V] = n[0] * u_g[V] - n[1] * u_g[U];
    u[P] = u_g[P];
}

static void rotate_fluxes_l2g(const double *n, const double *f, double *f_g)
{
    f_g[D] = f[D];
    f_g[DU] = n[0] * f[DU] - n[1] * f[DV];
    f_g[DV] = n[0] * f[DV] + n[1] * f[DU];
    f_g[DE] = f[DE];
}

static void physical_flux(const double *u, double *flux)
{
    flux[D] = u[DU];
    flux[DU] = u[DU] * u[U] + u[P];
    flux[DV] = u[DU] * u[V];
    flux[DE] = (u[DE] + u[P]) * u[U];
}
