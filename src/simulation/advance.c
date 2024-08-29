#include "advance.h"

#include <assert.h>
#include <float.h>

#include "gmres.h"
#include "teal/sync.h"
#include "teal/utils.h"

double explicit_euler(Simulation *sim, double max_dt)
{
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    const long n_cons = sim->eqns->n_cons;
    const long n_vars = sim->eqns->n_vars;
    const double(*dudt)[n_cons] = (void *)sim->eqns->vars.dudt;
    alias(update, sim->eqns->update.advance);
    double(*u)[n_vars] = (void *)sim->eqns->vars.u;

    const double dt0 = sim->cfl * equations_timestep(sim->eqns);
    const double dt = min(max_dt, dt0);

    equations_derivative(sim->eqns, sim->time);
    for (long i = 0; i < n_inner_cells; ++i) {
        for (long v = 0; v < n_cons; ++v) u[i][v] += dudt[i][v] * dt;
        update(sim->eqns, u[i]);
    }

    sim->time += dt;

    return dt0;
}

double lserk(Simulation *sim, double max_dt)
{
    // van Leer, 1989
    static const double alpha[3][6][7] = {
        {
            {0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.3333, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1481, 0.4000, 1.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.0833, 0.2069, 0.4265, 1.0000, 0.0000, 0.0000},
            {0.0000, 0.0533, 0.1263, 0.2375, 0.4414, 1.0000, 0.0000},
            {0.0000, 0.0370, 0.0851, 0.1521, 0.2562, 0.4512, 1.0000},
        },
        {
            {0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},  // first order
            {0.0000, 0.4242, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1918, 0.4929, 1.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1084, 0.2602, 0.5052, 1.0000, 0.0000, 0.0000},
            {0.0000, 0.0695, 0.1602, 0.2898, 0.5060, 1.0000, 0.0000},
            {0.0000, 0.0482, 0.1085, 0.1885, 0.3050, 0.5063, 1.0000},
        },
        {
            {0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},  // first order
            {0.0000, 0.6612, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.2884, 0.5010, 1.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1666, 0.3027, 0.5275, 1.0000, 0.0000, 0.0000},
            {0.0000, 0.1067, 0.1979, 0.3232, 0.5201, 1.0000, 0.0000},
            {0.0000, 0.0742, 0.1393, 0.2198, 0.3302, 0.5181, 1.0000},
        },
    };

    const long time_order = sim->time_order;
    const long n_stages = sim->n_stages;
    const double *a = alpha[time_order - 1][n_stages - 1];

    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    const long n_cons = sim->eqns->n_cons;
    const long n_vars = sim->eqns->n_vars;
    const double(*dudt)[n_cons] = (void *)sim->eqns->vars.dudt;
    alias(update, sim->eqns->update.advance);
    double(*u0)[n_cons] = (void *)sim->advance.buf[0];
    double(*u)[n_vars] = (void *)sim->eqns->vars.u;

    const double dt0 = sim->cfl * equations_timestep(sim->eqns);
    const double dt = min(max_dt, dt0);

    for (long i = 0; i < n_inner_cells; ++i)
        for (long v = 0; v < n_cons; ++v) u0[i][v] = u[i][v];

    for (long k = 0; k < n_stages; ++k) {
        equations_derivative(sim->eqns, sim->time + a[k] * dt);
        for (long i = 0; i < n_inner_cells; ++i) {
            for (long v = 0; v < n_cons; ++v) u[i][v] = u0[i][v] + a[k + 1] * dudt[i][v] * dt;
            update(sim->eqns, u[i]);
        }
    }

    sim->time += dt;

    return dt0;
}

double implicit_euler(Simulation *sim, double max_dt)
{
    // https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.KrylovJacobian.html
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    const long n_cons = sim->eqns->n_cons;
    const long nn = n_inner_cells * n_cons;
    const long n_vars = sim->eqns->n_vars;
    const long max_newton = 100 * (n_inner_cells + 1);
    const double(*dudt)[n_cons] = (void *)sim->eqns->vars.dudt;
    const double(*dx)[n_cons] = (void *)sim->advance.buf[4];
    alias(update, sim->eqns->update.advance);
    double(*u0)[n_cons] = (void *)sim->advance.buf[0];
    double(*xk)[n_cons] = (void *)sim->advance.buf[1];
    double(*fk)[n_cons] = (void *)sim->advance.buf[2];
    double(*rk)[n_cons] = (void *)sim->advance.buf[3];
    double(*u)[n_vars] = (void *)sim->eqns->vars.u;

    const double dt0 = sim->cfl * equations_timestep(sim->eqns);
    const double dt = min(max_dt, dt0);

    sim->time += dt;

    equations_derivative(sim->eqns, sim->time);
    for (long i = 0; i < n_inner_cells; ++i) {
        for (long v = 0; v < n_cons; ++v) {
            xk[i][v] = u0[i][v] = u[i][v];
            fk[i][v] = -dt * dudt[i][v];
            rk[i][v] = dudt[i][v];
        }
    }
    double norm_fk = sync_norm(*fk, nn);

    const double rdiff = 1e-6;
    double omega = rdiff * max(1.0, sync_norm(*xk, nn)) / max(1.0, norm_fk);

    const double tol_newton = sim->tol_newton * norm_fk;
    long n;
    for (n = 0; n < max_newton && norm_fk > tol_newton; ++n) {
        gmres_m(sim, norm_fk, omega, dt);

        for (long i = 0; i < n_inner_cells; ++i) {
            for (long v = 0; v < n_cons; ++v) u[i][v] = xk[i][v] += dx[i][v];
            update(sim->eqns, u[i]);
        }
        equations_derivative(sim->eqns, sim->time);
        for (long i = 0; i < n_inner_cells; ++i) {
            for (long v = 0; v < n_cons; ++v) {
                fk[i][v] = u[i][v] - u0[i][v] - dt * dudt[i][v];
                rk[i][v] = dudt[i][v];
            }
        }
        norm_fk = sync_norm(*fk, nn);

        omega = rdiff * max(1.0, sync_norm(*xk, nn)) / max(1.0, norm_fk);
    }
    sim->iter_newton += n;

    return dt0;
}
