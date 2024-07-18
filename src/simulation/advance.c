#include "advance.h"

#include "core/memory.h"
#include "core/utils.h"
#include "equations.h"

double euler(Simulation *sim, const double max_dt)
{
    const long n_vars = sim->eqns->n_vars;
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    const double(*dudt)[n_vars] = (void *)sim->eqns->vars.dudt;
    ALIAS(update, sim->eqns->advance);
    double(*u)[n_vars] = (void *)sim->eqns->vars.u;

    const double dt0 = sim->cfl * equations_timestep(sim->eqns);
    const double dt = min(max_dt, dt0);

    equations_derivative(sim->eqns, sim->time);
    for (long i = 0; i < n_inner_cells; ++i) {
        for (long v = 0; v < n_vars; ++v) u[i][v] += dudt[i][v] * dt;
        update(sim->eqns, u[i]);
    }

    sim->time += dt;

    return dt0;
}

double lserk(Simulation *sim, const double max_dt)
{
    // van Leer, 1989
    static const double ALPHA[3][6][7] = {
        {
            {0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.3333, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1481, 0.4000, 1.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.0833, 0.2069, 0.4265, 1.0000, 0.0000, 0.0000},
            {0.0000, 0.0533, 0.1263, 0.2375, 0.4414, 1.0000, 0.0000},
            {0.0000, 0.0370, 0.0851, 0.1521, 0.2562, 0.4512, 1.0000},
        },
        {
            {0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},  // not second order
            {0.0000, 0.4242, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1918, 0.4929, 1.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1084, 0.2602, 0.5052, 1.0000, 0.0000, 0.0000},
            {0.0000, 0.0695, 0.1602, 0.2898, 0.5060, 1.0000, 0.0000},
            {0.0000, 0.0482, 0.1085, 0.1885, 0.3050, 0.5063, 1.0000},
        },
        {
            {0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},  // not third order
            {0.0000, 0.6612, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000},  // maybe third order
            {0.0000, 0.2884, 0.5010, 1.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1666, 0.3027, 0.5275, 1.0000, 0.0000, 0.0000},
            {0.0000, 0.1067, 0.1979, 0.3232, 0.5201, 1.0000, 0.0000},
            {0.0000, 0.0742, 0.1393, 0.2198, 0.3302, 0.5181, 1.0000},
        },
    };
    const long time_order = sim->time_order;
    const long n_stages = sim->n_stages;
    const double *alpha = ALPHA[time_order - 1][n_stages - 1];

    const long n_vars = sim->eqns->n_vars;
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    const double(*dudt)[n_vars] = (void *)sim->eqns->vars.dudt;
    ALIAS(update, sim->eqns->advance);
    double(*u0)[n_vars] = (void *)sim->advance.buf;
    double(*u)[n_vars] = (void *)sim->eqns->vars.u;

    const double dt0 = sim->cfl * equations_timestep(sim->eqns);
    const double dt = min(max_dt, dt0);

    memory_copy(u0, u, n_inner_cells, sizeof(*u));
    for (long k = 0; k < n_stages; ++k) {
        equations_derivative(sim->eqns, sim->time + alpha[k] * dt);
        for (long i = 0; i < n_inner_cells; ++i) {
            for (long v = 0; v < n_vars; ++v) u[i][v] = u0[i][v] + alpha[k + 1] * dudt[i][v] * dt;
            update(sim->eqns, u[i]);
        }
    }

    sim->time += dt;

    return dt0;
}
