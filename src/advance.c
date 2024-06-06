#include "advance.h"

#include "equations.h"

void advance_euler(Simulation *sim, const double dt)
{
    const long n_vars = sim->eqns->vars.n_fields;
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    const DERIVS(dudt, sim->eqns->vars);
    FIELDS(u, sim->eqns->vars);

    equations_time_derivative(sim->eqns, sim->time);
    for (long i = 0; i < n_inner_cells; ++i) {
        for (long v = 0; v < n_vars; ++v) u[i][v] += dudt[i][v] * dt;
        sim->eqns->update(u[i]);
    }
}

void advance_lserk(Simulation *sim, const double dt)
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
    const double *a = alpha[time_order - 1][n_stages - 1];

    const long n_vars = sim->eqns->vars.n_fields;
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    const DERIVS(dudt, sim->eqns->vars);
    FIELDS(u, sim->eqns->vars);
    double(*u0)[n_vars] = TCAST(u0, sim->buf);

    for (long i = 0; i < n_inner_cells; ++i)
        for (long v = 0; v < n_vars; ++v) u0[i][v] = u[i][v];

    for (long k = 0; k < n_stages; ++k) {
        equations_time_derivative(sim->eqns, sim->time + a[k] * dt);
        for (long i = 0; i < n_inner_cells; ++i) {
            for (long v = 0; v < n_vars; ++v) u[i][v] = u0[i][v] + a[k + 1] * dudt[i][v] * dt;
            sim->eqns->update(u[i]);
        }
    }
}
