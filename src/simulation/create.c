#include <float.h>
#include <limits.h>
#include <string.h>

#include "advance.h"
#include "core/memory.h"
#include "core/utils.h"
#include "simulation.h"
#include "teal.h"

Simulation simulation_create(Equations *eqns, const char *prefix)
{
    Simulation sim = {0};

    sim.eqns = eqns;
    sim.prefix = prefix;

    simulation_set_time_order(&sim, eqns->space_order, eqns->space_order + 1);
    simulation_set_cfl(&sim, 0.99);
    simulation_set_max_time(&sim, DBL_MAX);
    simulation_set_output_time(&sim, DBL_MAX);
    simulation_set_max_iter(&sim, LONG_MAX);
    simulation_set_output_iter(&sim, LONG_MAX);
    simulation_set_abort(&sim, -1, 0);

    if (teal.restart) simulation_restart(&sim, teal.restart);

    return sim;
}

void simulation_set_time_order(Simulation *sim, long time_order, long n_stages)
{
    if (time_order < 1 || 3 < time_order) error("unsupported time order '%ld'", time_order);
    if (n_stages < 1 || 6 < n_stages) error("unsupported number of stages '%ld'", n_stages);
    sim->time_order = time_order;
    sim->n_stages = n_stages;

    const long n_vars = sim->eqns->n_vars;
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    memory_free(&sim->advance.buf);
    if (time_order == 1 && n_stages == 1) {
        strcpy(sim->advance.name, "euler");
        sim->advance.func = euler;
    }
    else {
        strcpy(sim->advance.name, "lserk");
        sim->advance.func = lserk;
        sim->advance.buf = memory_calloc(n_inner_cells * n_vars, sizeof(*sim->advance.buf));
    }
}

void simulation_set_cfl(Simulation *sim, double cfl)
{
    sim->cfl = cfl;
}

void simulation_set_max_time(Simulation *sim, double max_time)
{
    sim->max_time = sim->time + max_time;
}

void simulation_set_output_time(Simulation *sim, double output_time)
{
    sim->output_time = output_time;
}

void simulation_set_max_iter(Simulation *sim, long max_iter)
{
    sim->max_iter = sim->iter + max_iter;
}

void simulation_set_output_iter(Simulation *sim, long output_iter)
{
    sim->output_iter = output_iter;
}

void simulation_set_abort(Simulation *sim, long variable, double residual)
{
    sim->abort_variable = variable;
    sim->abort_residual = residual;
}
