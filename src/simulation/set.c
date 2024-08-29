#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "advance.h"
#include "free.h"
#include "simulation.h"
#include "teal/memory.h"

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
    sim->max_iter = max_iter;
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

void simulation_set_advance(Simulation *sim, const char *name)
{
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    const long n_cons = sim->eqns->n_cons;
    const long nn = n_inner_cells * n_cons;
    const long space_order = sim->eqns->space_order;

    sim->time_order = 0;
    sim->n_stages = 0;
    sim->tol_newton = 0;
    sim->dim_krylov = 0;

    advance_free(sim);
    strcpy(sim->advance.name, name);
    if (!strcmp(name, "explicit euler"))
        sim->advance.method = explicit_euler;
    else if (!strcmp(name, "lserk")) {
        sim->time_order = space_order;
        sim->n_stages = space_order + 1;
        sim->advance.method = lserk;
        sim->advance.buf = memory_calloc(1, sizeof(*sim->advance.buf));
        sim->advance.buf[0] = memory_calloc(nn, sizeof(**sim->advance.buf));
    }
    else if (!strcmp(name, "implicit euler")) {
        sim->tol_newton = 0.1;
        sim->tol_krylov = 0.3;
        sim->dim_krylov = 32;
        sim->advance.method = implicit_euler;
        sim->advance.buf = memory_calloc(12, sizeof(*sim->advance.buf));

        const long mm = sim->dim_krylov + 1;
        for (long i = 0; i < 6; ++i)
            sim->advance.buf[i] = memory_calloc(nn, sizeof(**sim->advance.buf));
        sim->advance.buf[6] = memory_calloc(nn * mm, sizeof(**sim->advance.buf));
        sim->advance.buf[7] = memory_calloc(mm * mm, sizeof(**sim->advance.buf));
        for (long i = 8; i < 12; ++i)
            sim->advance.buf[i] = memory_calloc(mm, sizeof(**sim->advance.buf));
    }
    else
        abort();
}

void simulation_set_cfl(Simulation *sim, double cfl)
{
    sim->cfl = cfl;
}

void simulation_set_time_order(Simulation *sim, long time_order)
{
    assert(1 <= time_order && time_order <= 3);
    sim->time_order = time_order;
}

void simulation_set_stage_count(Simulation *sim, long n_stages)
{
    assert(1 <= n_stages && n_stages <= 6);
    sim->n_stages = n_stages;
}

void simulation_set_newton_tolerance(Simulation *sim, double tol_newton)
{
    assert(tol_newton > 0);
    sim->tol_newton = tol_newton;
}

void simulation_set_krylov_tolerance(Simulation *sim, double tol_krylov)
{
    assert(tol_krylov > 0);
    sim->tol_krylov = tol_krylov;
}

void simulation_set_krylov_dimension(Simulation *sim, long dim_krylov)
{
    assert(dim_krylov > 0);
    sim->dim_krylov = dim_krylov;
}
