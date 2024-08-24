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

    simulation_set_explicit(&sim, eqns->space_order, eqns->space_order + 1);

    simulation_set_cfl(&sim, 0.99);
    simulation_set_max_time(&sim, DBL_MAX);
    simulation_set_output_time(&sim, DBL_MAX);
    simulation_set_max_iter(&sim, LONG_MAX);
    simulation_set_output_iter(&sim, LONG_MAX);
    simulation_set_abort(&sim, -1, 0);

    if (teal.restart) simulation_restart(&sim, teal.restart);

    return sim;
}

void simulation_set_explicit(Simulation *sim, long time_order, long n_stages)
{
    if (time_order < 1 || 3 < time_order) error("unsupported time order '%ld'", time_order);
    if (n_stages < 1 || 6 < n_stages) error("unsupported number of stages '%ld'", n_stages);
    sim->time_order = time_order;
    sim->n_stages = n_stages;

    sim->n_newton = 0;
    sim->n_krylov = 0;
    sim->tol_newton = 0;
    sim->tol_krylov = 0;

    const long n_cons = sim->eqns->n_cons;
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    memory_free(&sim->advance.u0);
    memory_free(&sim->advance.xk);
    memory_free(&sim->advance.f0);
    memory_free(&sim->advance.fk);
    memory_free(&sim->advance.rk);
    memory_free(&sim->advance.dx);
    memory_free(&sim->advance.V);
    memory_free(&sim->advance.g);
    memory_free(&sim->advance.w);
    memory_free(&sim->advance.H);
    memory_free(&sim->advance.s);
    memory_free(&sim->advance.c);
    memory_free(&sim->advance.y);
    if (time_order == 1 && n_stages == 1) {
        strcpy(sim->advance.name, "explicit euler");
        sim->advance.func = explicit_euler;
    }
    else {
        strcpy(sim->advance.name, "lserk");
        sim->advance.func = lserk;
        sim->advance.u0 = memory_calloc(n_inner_cells * n_cons, sizeof(*sim->advance.u0));
    }
}

void simulation_set_implicit(Simulation *sim, long n_newton, long n_krylov, double tol_newton,
                             double tol_krylov)
{
    if (n_newton <= 0) error("unsupported number of Newton iterations '%ld'", n_newton);
    if (n_krylov <= 0) error("unsupported number of Krylov iterations '%ld'", n_krylov);
    if (tol_newton < 0) error("unsupported Newton tolerance '%g'", tol_newton);
    if (tol_krylov < 0) error("unsupported Krylov tolerance '%g'", tol_krylov);
    sim->n_newton = n_newton;
    sim->n_krylov = n_krylov;
    sim->tol_newton = tol_newton;
    sim->tol_krylov = tol_krylov;

    sim->time_order = 0;
    sim->n_stages = 0;

    const long n_cons = sim->eqns->n_cons;
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    memory_free(&sim->advance.u0);
    memory_free(&sim->advance.xk);
    memory_free(&sim->advance.f0);
    memory_free(&sim->advance.fk);
    memory_free(&sim->advance.rk);
    memory_free(&sim->advance.dx);
    memory_free(&sim->advance.V);
    memory_free(&sim->advance.g);
    memory_free(&sim->advance.w);
    memory_free(&sim->advance.H);
    memory_free(&sim->advance.s);
    memory_free(&sim->advance.c);
    memory_free(&sim->advance.y);
    strcpy(sim->advance.name, "implicit euler");
    sim->advance.func = implicit_euler;
    sim->advance.u0 = memory_calloc(n_inner_cells * n_cons, sizeof(*sim->advance.u0));
    sim->advance.xk = memory_calloc(n_inner_cells * n_cons, sizeof(sim->advance.xk));
    sim->advance.f0 = memory_calloc(n_inner_cells * n_cons, sizeof(sim->advance.f0));
    sim->advance.fk = memory_calloc(n_inner_cells * n_cons, sizeof(sim->advance.fk));
    sim->advance.rk = memory_calloc(n_inner_cells * n_cons, sizeof(sim->advance.rk));
    sim->advance.dx = memory_calloc(n_inner_cells * n_cons, sizeof(sim->advance.dx));
    sim->advance.V = memory_calloc((n_krylov + 1) * n_inner_cells * n_cons, sizeof(sim->advance.V));
    sim->advance.g = memory_calloc(n_krylov + 1, sizeof(sim->advance.g));
    sim->advance.w = memory_calloc(n_inner_cells * n_cons, sizeof(sim->advance.w));
    sim->advance.H = memory_calloc((n_krylov + 1) * n_krylov, sizeof(sim->advance.H));
    sim->advance.s = memory_calloc(n_krylov, sizeof(sim->advance.s));
    sim->advance.c = memory_calloc(n_krylov, sizeof(sim->advance.c));
    sim->advance.y = memory_calloc(n_krylov, sizeof(sim->advance.y));
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
