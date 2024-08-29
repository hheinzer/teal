#pragma once

#include "equations.h"
#include "teal.h"

typedef struct Simulation Simulation;

typedef double Advance(Simulation *sim, double max_dt);

struct Simulation {
    Equations *eqns;
    String prefix;

    double time, max_time, output_time;
    long iter, max_iter, output_iter, output_count;

    long abort_variable;
    double abort_residual;

    struct {
        String name;
        Advance *method;
        double **buf;
    } advance;

    double cfl;
    long time_order;
    long n_stages;
    double tol_newton, tol_krylov;
    long dim_krylov;
    long iter_newton, iter_krylov;
};

Simulation *simulation_create(Equations *eqns, const char *prefix);

void simulation_set_max_time(Simulation *sim, double max_time);

void simulation_set_output_time(Simulation *sim, double output_time);

void simulation_set_max_iter(Simulation *sim, long max_iter);

void simulation_set_output_iter(Simulation *sim, long output_iter);

void simulation_set_abort(Simulation *sim, long variable, double residual);

void simulation_set_advance(Simulation *sim, const char *name);

void simulation_set_cfl(Simulation *sim, double cfl);

void simulation_set_time_order(Simulation *sim, long time_order);

void simulation_set_stage_count(Simulation *sim, long n_stages);

void simulation_set_newton_tolerance(Simulation *sim, double tol_newton);

void simulation_set_krylov_tolerance(Simulation *sim, double tol_krylov);

void simulation_set_krylov_dimension(Simulation *sim, long dim_krylov);

void simulation_print(const Simulation *sim);

void simulation_run(Simulation *sim);

double simulation_error(const Simulation *sim, long variable);

void simulation_free(Simulation **sim);
