#ifndef SIMULATION_H
#define SIMULATION_H

#include <stdio.h>

#include "equations.h"

typedef struct Simulation Simulation;
typedef void Advance(Simulation *sim, const double dt);
struct Simulation {
    long time_order, n_stages;
    double time, max_time, output_time;
    long iter, max_iter, output_iter, output_count;
    double cfl;

    long abort_variable;
    double abort_residual;
    FILE *residual_file;

    Advance *advance;
    double *buf;

    char *prefix;
    Equations *eqns;
};

Simulation simulation_create(const char *prefix, Equations *eqns);

void simulation_free(Simulation *sim);

void simulation_set_time_order(Simulation *sim, const long time_order, const long n_stages);

void simulation_print(const Simulation *sim);

void simulation_run(Simulation *sim);

void simulation_error(const Simulation *sim, const long i_vars, const long i_user);

#endif
