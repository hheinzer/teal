#ifndef SIMULATION_H
#define SIMULATION_H

#include "equations.h"

typedef struct Simulation Simulation;

typedef double Advance(Simulation *sim, double max_dt);

struct Simulation {
    Equations *eqns;
    const char *prefix;

    long time_order, n_stages;
    struct {
        char name[NAMELEN];
        Advance *func;
        double *buf;
    } advance;

    double cfl;

    double time, max_time, output_time;
    long iter, max_iter, output_iter, output_count;

    long abort_variable;
    double abort_residual;
};

Simulation simulation_create(Equations *eqns, const char *prefix);

void simulation_free(Simulation *sim);

void simulation_set_time_order(Simulation *sim, long time_order, long n_stages);

void simulation_set_cfl(Simulation *sim, double cfl);

void simulation_set_max_time(Simulation *sim, double time);

void simulation_set_output_time(Simulation *sim, double time);

void simulation_set_max_iter(Simulation *sim, long iter);

void simulation_set_output_iter(Simulation *sim, long iter);

void simulation_set_abort(Simulation *sim, long variable, double residual);

void simulation_restart(Simulation *sim, const char *fname);

void simulation_print(const Simulation *sim);

void simulation_run(Simulation *sim);

void simulation_error(const Simulation *sim, long ivars, long iuser);

#endif
