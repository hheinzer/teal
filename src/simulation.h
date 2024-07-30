#pragma once

#include "equations.h"

typedef struct Simulation Simulation;

/* Advance 'sim' in time, using a maximum time step of 'max_dt'. */
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

/* Create simulation of 'eqns' with the output 'prefix'. */
Simulation simulation_create(Equations *eqns, const char *prefix);

/* Deallocate memory of 'sim' and set all fields to 0. */
void simulation_free(Simulation *sim);

/* Set the 'time_order' and 'n_stages' of 'sim'. */
void simulation_set_time_order(Simulation *sim, long time_order, long n_stages);

/* Set the 'cfl' number of 'sim'. */
void simulation_set_cfl(Simulation *sim, double cfl);

/* Set the 'max_time' of 'sim'. If a teal.restart is set, increase restart time by 'max_time'. */
void simulation_set_max_time(Simulation *sim, double max_time);

/* Set the 'output_time' of 'sim'. */
void simulation_set_output_time(Simulation *sim, double output_time);

/* Set the 'max_iter' of 'sim'. If a teal.restart is set, increase restart iter by 'max_iter'. */
void simulation_set_max_iter(Simulation *sim, long max_iter);

/* Set the 'output_iter' of 'sim'. */
void simulation_set_output_iter(Simulation *sim, long output_iter);

/* Set the abort 'variable' and 'residual' of 'sim'. */
void simulation_set_abort(Simulation *sim, long variable, double residual);

/* Restart 'sim' from file 'fname'. */
void simulation_restart(Simulation *sim, const char *fname);

/* Print a summary of 'sim'. */
void simulation_print(const Simulation *sim);

/* Run 'sim'. */
void simulation_run(Simulation *sim);

/* Compute 'sim' error of solution variable 'ivars' and user variable 'iuser'. */
double simulation_error(const Simulation *sim, long ivars, long iuser);
