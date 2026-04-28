#pragma once

#include "equations.h"

// Advance the solution by one time step; returns the step size taken.
typedef double Advance(const Equations *eqns, double *time, void *residual, double max_step,
                       double courant, const void *context);

typedef struct {
    double max;
    double out;
} SimulationTime;

typedef struct {
    int max;
    int out;
} SimulationIter;

typedef struct {
    String condition;
    int variable;
    double threshold;
} SimulationTermination;

typedef struct {
    int time_order;
    int num_stages;
} RungeKutta;

typedef struct {
    double newton_tolerance;
    double krylov_tolerance;
    int krylov_dimension;
} NewtonKrylov;

typedef struct {
    String name;
    double courant;
    void *context;
    Advance *method;
} SimulationAdvance;

typedef struct {
    const Equations *eqns;
    String prefix;
    SimulationTime time;
    SimulationIter iter;
    SimulationTermination termination;
    SimulationAdvance advance;
} Simulation;

// Create a simulation driver bound to an equation system with an output file prefix.
Simulation *simulation_create(const Equations *eqns, const char *prefix);

// Set the maximum simulation time (physical).
void simulation_set_max_time(Simulation *sim, double time);

// Set the output interval in simulation time (physical).
void simulation_set_out_time(Simulation *sim, double time);

// Set the maximum number of time steps.
void simulation_set_max_iter(Simulation *sim, int iter);

// Set the output interval in time steps.
void simulation_set_out_iter(Simulation *sim, int iter);

// Set a convergence termination condition on a residual variable.
void simulation_set_termination(Simulation *sim, const char *condition, double threshold);

// Select the time advancement scheme by name with a Courant number.
void simulation_set_advance(Simulation *sim, const char *name, double courant, const void *context);

// Print a summary of the simulation configuration.
void simulation_summary(const Simulation *sim);

// Run the simulation to completion and return the final time.
double simulation_run(Simulation *sim);

// Compute the L2 error norm against the reference solution at the given time.
void simulation_error(const Simulation *sim, double time, void *norm);

// Release all simulation resources.
void simulation_destroy(Simulation *sim);
