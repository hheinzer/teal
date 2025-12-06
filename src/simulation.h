// Simulation couples an equation system to a time-integration loop with output/termination control.
//
// Components:
// - Time / iteration limits: stop and write every `out` interval up to `max`
// - Termination: optional residual criterion on a specific variable/component
// - Advance: explicit/implicit time integrators with Courant-like scaling and optional context
// - Output: optional VTKHDF dumps of mesh/solution at restart and on each write interval
#pragma once

#include "equations.h"

// Advance solution in time by at most `max_step` using a Courant-like scale factor; update `time`,
// write residual norm into `residual_`, and return the un-clipped step size (courant * dt_min).
typedef scalar Advance(const Equations *eqns, scalar *time, void *residual_, scalar max_step,
                       scalar courant, const void *ctx_);

typedef struct {
    scalar max;
    scalar out;
} SimulationTime;

typedef struct {
    long max;
    long out;
} SimulationIter;

typedef struct {
    const char *condition;  // variable name or "maximum" for global max residual
    long variable;          // flat index into residual array (-1 for global maximum)
    scalar threshold;
} SimulationTermination;

typedef struct {
    long time_order;
    long num_stages;
} RungeKutta;

typedef struct {
    scalar newton_tolerance;
    scalar krylov_tolerance;
    long krylov_dimension;
} NewtonKrylov;

typedef struct {
    Name name;
    scalar courant;
    const void *ctx;
    Advance *method;
} SimulationAdvance;

typedef struct {
    const Equations *eqns;
    const char *prefix;
    SimulationTime time;
    SimulationIter iter;
    SimulationTermination termination;
    SimulationAdvance advance;
} Simulation;

// Create a simulation bound to an equation system with default settings.
Simulation *simulation_create(const Equations *eqns, const char *prefix);

// Configure maximum time (physical).
void simulation_set_max_time(Simulation *sim, scalar time);

// Configure output time (physical).
void simulation_set_out_time(Simulation *sim, scalar time);

// Configure maximum iteration counts.
void simulation_set_max_iter(Simulation *sim, long iter);

// Configure output iteration counts.
void simulation_set_out_iter(Simulation *sim, long iter);

// Set convergence criterion: either "maximum" or a variable name (suffixed with -x/-y/-z).
void simulation_set_termination(Simulation *sim, const char *condition, scalar threshold);

// Select advance method, Courant factor, and optional context.
void simulation_set_advance(Simulation *sim, const char *name, scalar courant, const void *ctx_);

// Print a summary of the simulation configuration.
void simulation_summary(const Simulation *sim);

// Run the simulation loop; writes outputs per config and returns final time.
scalar simulation_run(Simulation *sim);

// Compute and print norms at a given time; pass a buffer or `0` to allocate temporaries.
void simulation_error(const Simulation *sim, scalar time, void *norm_);
