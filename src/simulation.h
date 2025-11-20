#pragma once

#include "equations.h"

typedef scalar Advance(const Equations *eqns, scalar *time, void *residual_, scalar courant,
                       scalar max_step, const void *ctx_);

typedef struct {
    scalar max;
    scalar output;
} SimulationTime;

typedef struct {
    long max;
    long output;
} SimulationIter;

typedef struct {
    const char *condition;
    long variable;
    scalar residual;
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
    void *ctx;
    Advance *method;
} SimulationAdvance;

typedef struct {
    const Equations *eqns;
    const char *prefix;
    scalar courant;
    SimulationTime time;
    SimulationIter iter;
    SimulationTermination termination;
    SimulationAdvance advance;
} Simulation;

Simulation *simulation_create(const Equations *eqns, const char *prefix);

void simulation_set_courant(Simulation *sim, scalar courant);

void simulation_set_max_time(Simulation *sim, scalar time);

void simulation_set_output_time(Simulation *sim, scalar time);

void simulation_set_max_iter(Simulation *sim, long iter);

void simulation_set_output_iter(Simulation *sim, long iter);

void simulation_set_termination(Simulation *sim, const char *condition, scalar residual);

void simulation_set_advance(Simulation *sim, const char *name, const void *ctx_);

void simulation_summary(const Simulation *sim);

scalar simulation_run(Simulation *sim);

void simulation_error(const Simulation *sim, scalar time, void *norm_);
