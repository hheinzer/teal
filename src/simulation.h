#pragma once

#include "equations.h"

typedef struct {
    scalar max;
    scalar output;
} SimulationTime;

typedef struct {
    number max;
    number output;
} SimulationIter;

typedef struct {
    const char *condition;
    number variable;
    scalar residual;
} SimulationTermination;

typedef scalar Advance(const Equations *eqns, scalar *time, void *residual_, scalar courant,
                       scalar max_timestep, const void *context_);

typedef struct {
    Name name;
    const void *context_;
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

void simulation_set_max_iter(Simulation *sim, number iter);

void simulation_set_output_iter(Simulation *sim, number iter);

void simulation_set_termination(Simulation *sim, const char *condition, scalar residual);

void simulation_set_advance(Simulation *sim, const char *name);

void simulation_summary(const Simulation *sim);

scalar simulation_run(Simulation *sim);

scalar simulation_error(const Simulation *sim, scalar time);
