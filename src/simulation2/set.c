#include <assert.h>
#include <string.h>

#include "advance.h"
#include "simulation2.h"
#include "teal2.h"

void simulation2_set_max_time(Simulation *sim, double time)
{
    assert(sim && time > 0);
    sim->time.max = time;
}

void simulation2_set_out_time(Simulation *sim, double time)
{
    assert(sim && time > 0);
    sim->time.out = time;
}

void simulation2_set_max_iter(Simulation *sim, int iter)
{
    assert(sim && iter > 0);
    sim->iter.max = iter;
}

void simulation2_set_out_iter(Simulation *sim, int iter)
{
    assert(sim && iter > 0);
    sim->iter.out = iter;
}

static int component_index(char chr)
{
    switch (chr) {
        case 'x': return 0;
        case 'y': return 1;
        case 'z': return 2;
        default: teal2_error("invalid component index (%c)", chr);
    }
}

void simulation2_set_termination(Simulation *sim, const char *condition, double threshold)
{
    assert(sim && condition && threshold > 0);

    strcpy(sim->termination.condition, condition);
    sim->termination.threshold = threshold;

    if (!strcmp(condition, "maximum")) {
        sim->termination.variable = -1;
        return;
    }

    int idx = 0;
    size_t len = strlen(condition);
    if (len > 2 && condition[len - 2] == '-') {
        idx = component_index(condition[len - 1]);
        len -= 2;
    }

    int num = sim->eqns->conserved.num;
    int *dimension = sim->eqns->conserved.dimension;
    String *name = sim->eqns->conserved.name;

    for (int i = 0, j = 0; i < num; j += dimension[i++]) {
        if (!strncmp(name[i], condition, len)) {
            sim->termination.variable = j + idx;
            return;
        }
    }
    teal2_error("invalid termination condition (%s)", condition);
}

void simulation2_set_advance(Simulation *sim, const char *name, double courant,
                             const void *context_)
{
    assert(sim && name && courant > 0);

    strcpy(sim->advance.name, name);
    sim->advance.courant = courant;

    teal2_free(sim->advance.context);
    sim->advance.context = 0;

    if (!strcmp(name, "euler")) {
        sim->advance.method = euler2;
        return;
    }
    if (!strcmp(name, "lserk")) {
        RungeKutta *context = teal2_calloc(1, sizeof(*context));
        if (context_) {
            memcpy(context, context_, sizeof(*context));
        }
        else {
            context->time_order = sim->eqns->space_order;
            context->num_stages = sim->eqns->space_order + 1;
        }
        if (!(1 <= context->time_order && context->time_order <= 3)) {
            teal2_error("invalid time order (%d)", context->time_order);
        }
        if (!(1 <= context->num_stages && context->num_stages <= 6)) {
            teal2_error("invalid number of stages (%d)", context->num_stages);
        }
        sim->advance.context = context;
        sim->advance.method = lserk2;
        return;
    }
    if (!strcmp(name, "implicit euler")) {
        NewtonKrylov *context = teal2_calloc(1, sizeof(*context));
        if (context_) {
            memcpy(context, context_, sizeof(*context));
        }
        else {
            context->newton_tolerance = 0.1;
            context->krylov_tolerance = 0.3;
            context->krylov_dimension = 32;
        }
        if (context->newton_tolerance <= 0) {
            teal2_error("invalid newton tolerance (%g)", context->newton_tolerance);
        }
        if (context->krylov_tolerance <= 0) {
            teal2_error("invalid krylov tolerance (%g)", context->krylov_tolerance);
        }
        if (context->krylov_dimension <= 0) {
            teal2_error("invalid krylov dimension (%d)", context->krylov_dimension);
        }
        sim->advance.context = context;
        sim->advance.method = implicit_euler2;
        return;
    }
    teal2_error("invalid advance method (%s)", name);
}
