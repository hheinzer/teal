#include <string.h>

#include "advance.h"
#include "simulation.h"
#include "teal/arena.h"
#include "teal/assert.h"
#include "teal/utils.h"

void simulation_set_courant(Simulation *sim, scalar courant)
{
    assert(sim && courant > 0);
    sim->courant = courant;
}

void simulation_set_max_time(Simulation *sim, scalar time)
{
    assert(sim && time > 0);
    sim->time.max = time;
}

void simulation_set_output_time(Simulation *sim, scalar time)
{
    assert(sim && time > 0);
    sim->time.output = time;
}

void simulation_set_max_iter(Simulation *sim, number iter)
{
    assert(sim && iter > 0);
    sim->iter.max = iter;
}

void simulation_set_output_iter(Simulation *sim, number iter)
{
    assert(sim && iter > 0);
    sim->iter.output = iter;
}

static number component_index(char chr)
{
    switch (chr) {
        case 'x': return 0;
        case 'y': return 1;
        case 'z': return 2;
        default: error("invalid component index -- '%c'", chr);
    }
}

void simulation_set_termination(Simulation *sim, const char *condition, scalar residual)
{
    assert(sim && condition && residual > 0);

    sim->termination.condition = condition;
    sim->termination.residual = residual;
    if (!strcmp(condition, "maximum")) {
        sim->termination.variable = -1;
        return;
    }

    number idx = 0;
    int len = strlen(condition);
    if (condition[len - 1] == '-') {
        idx = component_index(condition[len]);
        len -= 2;
    }
    else if (condition[len - 2] == '-') {
        idx = (component_index(condition[len - 1]) * 3) + component_index(condition[len]);
        len -= 3;
    }

    number num = sim->eqns->variables.num;
    Name *name = sim->eqns->variables.name;
    Type *type = sim->eqns->variables.type;

    for (number j = 0, i = 0; i < num; j += type[i++]) {
        if (!strncmp(name[i], condition, len)) {
            sim->termination.variable = j + idx;
            return;
        }
    }
    error("invalid variable -- '%.*s'", len, condition);
}

void simulation_set_advance(Simulation *sim, const char *name)
{
    assert(sim && name);
    strcpy(sim->advance.name, name);
    if (!strcmp(name, "euler")) {
        sim->advance.method = euler;
        return;
    }
    if (!strcmp(name, "midpoint")) {
        sim->advance.method = midpoint;
        return;
    }
    if (!strcmp(name, "heun")) {
        sim->advance.method = heun;
        return;
    }
    if (!strcmp(name, "ralston")) {
        sim->advance.method = ralston;
        return;
    }
    if (!strcmp(name, "ssprk2")) {
        sim->advance.method = ssprk2;
        return;
    }
    if (!strcmp(name, "ssprk3")) {
        sim->advance.method = ssprk3;
        return;
    }
    if (!strcmp(name, "rk3")) {
        sim->advance.method = rk3;
        return;
    }
    if (!strcmp(name, "rk4")) {
        sim->advance.method = rk4;
        return;
    }
    if (!strncmp(name, "lserk", 5)) {
        Lserk *context = arena_malloc(1, sizeof(*context));
        context->time_order = name[5] - '0';
        context->num_stages = name[6] - '0';
        assert(1 <= context->time_order && context->time_order <= 3);
        assert(1 <= context->num_stages && context->num_stages <= 6);
        sim->advance.context_ = context;
        sim->advance.method = lserk;
        return;
    }
    error("invalid advance method -- '%s'", name);
}
