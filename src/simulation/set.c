#include <string.h>

#include "advance.h"
#include "simulation.h"
#include "teal/assert.h"
#include "teal/utils.h"

void simulation_set_courant(Simulation *sim, scalar courant)
{
    assert(sim && courant > 0);
    sim->courant = courant;
}

void simulation_set_time_order(Simulation *sim, number time_order)
{
    assert(sim && 1 <= time_order && time_order <= 3);
    sim->time_order = time_order;
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
    number len = strlen(condition);
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
    error("invalid variable -- '%.*s'", (int)len, condition);
}

void simulation_set_advance(Simulation *sim, const char *name)
{
    assert(sim && name);
    strcpy(sim->advance.name, name);
    if (!strcmp(name, "explicit euler")) {
        sim->time_order = 1;
        sim->advance.method = explicit_euler;
        return;
    }
    error("invalid advance method -- '%s'", name);
}
