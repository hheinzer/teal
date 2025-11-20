#include <assert.h>
#include <string.h>

#include "advance.h"
#include "simulation.h"
#include "teal/arena.h"
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

void simulation_set_max_iter(Simulation *sim, long iter)
{
    assert(sim && iter > 0);
    sim->iter.max = iter;
}

void simulation_set_output_iter(Simulation *sim, long iter)
{
    assert(sim && iter > 0);
    sim->iter.output = iter;
}

static long component_index(char chr)
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

    long idx = 0;
    int len = strlen(condition);
    if (condition[len - 1] == '-') {
        idx = component_index(condition[len]);
        len -= 2;
    }
    else if (condition[len - 2] == '-') {
        idx = (component_index(condition[len - 1]) * 3) + component_index(condition[len]);
        len -= 3;
    }

    long num = sim->eqns->variables.num;
    long *size = sim->eqns->variables.size;
    Name *name = sim->eqns->variables.name;

    for (long j = 0, i = 0; i < num; j += size[i++]) {
        if (!strncmp(name[i], condition, len)) {
            sim->termination.variable = j + idx;
            return;
        }
    }
    error("invalid variable -- '%.*s'", len, condition);
}

void simulation_set_advance(Simulation *sim, const char *name, const void *ctx_)
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
    if (!strcmp(name, "lserk")) {
        const RungeKutta *ctx = ctx_;
        if (!(1 <= ctx->time_order && ctx->time_order <= 3)) {
            error("invalid time order -- '%ld'", ctx->time_order);
        }
        if (!(1 <= ctx->num_stages && ctx->num_stages <= 6)) {
            error("invalid long of stages -- '%ld'", ctx->num_stages);
        }
        sim->advance.ctx = arena_memdup(ctx, 1, sizeof(*ctx));
        sim->advance.method = lserk;
        return;
    }
    if (!strcmp(name, "implicit euler")) {
        const NewtonKrylov *ctx = ctx_;
        if (ctx->newton_tolerance <= 0) {
            error("invalid newton tolerance -- '%g'", ctx->newton_tolerance);
        }
        if (ctx->krylov_tolerance <= 0) {
            error("invalid krylov tolerance -- '%g'", ctx->krylov_tolerance);
        }
        if (ctx->krylov_dimension <= 0) {
            error("invalid krylov dimension -- '%ld'", ctx->krylov_dimension);
        }
        sim->advance.ctx = arena_memdup(ctx, 1, sizeof(*ctx));
        sim->advance.method = implicit_euler;
        return;
    }
    error("invalid advance method -- '%s'", name);
}
