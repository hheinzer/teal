#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "simulation.h"
#include "teal/arena.h"
#include "teal/utils.h"

void simulation_error(const Simulation *sim, scalar time, void *norm_)
{
    assert(sim);
    Arena save = arena_save();

    long num = sim->eqns->variables.num;
    long *size = sim->eqns->variables.size;
    Name *name = sim->eqns->variables.name;

    int width = 0;
    for (long i = 0; i < num; i++) {
        width = lmax(width, strlen(name[i]));
    }

    scalar *norm = norm_;
    if (!norm) {
        long stride = sim->eqns->variables.stride;
        norm = arena_malloc(stride, sizeof(*norm));
    }
    equations_norm(sim->eqns, time, norm);

    println("Simulation error");
    for (long j = 0, i = 0; i < num; j += size[i++]) {
        long pos = 0;
        char buf[128];
        for (long k = 0; k < size[i]; k++) {
            pos += sprintf(&buf[pos], " %g", norm[j + k]);
        }
        println("\t %-*s :%s", width, name[i], buf);
    }

    arena_load(save);
}
