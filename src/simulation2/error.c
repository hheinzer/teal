#include <assert.h>
#include <stdio.h>

#include "equations2.h"
#include "simulation2.h"
#include "teal2.h"

void simulation2_error(const Simulation *sim, double time, void *norm_)
{
    assert(sim);

    int num = sim->eqns->primitive.num;
    int *dimension = sim->eqns->primitive.dimension;
    String *name = sim->eqns->primitive.name;

    double *norm = norm_;
    if (!norm) {
        norm = teal2_calloc(sim->eqns->primitive.stride, (int)sizeof(*norm));
    }
    equations2_norm(sim->eqns, norm, time);

    teal2_print("Simulation error");
    for (int i = 0, j = 0; i < num; j += dimension[i++]) {
        int pos = 0;
        char buf[128];
        for (int k = 0; k < dimension[i]; k++) {
            pos += sprintf(&buf[pos], " %g", norm[j + k]);
        }
        teal2_print("\t %-20s :%s", name[i], buf);
    }

    if (!norm_) {
        teal2_free(norm);
    }
}
