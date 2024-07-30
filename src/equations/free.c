#include <stdlib.h>

#include "equations.h"

void equations_free(Equations *eqns)
{
    free(eqns->scalar.name);
    free(eqns->scalar.value);

    free(eqns->vars.name);
    free(eqns->vars.u);
    free(eqns->vars.dudx);
    free(eqns->vars.dudt);
    free(eqns->vars.dim);

    free(eqns->limiter.eps2);

    free(eqns->user.name);
    free(eqns->user.dim);

    free(eqns->bc.name);
    free(eqns->bc.state);
    free(eqns->bc.apply);
    free(eqns->bc.func);

    free(eqns->sync.buf);
    free(eqns->sync.recv);
    free(eqns->sync.send);

    *eqns = (Equations){0};
}
