#include "equations.h"
#include "teal/memory.h"

void equations_free(Equations **eqns)
{
    memory_free(&(*eqns)->timestep.value);

    memory_free(&(*eqns)->vars.name);
    memory_free(&(*eqns)->vars.u);
    memory_free(&(*eqns)->vars.dudx);
    memory_free(&(*eqns)->vars.dudt);
    memory_free(&(*eqns)->vars.dim);

    memory_free(&(*eqns)->scalar.name);
    memory_free(&(*eqns)->scalar.value);

    memory_free(&(*eqns)->user.name);
    memory_free(&(*eqns)->user.dim);

    memory_free(&(*eqns)->limiter.eps2);

    memory_free(&(*eqns)->bc.name);
    memory_free(&(*eqns)->bc.state);
    memory_free(&(*eqns)->bc.apply);
    memory_free(&(*eqns)->bc.compute);

    memory_free(&(*eqns)->sync.buf);
    memory_free(&(*eqns)->sync.recv);
    memory_free(&(*eqns)->sync.send);

    memory_free(eqns);
}
