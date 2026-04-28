#include <assert.h>

#include "equations.h"
#include "teal.h"

static void destroy_variables(EquationsVariables *variables)
{
    teal_free(variables->name);
    teal_free(variables->dimension);
    teal_free(variables->data);
}

void equations_destroy(Equations *eqns)
{
    assert(eqns);

    destroy_variables(&eqns->primitive);
    destroy_variables(&eqns->conserved);
    destroy_variables(&eqns->reference);

    teal_free(eqns->properties.name);
    teal_free(eqns->properties.data);

    teal_free(eqns->boundary.name);
    teal_free(eqns->boundary.compute);
    teal_free(eqns->boundary.condition);
    teal_free(eqns->boundary.context);

    teal_free(eqns->limiter.parameter);

    teal_free(eqns);
}
