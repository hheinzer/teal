#include <assert.h>

#include "equations2.h"
#include "teal2.h"

static void destroy_variables(EquationsVariables *variables)
{
    teal2_free(variables->name);
    teal2_free(variables->dimension);
    teal2_free(variables->data);
}

void equations2_destroy(Equations *eqns)
{
    assert(eqns);

    destroy_variables(&eqns->primitive);
    destroy_variables(&eqns->conserved);
    destroy_variables(&eqns->reference);

    teal2_free(eqns->properties.name);
    teal2_free(eqns->properties.data);

    teal2_free(eqns->boundary.name);
    teal2_free(eqns->boundary.compute);
    teal2_free(eqns->boundary.condition);
    teal2_free(eqns->boundary.context);

    teal2_free(eqns->limiter.parameter);

    teal2_free(eqns);
}
