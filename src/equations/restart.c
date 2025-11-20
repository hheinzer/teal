#include <assert.h>

#include "equations.h"
#include "teal/option.h"

void equations_restart(const Equations *eqns, scalar *time, long *index)
{
    assert(eqns && time && index);
    if (!option.restart) {
        *time = 0;
        *index = 0;
        return;
    }
    // TODO
    (void)eqns;
}
