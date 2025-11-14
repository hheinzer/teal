#include "euler.h"
#include "teal/vector.h"

void euler_conserved(void *variable_, const scalar *property)
{
    Euler *variable = variable_;
    scalar gamma = property[0];

    variable->momentum = vector_mul(variable->density, variable->velocity);
    variable->energy = (variable->pressure / (gamma - 1)) +
                       (vector_dot(variable->momentum, variable->velocity) / 2);
}

void euler_primitive(void *variable_, const scalar *property)
{
    Euler *variable = variable_;
    scalar gamma = property[0];

    variable->velocity = vector_div(variable->momentum, variable->density);
    variable->pressure =
        (gamma - 1) * (variable->energy - (vector_dot(variable->momentum, variable->velocity) / 2));
}
