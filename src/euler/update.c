#include <math.h>

#include "euler.h"

void euler_conserved(void *variable_, const scalar *property)
{
    Euler *variable = variable_;
    scalar gamma = property[0];

    variable->momentum.x = variable->density * variable->velocity.x;
    variable->momentum.y = variable->density * variable->velocity.y;
    variable->momentum.z = variable->density * variable->velocity.z;
    variable->energy =
        (variable->pressure / (gamma - 1)) + (((variable->momentum.x * variable->velocity.x) +
                                               (variable->momentum.y * variable->velocity.y) +
                                               (variable->momentum.z * variable->velocity.z)) /
                                              2);
}

void euler_primitive(void *variable_, const scalar *property)
{
    Euler *variable = variable_;
    scalar gamma = property[0];

    static const scalar min_density = 1e-8;
    variable->density = fmax(min_density, variable->density);

    scalar factor = 1 / variable->density;
    variable->velocity.x = factor * variable->momentum.x;
    variable->velocity.y = factor * variable->momentum.y;
    variable->velocity.z = factor * variable->momentum.z;
    variable->pressure =
        (gamma - 1) * (variable->energy - ((variable->momentum.x * variable->velocity.x) +
                                           (variable->momentum.y * variable->velocity.y) +
                                           (variable->momentum.z * variable->velocity.z)) /
                                              2);

    static const scalar min_pressure = 1e-8;
    variable->pressure = fmax(min_pressure, variable->pressure);
}
