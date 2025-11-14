#include <math.h>

#include "euler.h"
#include "teal/assert.h"
#include "teal/vector.h"

static void conserved(void *variable_, const scalar *property)
{
    Euler *variable = variable_;
    scalar gamma = property[0];

    variable->momentum = vector_mul(variable->density, variable->velocity);
    variable->energy = (variable->pressure / (gamma - 1)) +
                       (vector_dot(variable->momentum, variable->velocity) / 2);
}

static void primitive(void *variable_, const scalar *property)
{
    Euler *variable = variable_;
    scalar gamma = property[0];

    variable->velocity = vector_div(variable->momentum, variable->density);
    variable->pressure =
        (gamma - 1) * (variable->energy - (vector_dot(variable->momentum, variable->velocity) / 2));
}

static scalar timestep(const void *variable_, const scalar *property, scalar volume,
                       vector projection)
{
    const Euler *variable = variable_;
    scalar gamma = property[0];

    scalar speed_of_sound = sqrt(gamma * variable->pressure / variable->density);
    scalar lambda_c_x = (fabs(variable->velocity.x) + speed_of_sound) * projection.x;
    scalar lambda_c_y = (fabs(variable->velocity.y) + speed_of_sound) * projection.y;
    scalar lambda_c_z = (fabs(variable->velocity.z) + speed_of_sound) * projection.z;
    return volume / (lambda_c_x + lambda_c_y + lambda_c_z);
}

Equations *euler_create(const Mesh *mesh)
{
    assert(mesh);

    Equations *eqns = equations_create(mesh, "euler");

    const char *variable[5] = {"density", "momentum", "energy", "velocity", "pressure"};
    Type type[5] = {SCALAR, VECTOR, SCALAR, VECTOR, SCALAR};
    equations_create_variables(eqns, variable, type, conserved, primitive, 3, 5);

    const char *property = "heat capacity ratio";
    const scalar gamma = 1.4;  // air
    equations_create_properties(eqns, &property, &gamma, 1);

    equations_set_space_order(eqns, 1);
    equations_set_timestep(eqns, timestep);
    equations_set_boundary_select(eqns, euler_boundary);
    equations_set_convective_select(eqns, euler_convective);
    equations_set_limiter(eqns, "minmod", 0);
    equations_set_convective_flux(eqns, "godunov");

    return eqns;
}
