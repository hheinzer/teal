#include <assert.h>
#include <math.h>

#include "euler2.h"

double euler2_timestep(const void *primitive_, const double *property, double volume,
                       Vector projection)
{
    const EulerPrimitive *primitive = primitive_;
    double gamma = property[EULER_HEAT_CAPACITY_RATIO];

    double speed_of_sound = sqrt(gamma * primitive->pressure / primitive->density);

    double lambda_c_x = (fabs(primitive->velocity.x) + speed_of_sound) * projection.x;
    double lambda_c_y = (fabs(primitive->velocity.y) + speed_of_sound) * projection.y;
    double lambda_c_z = (fabs(primitive->velocity.z) + speed_of_sound) * projection.z;

    return volume / (lambda_c_x + lambda_c_y + lambda_c_z);
}

void euler2_primitive(void *primitive_, const void *conserved_, const double *property)
{
    EulerPrimitive *primitive = primitive_;
    const EulerConserved *conserved = conserved_;
    double gamma = property[EULER_HEAT_CAPACITY_RATIO];

    primitive->density = fmax(1e-8, conserved->density);
    primitive->velocity = vector2_div(conserved->momentum, primitive->density);
    primitive->pressure =
        fmax(1e-8, (gamma - 1) * (conserved->energy -
                                  (vector2_dot(conserved->momentum, primitive->velocity) / 2)));
}

void euler2_conserved(void *conserved_, const void *primitive_, const double *property)
{
    EulerConserved *conserved = conserved_;
    const EulerPrimitive *primitive = primitive_;
    double gamma = property[EULER_HEAT_CAPACITY_RATIO];

    conserved->density = primitive->density;
    conserved->momentum = vector2_mul(conserved->density, primitive->velocity);
    conserved->energy = (primitive->pressure / (gamma - 1)) +
                        (vector2_dot(conserved->momentum, primitive->velocity) / 2);
}

Equations *euler2_create(const Mesh *mesh)
{
    assert(mesh);

    Equations *eqns =
        equations2_create(mesh, "euler", euler2_timestep, euler2_convective, 0, euler2_boundary);

    const char *primitive[] = {"density", "velocity", "pressure"};
    equations2_create_primitive(eqns, primitive, (int[]){1, 3, 1}, euler2_primitive, 3);

    const char *conserved[] = {"_density", "momentum", "energy"};
    equations2_create_conserved(eqns, conserved, (int[]){1, 3, 1}, euler2_conserved, 3);

    const char *property[] = {"heat capacity ratio"};
    equations2_create_properties(eqns, property, (double[]){1.4}, 1);

    equations2_set_space_order(eqns, 2);
    equations2_set_convective_flux(eqns, "hllc");
    equations2_set_limiter(eqns, "minmod", 0);

    return eqns;
}
