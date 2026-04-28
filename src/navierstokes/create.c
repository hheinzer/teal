#include <assert.h>
#include <math.h>

#include "euler.h"
#include "navierstokes.h"

double navierstokes_timestep(const void *primitive_, const double *property, double volume,
                             Vector projection)
{
    const NavierStokesPrimitive *primitive = primitive_;
    double gamma = property[NAVIERSTOKES_HEAT_CAPACITY_RATIO];
    double viscosity = property[NAVIERSTOKES_DYNAMIC_VISCOSITY];
    double prandtl = property[NAVIERSTOKES_PRANDTL];

    double speed_of_sound = sqrt(gamma * primitive->pressure / primitive->density);

    double lambda_c_x = (fabs(primitive->velocity.x) + speed_of_sound) * projection.x;
    double lambda_c_y = (fabs(primitive->velocity.y) + speed_of_sound) * projection.y;
    double lambda_c_z = (fabs(primitive->velocity.z) + speed_of_sound) * projection.z;

    double lambda_v = fmax(4 * viscosity / (3 * primitive->density),
                           gamma * viscosity / (primitive->density * prandtl));

    return volume / (lambda_c_x + lambda_c_y + lambda_c_z +
                     (4 * lambda_v * vector_norm2(projection) / volume));
}

Equations *navierstokes_create(const Mesh *mesh)
{
    assert(mesh);

    Equations *eqns =
        equations_create(mesh, "navier-stokes", navierstokes_timestep, euler_convective,
                         navierstokes_viscous, navierstokes_boundary);

    const char *primitive[] = {"density", "velocity", "pressure"};
    equations_create_primitive(eqns, primitive, (int[]){1, 3, 1}, euler_primitive, 3);

    const char *conserved[] = {"_density", "momentum", "energy"};
    equations_create_conserved(eqns, conserved, (int[]){1, 3, 1}, euler_conserved, 3);

    const char *property[] = {"heat capacity ratio", "dynamic viscosity", "prandtl number"};
    equations_create_properties(eqns, property, (double[]){1.4, 0, 0.71}, 3);

    equations_set_space_order(eqns, 2);
    equations_set_convective_flux(eqns, "hllc");
    equations_set_viscous_flux(eqns, "central");
    equations_set_limiter(eqns, "minmod", 0);

    return eqns;
}
