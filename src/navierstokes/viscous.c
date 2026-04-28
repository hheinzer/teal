#include <string.h>

#include "navierstokes.h"

typedef struct {
    Vector density;
    Matrix velocity;
    Vector pressure;
} Gradient;

static void central(void *flux_, const void *primitive_, const void *gradient_,
                    const double *property, Matrix basis)
{
    NavierStokesConserved *flux = flux_;
    const NavierStokesPrimitive *primitive = primitive_;
    const Gradient *gradient = gradient_;
    double gamma = property[NAVIERSTOKES_HEAT_CAPACITY_RATIO];
    double viscosity = property[NAVIERSTOKES_DYNAMIC_VISCOSITY];
    double prandtl = property[NAVIERSTOKES_PRANDTL];

    Matrix tau;
    double divergence = gradient->velocity.x.x + gradient->velocity.y.y + gradient->velocity.z.z;
    tau.x.x = 2 * viscosity * (gradient->velocity.x.x - (divergence / 3));
    tau.y.y = 2 * viscosity * (gradient->velocity.y.y - (divergence / 3));
    tau.z.z = 2 * viscosity * (gradient->velocity.z.z - (divergence / 3));
    tau.y.x = tau.x.y = viscosity * (gradient->velocity.x.y + gradient->velocity.y.x);
    tau.z.x = tau.x.z = viscosity * (gradient->velocity.x.z + gradient->velocity.z.x);
    tau.z.y = tau.y.z = viscosity * (gradient->velocity.y.z + gradient->velocity.z.y);

    double factor = -viscosity * gamma / ((gamma - 1) * prandtl * sq(primitive->density));
    Vector heat_flux =
        vector_mul(factor, vector_sub(vector_mul(primitive->density, gradient->pressure),
                                      vector_mul(primitive->pressure, gradient->density)));
    Vector theta = vector_sub(matrix_vector(tau, primitive->velocity), heat_flux);

    flux->density = 0;
    flux->momentum = matrix_vector(tau, basis.x);
    flux->energy = vector_dot(theta, basis.x);
}

Viscous *navierstokes_viscous(const char *name)
{
    if (!strcmp(name, "central")) {
        return central;
    }
    teal_error("invalid viscous flux (%s)", name);
}
