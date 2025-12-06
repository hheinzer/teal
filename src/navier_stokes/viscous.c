#include <string.h>

#include "navier_stokes.h"
#include "teal/utils.h"

typedef struct {
    scalar density;
    vector momentum;
    scalar energy;
} Conserved;

typedef struct {
    vector density;
    matrix momentum __attribute((unused));
    vector energy __attribute((unused));
    matrix velocity;
    vector pressure;
} Gradient;

static void central(void *flux_, const void *variable_, const void *gradient_,
                    const scalar *property, const Basis *basis)
{
    Conserved *flux = flux_;
    const NavierStokes *variable = variable_;
    const Gradient *gradient = gradient_;
    scalar gamma = property[NAVIER_STOKES_HEAT_CAPACITY_RATIO];
    scalar viscosity = property[NAVIER_STOKES_DYNAMIC_VISCOSITY];
    scalar prandtl = property[NAVIER_STOKES_PRANDTL];
    vector normal = basis->normal;

    scalar divergence = gradient->velocity.x.x + gradient->velocity.y.y + gradient->velocity.z.z;
    scalar tau_xx = 2 * viscosity * (gradient->velocity.x.x - (divergence / 3));
    scalar tau_yy = 2 * viscosity * (gradient->velocity.y.y - (divergence / 3));
    scalar tau_zz = 2 * viscosity * (gradient->velocity.z.z - (divergence / 3));
    scalar tau_xy = viscosity * (gradient->velocity.x.y + gradient->velocity.y.x);
    scalar tau_xz = viscosity * (gradient->velocity.x.z + gradient->velocity.z.x);
    scalar tau_yz = viscosity * (gradient->velocity.y.z + gradient->velocity.z.y);

    scalar factor = -viscosity * gamma / ((gamma - 1) * prandtl * sq(variable->density));
    scalar heat_flux_x = factor * ((gradient->pressure.x * variable->density) -
                                   (variable->pressure * gradient->density.x));
    scalar heat_flux_y = factor * ((gradient->pressure.y * variable->density) -
                                   (variable->pressure * gradient->density.y));
    scalar heat_flux_z = factor * ((gradient->pressure.z * variable->density) -
                                   (variable->pressure * gradient->density.z));

    scalar theta_x = (variable->velocity.x * tau_xx) + (variable->velocity.y * tau_xy) +
                     (variable->velocity.z * tau_xz) - heat_flux_x;
    scalar theta_y = (variable->velocity.x * tau_xy) + (variable->velocity.y * tau_yy) +
                     (variable->velocity.z * tau_yz) - heat_flux_y;
    scalar theta_z = (variable->velocity.x * tau_xz) + (variable->velocity.y * tau_yz) +
                     (variable->velocity.z * tau_zz) - heat_flux_z;

    flux->density = 0;
    flux->momentum.x = (normal.x * tau_xx) + (normal.y * tau_xy) + (normal.z * tau_xz);
    flux->momentum.y = (normal.x * tau_xy) + (normal.y * tau_yy) + (normal.z * tau_yz);
    flux->momentum.z = (normal.x * tau_xz) + (normal.y * tau_yz) + (normal.z * tau_zz);
    flux->energy = (normal.x * theta_x) + (normal.y * theta_y) + (normal.z * theta_z);
}

Viscous *navier_stokes_viscous(const char *name)
{
    if (!strcmp(name, "central")) {
        return central;
    }
    error("invalid viscous flux (%s)", name);
}
