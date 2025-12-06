#include <math.h>

#include "navier_stokes.h"
#include "teal/utils.h"

scalar navier_stokes_time_step(const void *variable_, const scalar *property, scalar volume,
                               vector projection)
{
    const NavierStokes *variable = variable_;
    scalar gamma = property[NAVIER_STOKES_HEAT_CAPACITY_RATIO];
    scalar viscosity = property[NAVIER_STOKES_DYNAMIC_VISCOSITY];
    scalar prandtl = property[NAVIER_STOKES_PRANDTL];

    scalar speed_of_sound = sqrt(gamma * variable->pressure / variable->density);
    scalar lambda_c_x = (fabs(variable->velocity.x) + speed_of_sound) * projection.x;
    scalar lambda_c_y = (fabs(variable->velocity.y) + speed_of_sound) * projection.y;
    scalar lambda_c_z = (fabs(variable->velocity.z) + speed_of_sound) * projection.z;

    scalar factor = fmax(4 * viscosity / (3 * variable->density),
                         gamma * viscosity / (variable->density * prandtl)) /
                    volume;
    scalar sum_lambda_v = sq(projection.x) + sq(projection.y) + sq(projection.z);
    return volume / (lambda_c_x + lambda_c_y + lambda_c_z + 4 * factor * sum_lambda_v);
}
