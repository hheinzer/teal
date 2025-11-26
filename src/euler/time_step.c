#include <math.h>

#include "euler.h"

scalar euler_time_step(const void *variable_, const scalar *property, scalar volume,
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
