#include "timestep.h"

#include <math.h>

#include "navierstokes.h"
#include "teal/utils.h"

double navierstokes_timestep(const Equations *eqns, const double *u, const Vector3d projection,
                             double volume)
{
    // Blazek 2015, eq. (6.22)
    const double gamma = eqns->scalar.value[GAMMA];
    const double mu = eqns->scalar.value[MU];
    const double prandtl = eqns->scalar.value[PRANDTL];

    const double c = sqrt(gamma * u[P] / u[D]);
    const double lambdaCx = (fabs(u[U]) + c) * projection[X];
    const double lambdaCy = (fabs(u[V]) + c) * projection[Y];
    const double lambdaCz = (fabs(u[W]) + c) * projection[Z];

    const double fac = max(4 * (mu) / (3 * u[D]), gamma / u[D] * (mu / prandtl)) / volume;
    const double lambdaVx = sq(projection[X]);
    const double lambdaVy = sq(projection[Y]);
    const double lambdaVz = sq(projection[Z]);

    return volume / (lambdaCx + lambdaCy + lambdaCz + 4 * fac * (lambdaVx + lambdaVy + lambdaVz));
}
