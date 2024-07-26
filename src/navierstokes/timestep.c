#include "timestep.h"

#include <math.h>

#include "core/utils.h"
#include "navierstokes.h"

#define C 4  // constant for central spatial discretization

double timestep_visc(const Equations *eqns, const double *u, const double *projection,
                     double volume)
{
    // Blazek 2015, eq. (6.22)
    const double gamma = eqns->scalar.value[GAMMA];
    const double prandtl = eqns->scalar.value[PRANDTL];
    const double mu = eqns->scalar.value[MU];
    const double c = sqrt(gamma * u[P] / u[D]);
    const double lambdaCx = (fabs(u[U]) + c) * projection[X];
    const double lambdaCy = (fabs(u[V]) + c) * projection[Y];
    const double lambdaCz = (fabs(u[W]) + c) * projection[Z];
    const double fac = max(4 * (mu) / (3 * u[D]), gamma / u[D] * (mu / prandtl)) / volume;
    const double lambdaVx = sq(projection[X]);
    const double lambdaVy = sq(projection[Y]);
    const double lambdaVz = sq(projection[Z]);
    return volume / (lambdaCx + lambdaCy + lambdaCz + C * fac * (lambdaVx + lambdaVy + lambdaVz));
}
