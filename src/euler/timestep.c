#include "timestep.h"

#include <math.h>

#include "euler.h"

double timestep_conv(const Equations *eqns, const double *u, const double *projection,
                     double volume)
{
    // Blazek 2015, eq. (6.22)
    const double gamma = eqns->scalar.value[GAMMA];
    const double c = sqrt(gamma * u[P] / u[D]);
    const double lambdaCx = (fabs(u[U]) + c) * projection[X];
    const double lambdaCy = (fabs(u[V]) + c) * projection[Y];
    const double lambdaCz = (fabs(u[W]) + c) * projection[Z];
    return volume / (lambdaCx + lambdaCy + lambdaCz);
}
