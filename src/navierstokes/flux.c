#include "flux.h"

#include <stdlib.h>
#include <string.h>

#include "navierstokes.h"
#include "teal/utils.h"

static ViscFlux central;

ViscFlux *navierstokes_visc_flux(const char *name)
{
    if (!strcmp(name, "central")) return central;
    abort();
}

static void central(const Equations *eqns, const Matrix3d b, const double *um,
                    const Vector3d *dudxm, double *f)
{
    const double *n = b[X];
    const double gamma = eqns->scalar.value[GAMMA];
    const double mu = eqns->scalar.value[MU];
    const double prandtl = eqns->scalar.value[PRANDTL];

    const double divU = dudxm[U][X] + dudxm[V][Y] + dudxm[W][Z];
    const double tauXX = 2 * mu * (dudxm[U][X] - divU / 3);
    const double tauYY = 2 * mu * (dudxm[V][Y] - divU / 3);
    const double tauZZ = 2 * mu * (dudxm[W][Z] - divU / 3);
    const double tauXY = mu * (dudxm[U][Y] + dudxm[V][X]);
    const double tauXZ = mu * (dudxm[U][Z] + dudxm[W][X]);
    const double tauYZ = mu * (dudxm[V][Z] + dudxm[W][Y]);

    const double fac = -mu * gamma / ((gamma - 1) * prandtl * sq(um[D]));
    const double qX = fac * (dudxm[P][X] * um[D] - um[P] * dudxm[D][X]);
    const double qY = fac * (dudxm[P][Y] * um[D] - um[P] * dudxm[D][Y]);
    const double qZ = fac * (dudxm[P][Z] * um[D] - um[P] * dudxm[D][Z]);

    const double thetaX = um[U] * tauXX + um[V] * tauXY + um[W] * tauXZ - qX;
    const double thetaY = um[U] * tauXY + um[V] * tauYY + um[W] * tauYZ - qY;
    const double thetaZ = um[U] * tauXZ + um[V] * tauYZ + um[W] * tauZZ - qZ;

    f[D] = 0;
    f[DU] = n[X] * tauXX + n[Y] * tauXY + n[Z] * tauXZ;
    f[DV] = n[X] * tauXY + n[Y] * tauYY + n[Z] * tauYZ;
    f[DW] = n[X] * tauXZ + n[Y] * tauYZ + n[Z] * tauZZ;
    f[DE] = n[X] * thetaX + n[Y] * thetaY + n[Z] * thetaZ;
}
