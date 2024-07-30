#include <string.h>

#include "core/memory.h"
#include "core/sync.h"
#include "core/utils.h"
#include "navierstokes.h"

void navierstokes_body_force(const Equations *eqns, double *force)
{
    const long n_vars = eqns->n_vars;
    const long n_entities = eqns->mesh->n_entities;
    const ALIAS(j_face, eqns->mesh->entity.j_face);
    const ALIAS(cell, eqns->mesh->face.cell);
    const ALIAS(n, eqns->mesh->face.normal);
    const ALIAS(area, eqns->mesh->face.area);
    const double(*u)[n_vars] = (void *)eqns->vars.u;
    const double(*dudx)[n_vars][N_DIMS] = (void *)eqns->vars.dudx;
    const double mu = eqns->scalar.value[MU];
    const double volume = eqns->mesh->volume;

    memory_setzero(force, N_DIMS, sizeof(*force));
    for (long e = 0; e < n_entities; ++e) {
        if (strcmp(eqns->bc.name[e], "wall")) continue;
        for (long j = j_face[e]; j < j_face[e + 1]; ++j) {
            const long l = cell[j][L], r = cell[j][R];

            const double pw = 0.5 * (u[l][P] + u[r][P]);
            const double dudxw[N_DIMS][N_DIMS] = {
                {
                    0.5 * (dudx[l][U][X] + dudx[r][U][X]),
                    0.5 * (dudx[l][U][Y] + dudx[r][U][Y]),
                    0.5 * (dudx[l][U][Z] + dudx[r][U][Z]),
                },
                {
                    0.5 * (dudx[l][V][X] + dudx[r][V][X]),
                    0.5 * (dudx[l][V][Y] + dudx[r][V][Y]),
                    0.5 * (dudx[l][V][Z] + dudx[r][V][Z]),
                },
                {
                    0.5 * (dudx[l][W][X] + dudx[r][W][X]),
                    0.5 * (dudx[l][W][Y] + dudx[r][W][Y]),
                    0.5 * (dudx[l][W][Z] + dudx[r][W][Z]),
                },
            };

            const double divU = dudxw[X][X] + dudxw[Y][Y] + dudxw[Z][Z];
            const double tauXX = 2 * mu * (dudxw[X][X] - divU / 3);
            const double tauYY = 2 * mu * (dudxw[Y][Y] - divU / 3);
            const double tauZZ = 2 * mu * (dudxw[Z][Z] - divU / 3);
            const double tauXY = mu * (dudxw[X][Y] + dudxw[Y][X]);
            const double tauXZ = mu * (dudxw[X][Z] + dudxw[Z][X]);
            const double tauYZ = mu * (dudxw[Y][Z] + dudxw[Z][Y]);

            force[X] +=
                (pw * n[j][X] - (n[j][X] * tauXX + n[j][Y] * tauXY + n[j][Z] * tauXZ)) * area[j];
            force[Y] +=
                (pw * n[j][Y] - (n[j][X] * tauXY + n[j][Y] * tauYY + n[j][Z] * tauYZ)) * area[j];
            force[Z] +=
                (pw * n[j][Z] - (n[j][X] * tauXZ + n[j][Y] * tauYZ + n[j][Z] * tauZZ)) * area[j];
        }
    }
    for (long d = 0; d < N_DIMS; ++d) force[d] = sync_sum(force[d]) / volume;
}
