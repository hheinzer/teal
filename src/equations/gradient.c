#include "equations.h"
#include "sync.h"
#include "teal/memory.h"
#include "teal/utils.h"

void equations_gradient(Equations *eqns)
{
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const long n_inner_faces = eqns->mesh->n_inner_faces;
    const long n_ghost_faces = eqns->mesh->n_ghost_faces;
    const long n_faces = eqns->mesh->n_faces;
    const long n_vars = eqns->n_vars;
    const alias(cell, eqns->mesh->face.cell);
    const alias(weight, eqns->mesh->face.weight);
    double(*u)[n_vars] = (void *)eqns->vars.u;
    Vector3d(*dudx)[n_vars] = (void *)eqns->vars.dudx;

    sync_begin(eqns, *u, n_vars);

    memory_setzero(dudx, n_inner_cells, sizeof(*dudx));

    for (long i = 0; i < n_inner_faces; ++i) {
        for (long v = 0; v < n_vars; ++v) {
            const double du = u[cell[i][R]][v] - u[cell[i][L]][v];
            for (long d = 0; d < N_DIMS; ++d) {
                dudx[cell[i][L]][v][d] += weight[i][d] * du;
                dudx[cell[i][R]][v][d] += weight[i][d] * du;
            }
        }
    }
    for (long i = n_inner_faces; i < n_inner_faces + n_ghost_faces; ++i) {
        for (long v = 0; v < n_vars; ++v) {
            const double du = u[cell[i][R]][v] - u[cell[i][L]][v];
            for (long d = 0; d < N_DIMS; ++d) {
                dudx[cell[i][L]][v][d] += weight[i][d] * du;
            }
        }
    }

    sync_wait(eqns);

    for (long i = n_inner_faces + n_ghost_faces; i < n_faces; ++i) {
        for (long v = 0; v < n_vars; ++v) {
            const double du = u[cell[i][R]][v] - u[cell[i][L]][v];
            for (long d = 0; d < N_DIMS; ++d) {
                dudx[cell[i][L]][v][d] += weight[i][d] * du;
                dudx[cell[i][R]][v][d] += weight[i][d] * du;
            }
        }
    }

    sync_end(eqns);
}
