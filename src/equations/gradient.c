#include "core/memory.h"
#include "core/sync.h"
#include "core/utils.h"
#include "equations.h"

void equations_gradient(Equations *eqns)
{
    const long n_vars = eqns->n_vars;
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const long n_inner_faces = eqns->mesh->n_inner_faces;
    const long n_bound_faces = eqns->mesh->n_bound_faces;
    const long n_faces = eqns->mesh->n_faces;
    const ALIAS(cell, eqns->mesh->face.cell);
    const ALIAS(weight, eqns->mesh->face.gradient_weight);
    double(*u)[n_vars] = (void *)eqns->vars.u;
    double(*dudx)[n_vars][N_DIMS] = (void *)eqns->vars.dudx;

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
    for (long i = n_inner_faces; i < n_inner_faces + n_bound_faces; ++i) {
        for (long v = 0; v < n_vars; ++v) {
            const double du = u[cell[i][R]][v] - u[cell[i][L]][v];
            for (long d = 0; d < N_DIMS; ++d) {
                dudx[cell[i][L]][v][d] += weight[i][d] * du;
            }
        }
    }

    sync_wait(eqns);

    for (long i = n_inner_faces + n_bound_faces; i < n_faces; ++i) {
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
