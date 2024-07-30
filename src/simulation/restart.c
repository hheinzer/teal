#include <mpi.h>
#include <stdlib.h>
#include <string.h>

#include "core/array.h"
#include "core/h5io.h"
#include "core/kdtree.h"
#include "core/memory.h"
#include "core/utils.h"
#include "simulation.h"
#include "teal.h"

static long read_restart(Simulation *sim, const char *fname, double (**center)[N_DIMS],
                         double **vars);

static void compute_variables(const Kdtree *x2i, const double *vars, const double (*x)[N_DIMS],
                              double *u, long n_inner_cells, long n_vars);

void simulation_restart(Simulation *sim, const char *fname)
{
    const long n_vars = sim->eqns->n_vars;
    const long n_local_cells = sim->eqns->mesh->n_inner_cells;
    const long nlx = n_local_cells * N_DIMS;
    const long nlu = n_local_cells * n_vars;
    const ALIAS(lx, sim->eqns->mesh->cell.center);
    ALIAS(lu, sim->eqns->vars.u);

    const char *underscore = strrchr(fname, '_');
    if (underscore) sim->output_count = max(strtol(underscore + 1, 0, 10), 0);

    if (teal.rank == 0) {
        smart long *n_inner_cells = memory_calloc(teal.size, sizeof(*n_inner_cells));
        MPI_Gather(&n_local_cells, 1, MPI_LONG, n_inner_cells, 1, MPI_LONG, 0, teal.comm);
        const long n_global_cells = array_sum(n_inner_cells, teal.size);

        smart double(*gx)[N_DIMS] = memory_calloc(n_global_cells, sizeof(*gx));
        smart int *counts = memory_calloc(teal.size, sizeof(*counts));
        smart int *displs = memory_calloc(teal.size, sizeof(*displs));
        for (long rank = 0; rank < teal.size; ++rank) {
            counts[rank] = n_inner_cells[rank] * N_DIMS;
            if (rank > 0) displs[rank] = displs[rank - 1] + counts[rank - 1];
        }
        MPI_Gatherv(lx, nlx, MPI_DOUBLE, gx, counts, displs, MPI_DOUBLE, 0, teal.comm);

        smart double(*center)[N_DIMS], *vars;
        const long n_restart_cells = read_restart(sim, fname, &center, &vars);

        defer(kdtree_free) Kdtree x2i = kdtree_create(n_restart_cells, N_DIMS);
        for (long i = 0; i < n_restart_cells; ++i) kdtree_insert(&x2i, center[i], &i, 1);

        double *u = memory_calloc(n_global_cells * n_vars, sizeof(*u));
        compute_variables(&x2i, vars, gx, u, n_global_cells, n_vars);

        for (long rank = 0; rank < teal.size; ++rank) {
            counts[rank] = n_inner_cells[rank] * n_vars;
            if (rank > 0) displs[rank] = displs[rank - 1] + counts[rank - 1];
        }
        MPI_Scatterv(u, counts, displs, MPI_DOUBLE, lu, nlu, MPI_DOUBLE, 0, teal.comm);
    }
    else {
        MPI_Gather(&n_local_cells, 1, MPI_LONG, 0, 0, MPI_LONG, 0, teal.comm);
        MPI_Gatherv(lx, nlx, MPI_DOUBLE, 0, 0, 0, MPI_DOUBLE, 0, teal.comm);
        MPI_Scatterv(0, 0, 0, MPI_DOUBLE, lu, nlu, MPI_DOUBLE, 0, teal.comm);
    }
}

static long read_restart(Simulation *sim, const char *fname, double (**center)[N_DIMS],
                         double **vars)
{
    hid_t file = h5_file_open(fname);
    hid_t vtkhdf = h5_group_open("VTKHDF", file);

    long n_inner_cells;
    h5_dataset_read("NumberOfCells", &n_inner_cells, vtkhdf);

    hid_t field = h5_group_open("FieldData", vtkhdf);
    h5_dataset_read("time", &sim->time, field);
    h5_dataset_read("iter", &sim->iter, field);
    h5_group_close(field);

    hid_t cell = h5_group_open("CellData", vtkhdf);
    *center = memory_calloc(n_inner_cells, sizeof(**center));
    h5_dataset_read("center", *center, cell);

    const long n_vars = sim->eqns->n_vars;
    const ALIAS(vdim, sim->eqns->vars.dim);
    const ALIAS(vname, sim->eqns->vars.name);
    smart double *buf = memory_calloc(n_inner_cells * N_DIMS, sizeof(*buf));
    *vars = memory_calloc(n_inner_cells * n_vars, sizeof(**vars));
    for (long n = 0, v = 0; v < n_vars; v += vdim[n], ++n) {
        h5_dataset_read(vname[v], buf, cell);
        for (long i = 0; i < n_inner_cells; ++i)
            for (long d = 0; d < vdim[n]; ++d) (*vars)[i * n_vars + v + d] = buf[i * vdim[n] + d];
    }

    h5_group_close(cell);
    h5_group_close(vtkhdf);
    h5_file_close(file);

    return n_inner_cells;
}

static void compute_variables(const Kdtree *x2i, const double *vars, const double (*x)[N_DIMS],
                              double *u, long n_inner_cells, long n_vars)
{
    for (long i = 0; i < n_inner_cells; ++i) {
        const KdtreeItem *item = kdtree_nearest(x2i, x[i], 1);
        if (!item) error("could not find nearest neighbor of cell '%ld'", i);
        for (long v = 0; v < n_vars; ++v) u[i * n_vars + v] = vars[*item->val * n_vars + v];
    }
}
