#include "restart.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "teal/h5io.h"
#include "teal/isclose.h"
#include "teal/kdtree.h"
#include "teal/memory.h"
#include "teal/sync.h"
#include "teal/utils.h"

static long read_restart(Simulation *sim, const char *fname, Vector3d **x, double **u);

void simulation_restart(Simulation *sim, const char *fname)
{
    if (!fname) return;

    const char *underscore = strrchr(fname, '_');
    if (underscore) sim->output_count = max(strtol(underscore + 1, 0, 10), 0);

    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    const long n_vars = sim->eqns->n_vars;
    const alias(x, sim->eqns->mesh->cell.center);
    double(*u)[n_vars] = (void *)sim->eqns->vars.u;

    const long n_global_cells = sync_sum(n_inner_cells);
    smart Vector3d *gx = (sync.rank == 0 ? memory_calloc(n_global_cells, sizeof(*gx)) : 0);
    sync_gatherv(*gx, *x, n_inner_cells * N_DIMS, 0);

    sim->time = 0;
    smart double(*gu)[n_vars] = (sync.rank == 0 ? memory_calloc(n_global_cells, sizeof(*gu)) : 0);
    if (sync.rank == 0) {
        smart Vector3d *rx;
        smart double *ru;
        const long n_restart_cells = read_restart(sim, fname, &rx, &ru);

        defer(kdtree_free) Kdtree *x2i = kdtree_create(n_restart_cells, N_DIMS);
        for (long i = 0; i < n_restart_cells; ++i) kdtree_insert(x2i, rx[i], &i, 1);

        for (long i = 0; i < n_global_cells; ++i) {
            const KdtreeItem *item = kdtree_nearest(x2i, gx[i], 1);
            assert(item && is_close(item->dist2, 0));
            for (long v = 0; v < n_vars; ++v) gu[i][v] = ru[*item->val * n_vars + v];
        }
    }

    sim->time = sync_sum(sim->time);
    sync_scatterv(*u, *gu, n_inner_cells * n_vars, 0);
}

static long read_restart(Simulation *sim, const char *fname, Vector3d **x, double **u)
{
    hid_t file = h5_file_open(fname);
    hid_t vtkhdf = h5_group_open("VTKHDF", file);

    long n_inner_cells;
    h5_dataset_read("NumberOfCells", &n_inner_cells, vtkhdf);

    hid_t field = h5_group_open("FieldData", vtkhdf);
    h5_dataset_read("time", &sim->time, field);
    h5_group_close(field);

    hid_t cell = h5_group_open("CellData", vtkhdf);
    *x = memory_calloc(n_inner_cells, sizeof(**x));
    h5_dataset_read("center", *x, cell);

    const long n_vars = sim->eqns->n_vars;
    const alias(vdim, sim->eqns->vars.dim);
    const alias(vname, sim->eqns->vars.name);
    smart double *buf = memory_calloc(n_inner_cells * N_DIMS, sizeof(*buf));
    *u = memory_calloc(n_inner_cells * n_vars, sizeof(**u));
    for (long n = 0, v = 0; v < n_vars; v += vdim[n], ++n) {
        h5_dataset_read(vname[v], buf, cell);
        for (long i = 0; i < n_inner_cells; ++i)
            for (long d = 0; d < vdim[n]; ++d) (*u)[i * n_vars + v + d] = buf[i * vdim[n] + d];
    }

    h5_group_close(cell);
    h5_group_close(vtkhdf);
    h5_file_close(file);

    return n_inner_cells;
}
