#include <stdlib.h>
#include <string.h>

#include "core/kdtree.h"
#include "core/memory.h"
#include "core/utils.h"
#include "mesh.h"
#include "teal.h"

static void find_entities(const Mesh *mesh, const char *entity0, const char *entity1, long *E);
static void compute_cell_coordinates(const Mesh *mesh, const long *E, long n_cells,
                                     double (*mean)[N_DIMS], double (*x)[n_cells][N_DIMS]);
static void compute_coordinate_to_cell(const Mesh *mesh, const long *E, long n_cells,
                                       const double (*x)[n_cells][N_DIMS], Kdtree *x2c);
static void compute_cell_to_cell(const Mesh *mesh, const long *E, const double (*mean)[N_DIMS],
                                 const Kdtree *x2c, long n_cells,
                                 const double (*x)[n_cells][N_DIMS], long *c2c);
static void make_cells_periodic(Mesh *mesh, const long *E, const long *c2c);
static void make_entities_periodic(Mesh *mesh, const long *E, const double (*mean)[N_DIMS]);

void mesh_periodic(Mesh *mesh, const char *entity0, const char *entity1)
{
    if (teal.rank != 0) return;

    long E[2];
    find_entities(mesh, entity0, entity1, E);
    const long n_cells = mesh->entity.j_cell[E[0] + 1] - mesh->entity.j_cell[E[0]];
    if (n_cells != mesh->entity.j_cell[E[1] + 1] - mesh->entity.j_cell[E[1]])
        error("entities '%s' and '%s' do not have the same number of cells", entity0, entity1);

    double mean[2][N_DIMS] = {0};
    cleanup double(*x)[n_cells][N_DIMS] = memory_calloc(2, sizeof(*x));
    compute_cell_coordinates(mesh, E, n_cells, mean, x);

    fcleanup(kdtree_free) Kdtree x2c = kdtree_create(2 * n_cells, N_DIMS);
    cleanup long *c2c = memory_calloc(mesh->n_cells, sizeof(*c2c));
    compute_coordinate_to_cell(mesh, E, n_cells, x, &x2c);
    compute_cell_to_cell(mesh, E, mean, &x2c, n_cells, x, c2c);

    make_cells_periodic(mesh, E, c2c);
    make_entities_periodic(mesh, E, mean);
}

static void find_entities(const Mesh *mesh, const char *entity0, const char *entity1, long *E)
{
    for (long e = 0; e < 2; ++e) E[e] = -1;
    for (long e = 0; e < mesh->n_entities; ++e) {
        if (!strcmp(mesh->entity.name[e], entity0)) E[0] = e;
        if (!strcmp(mesh->entity.name[e], entity1)) E[1] = e;
    }
    if (E[0] == -1) error("could not find entity '%s'", entity0);
    if (E[1] == -1) error("could not find entity '%s'", entity1);
}

static void compute_cell_coordinates(const Mesh *mesh, const long *E, long n_cells,
                                     double (*mean)[N_DIMS], double (*x)[n_cells][N_DIMS])
{
    const ALIAS(j_cell, mesh->entity.j_cell);
    const ALIAS(i_node, mesh->cell.i_node);
    const ALIAS(node, mesh->cell.node);
    const ALIAS(coord, mesh->node.coord);
    for (long e = 0; e < 2; ++e) {
        for (long n = 0, j = j_cell[E[e]]; j < j_cell[E[e] + 1]; ++j) {
            const long n_nodes = i_node[j + 1] - i_node[j];
            for (long i = i_node[j]; i < i_node[j + 1]; ++i)
                for (long d = 0; d < N_DIMS; ++d) x[e][n][d] += coord[node[i]][d] / n_nodes;
            for (long d = 0; d < N_DIMS; ++d) mean[e][d] += x[e][n][d] / n_cells;
            n += 1;
        }
    }
}

static void compute_coordinate_to_cell(const Mesh *mesh, const long *E, long n_cells,
                                       const double (*x)[n_cells][N_DIMS], Kdtree *x2c)
{
    const ALIAS(j_cell, mesh->entity.j_cell);
    for (long e = 0; e < 2; ++e)
        for (long n = 0, j = j_cell[E[e]]; j < j_cell[E[e] + 1]; ++j)
            kdtree_insert(x2c, x[e][n++], &j, 1);
}

static void compute_cell_to_cell(const Mesh *mesh, const long *E, const double (*mean)[N_DIMS],
                                 const Kdtree *x2c, long n_cells,
                                 const double (*x)[n_cells][N_DIMS], long *c2c)
{
    const ALIAS(j_cell, mesh->entity.j_cell);
    const ALIAS(i_cell, mesh->cell.i_cell);
    const ALIAS(cell, mesh->cell.cell);
    for (long j = 0; j < mesh->n_cells; ++j) c2c[j] = j;
    for (long e = 0; e < 2; ++e) {
        for (long n = 0, j0 = j_cell[E[e]]; j0 < j_cell[E[e] + 1]; ++j0) {
            double key[N_DIMS];
            for (long d = 0; d < N_DIMS; ++d)
                key[d] = x[e][n][d] - mean[e][d] + mean[(e + 1) % 2][d];
            cleanup KdtreeItem *item = kdtree_nearest(x2c, key, 1);
            if (!item || !isclose(item->dist2, 0))
                error("could not find periodic connection of cell '%ld'", j0);
            const long ghost = j0;
            const long inner = cell[i_cell[*item->val]];
            c2c[ghost] = inner;
            n += 1;
        }
    }
}

static void make_cells_periodic(Mesh *mesh, const long *E, const long *c2c)
{
    const ALIAS(j_cell, mesh->entity.j_cell);
    long *i_cell = memory_calloc(mesh->n_cells + 1, sizeof(*i_cell));
    long *cell = memory_calloc(mesh->n_cells * MAX_CELL_FACES, sizeof(*cell));
    for (long n = 0, e = 0; e < mesh->n_entities; ++e) {
        for (long j0 = j_cell[e]; j0 < j_cell[e + 1]; ++j0) {
            i_cell[n + 1] = i_cell[n];
            for (long i = mesh->cell.i_cell[j0]; i < mesh->cell.i_cell[j0 + 1]; ++i) {
                const long inner = c2c[mesh->cell.cell[i]];
                if (inner == j0) continue;  // filter out self connections
                cell[i_cell[n + 1]++] = inner;
            }
            if (e == E[0] || e == E[1]) cell[i_cell[n + 1]++] = c2c[j0];
            n += 1;
        }
    }
    free(mesh->cell.i_cell);
    free(mesh->cell.cell);
    mesh->cell.i_cell = i_cell;
    mesh->cell.cell = memory_realloc(cell, i_cell[mesh->n_cells], sizeof(*cell));
}

static void make_entities_periodic(Mesh *mesh, const long *E, const double (*mean)[N_DIMS])
{
    ALIAS(name, mesh->entity.name);
    for (long e = 0; e < 2; ++e) strcpy(name[E[e]], "periodic");

    ALIAS(offset, mesh->entity.offset);
    for (long d = 0; d < N_DIMS; ++d) {
        offset[E[0]][d] = mean[0][d] - mean[1][d];
        offset[E[1]][d] = mean[1][d] - mean[0][d];
    }
}
