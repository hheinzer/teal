#include "periodic.h"

#include <assert.h>

#include "find.h"
#include "teal/isclose.h"
#include "teal/kdtree.h"
#include "teal/memory.h"
#include "teal/sync.h"
#include "teal/utils.h"

static void compute_cell_coordinates(const Mesh *mesh, const Vector2l E, long n_cells,
                                     Vector3d mean[2], Vector3d x[2][n_cells]);

static void compute_coordinate_to_cell(const Mesh *mesh, const Vector2l E, long n_cells,
                                       const Vector3d x[2][n_cells], Kdtree *x2c);

static void compute_cell_to_cell(const Mesh *mesh, const Vector2l E, const Vector3d mean[2],
                                 const Kdtree *x2c, long n_cells, const Vector3d x[2][n_cells],
                                 long *c2c);

static void make_cells_periodic(Mesh *mesh, const Vector2l E, const long *c2c);

static void compute_entitiy_offsets(Mesh *mesh, const Vector2l E, const Vector3d mean[2]);

void mesh_periodic(Mesh *mesh, const char *entity_a, const char *entity_b)
{
    if (sync.rank != 0) return;

    const Vector2l E = {mesh_find_entity(mesh, entity_a), mesh_find_entity(mesh, entity_b)};
    const long n_cells = mesh->entity.j_cell[E[0] + 1] - mesh->entity.j_cell[E[0]];
    assert(n_cells == mesh->entity.j_cell[E[1] + 1] - mesh->entity.j_cell[E[1]]);

    Vector3d mean[2] = {0};
    smart Vector3d(*x)[n_cells] = memory_calloc(2, sizeof(*x));
    compute_cell_coordinates(mesh, E, n_cells, mean, x);

    defer(kdtree_free) Kdtree *x2c = kdtree_create(2 * n_cells, N_DIMS);
    smart long *c2c = memory_calloc(mesh->n_cells, sizeof(*c2c));
    compute_coordinate_to_cell(mesh, E, n_cells, x, x2c);
    compute_cell_to_cell(mesh, E, mean, x2c, n_cells, x, c2c);

    make_cells_periodic(mesh, E, c2c);
    compute_entitiy_offsets(mesh, E, mean);
}

void periodic_cell_to_node(const Mesh *mesh, Dict *periodic)
{
    const alias(i_node, mesh->cell.i_node);
    const alias(node, mesh->cell.node);
    const alias(i_cell, mesh->cell.i_cell);
    const alias(cell, mesh->cell.cell);
    for (long j = mesh->n_inner_cells; j < mesh->n_inner_cells + mesh->n_ghost_cells; ++j) {
        if (i_cell[j + 1] - i_cell[j] != N_SIDES) continue;
        const long n_nodes = i_node[j + 1] - i_node[j];
        assert(n_nodes >= N_DIMS);
        dict_insert(periodic, &cell[i_cell[j]], &node[i_node[j]], N_SIDES, n_nodes);
    }
}

void periodic_cell_to_cell(const Mesh *mesh, Dict *periodic)
{
    const alias(i_cell, mesh->cell.i_cell);
    const alias(cell, mesh->cell.cell);
    for (long j = mesh->n_inner_cells; j < mesh->n_inner_cells + mesh->n_ghost_cells; ++j) {
        if (i_cell[j + 1] - i_cell[j] != N_SIDES) continue;
        dict_insert(periodic, &cell[i_cell[j]], &j, N_SIDES, 1);
    }
}

static void compute_cell_coordinates(const Mesh *mesh, const Vector2l E, long n_cells,
                                     Vector3d mean[2], Vector3d x[2][n_cells])
{
    const alias(coord, mesh->node.coord);
    const alias(i_node, mesh->cell.i_node);
    const alias(node, mesh->cell.node);
    const alias(j_cell, mesh->entity.j_cell);
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

static void compute_coordinate_to_cell(const Mesh *mesh, const Vector2l E, long n_cells,
                                       const Vector3d x[2][n_cells], Kdtree *x2c)
{
    const alias(j_cell, mesh->entity.j_cell);
    for (long e = 0; e < 2; ++e)
        for (long n = 0, j = j_cell[E[e]]; j < j_cell[E[e] + 1]; ++j)
            kdtree_insert(x2c, x[e][n++], &j, 1);
}

static void compute_cell_to_cell(const Mesh *mesh, const Vector2l E, const Vector3d mean[2],
                                 const Kdtree *x2c, long n_cells, const Vector3d x[2][n_cells],
                                 long *c2c)
{
    const alias(i_cell, mesh->cell.i_cell);
    const alias(cell, mesh->cell.cell);
    const alias(j_cell, mesh->entity.j_cell);
    for (long j = 0; j < mesh->n_cells; ++j) c2c[j] = j;
    for (long e = 0; e < 2; ++e) {
        for (long n = 0, j0 = j_cell[E[e]]; j0 < j_cell[E[e] + 1]; ++j0) {
            Vector3d key;
            for (long d = 0; d < N_DIMS; ++d)
                key[d] = x[e][n][d] - mean[e][d] + mean[(e + 1) % 2][d];
            smart KdtreeItem *item = kdtree_nearest(x2c, key, 1);
            assert(item && is_close(item->dist2, 0));
            const long ghost = j0;
            const long inner = cell[i_cell[*item->val]];
            c2c[ghost] = inner;
            n += 1;
        }
    }
}

static void make_cells_periodic(Mesh *mesh, const Vector2l E, const long *c2c)
{
    const alias(j_cell, mesh->entity.j_cell);
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
    memory_free(&mesh->cell.i_cell);
    memory_free(&mesh->cell.cell);
    mesh->cell.i_cell = i_cell;
    mesh->cell.cell = memory_realloc(cell, i_cell[mesh->n_cells], sizeof(*cell));
}

static void compute_entitiy_offsets(Mesh *mesh, const Vector2l E, const Vector3d mean[2])
{
    alias(offset, mesh->entity.offset);
    for (long d = 0; d < N_DIMS; ++d) {
        offset[E[0]][d] = mean[0][d] - mean[1][d];
        offset[E[1]][d] = mean[1][d] - mean[0][d];
    }
}
