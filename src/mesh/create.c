#include <string.h>

#include "connectivity.h"
#include "core/array.h"
#include "core/memory.h"
#include "core/utils.h"
#include "mesh.h"
#include "teal.h"

static void create_nodes(Mesh *mesh, const double *x0, const double *x1, const long *n_cells,
                         const long *n_nodes);

static void create_cells(Mesh *mesh, const long *n_cells, const long *n_nodes);

static void create_entities(Mesh *mesh, const long *n_cells);

Mesh mesh_create(const double *x0, const double *x1, const long *n_cells)
{
    for (long d = 0; d < N_DIMS; ++d) {
        if (x0[d] >= x1[d]) error("illegal domain extent '%g' to '%g'", x0[d], x1[d]);
        if (n_cells[d] <= 0) error("illegal number of cells '%ld'", n_cells[d]);
    }

    Mesh mesh = {0};
    if (teal.rank != 0) return mesh;

    long n_nodes[N_DIMS];
    for (long d = 0; d < N_DIMS; ++d) n_nodes[d] = n_cells[d] + 1;

    create_nodes(&mesh, x0, x1, n_cells, n_nodes);
    create_cells(&mesh, n_cells, n_nodes);
    create_entities(&mesh, n_cells);

    connectivity_cells(&mesh);

    return mesh;
}

static void create_nodes(Mesh *mesh, const double *x0, const double *x1, const long *n_cells,
                         const long *n_nodes)
{
    mesh->n_inner_nodes = array_product(n_nodes, N_DIMS);
    mesh->n_nodes = mesh->n_inner_nodes;

    double(*coord)[N_DIMS] = memory_calloc(mesh->n_nodes, sizeof(*coord));
    for (long n = 0, k = 0; k < n_nodes[Z]; ++k) {
        for (long j = 0; j < n_nodes[Y]; ++j) {
            for (long i = 0; i < n_nodes[X]; ++i) {
                coord[n][X] = x0[X] + (x1[X] - x0[X]) * i / n_cells[X];
                coord[n][Y] = x0[Y] + (x1[Y] - x0[Y]) * j / n_cells[Y];
                coord[n][Z] = x0[Z] + (x1[Z] - x0[Z]) * k / n_cells[Z];
                n += 1;
            }
        }
    }
    mesh->node.coord = coord;
}

static void create_cells(Mesh *mesh, const long *n_cells, const long *n_nodes)
{
    mesh->n_inner_cells = array_product(n_cells, N_DIMS);
    for (long d0 = 0; d0 < N_DIMS; ++d0)
        for (long d1 = 0; d1 < N_DIMS; ++d1)
            if (d1 != d0) mesh->n_ghost_cells += n_cells[d0] * n_cells[d1];
    mesh->n_cells = mesh->n_inner_cells + mesh->n_ghost_cells;

    long *i_node = memory_calloc(mesh->n_cells + 1, sizeof(*i_node));
    long *node = memory_calloc(mesh->n_cells * MAX_CELL_NODES, sizeof(*node));
    long n = 0;
    for (long k = 0; k < n_cells[Z]; ++k) {
        for (long j = 0; j < n_cells[Y]; ++j) {
            for (long i = 0; i < n_cells[X]; ++i) {
                i_node[n + 1] = i_node[n];
                node[i_node[n + 1]++] = (i + 0) + n_nodes[X] * ((j + 0) + n_nodes[Y] * (k + 0));
                node[i_node[n + 1]++] = (i + 1) + n_nodes[X] * ((j + 0) + n_nodes[Y] * (k + 0));
                node[i_node[n + 1]++] = (i + 1) + n_nodes[X] * ((j + 1) + n_nodes[Y] * (k + 0));
                node[i_node[n + 1]++] = (i + 0) + n_nodes[X] * ((j + 1) + n_nodes[Y] * (k + 0));
                node[i_node[n + 1]++] = (i + 0) + n_nodes[X] * ((j + 0) + n_nodes[Y] * (k + 1));
                node[i_node[n + 1]++] = (i + 1) + n_nodes[X] * ((j + 0) + n_nodes[Y] * (k + 1));
                node[i_node[n + 1]++] = (i + 1) + n_nodes[X] * ((j + 1) + n_nodes[Y] * (k + 1));
                node[i_node[n + 1]++] = (i + 0) + n_nodes[X] * ((j + 1) + n_nodes[Y] * (k + 1));
                n += 1;
            }
        }
    }
    for (long i = 0; i <= n_cells[X]; i += n_cells[X]) {
        for (long k = 0; k < n_cells[Z]; ++k) {
            for (long j = 0; j < n_cells[Y]; ++j) {
                i_node[n + 1] = i_node[n];
                node[i_node[n + 1]++] = (i + 0) + n_nodes[X] * ((j + 0) + n_nodes[Y] * (k + 0));
                node[i_node[n + 1]++] = (i + 0) + n_nodes[X] * ((j + 1) + n_nodes[Y] * (k + 0));
                node[i_node[n + 1]++] = (i + 0) + n_nodes[X] * ((j + 0) + n_nodes[Y] * (k + 1));
                node[i_node[n + 1]++] = (i + 0) + n_nodes[X] * ((j + 1) + n_nodes[Y] * (k + 1));
                n += 1;
            }
        }
    }
    for (long j = 0; j <= n_cells[Y]; j += n_cells[Y]) {
        for (long k = 0; k < n_cells[Z]; ++k) {
            for (long i = 0; i < n_cells[X]; ++i) {
                i_node[n + 1] = i_node[n];
                node[i_node[n + 1]++] = (i + 0) + n_nodes[X] * ((j + 0) + n_nodes[Y] * (k + 0));
                node[i_node[n + 1]++] = (i + 1) + n_nodes[X] * ((j + 0) + n_nodes[Y] * (k + 0));
                node[i_node[n + 1]++] = (i + 0) + n_nodes[X] * ((j + 0) + n_nodes[Y] * (k + 1));
                node[i_node[n + 1]++] = (i + 1) + n_nodes[X] * ((j + 0) + n_nodes[Y] * (k + 1));
                n += 1;
            }
        }
    }
    for (long k = 0; k <= n_cells[Z]; k += n_cells[Z]) {
        for (long j = 0; j < n_cells[Y]; ++j) {
            for (long i = 0; i < n_cells[X]; ++i) {
                i_node[n + 1] = i_node[n];
                node[i_node[n + 1]++] = (i + 0) + n_nodes[X] * ((j + 0) + n_nodes[Y] * (k + 0));
                node[i_node[n + 1]++] = (i + 1) + n_nodes[X] * ((j + 0) + n_nodes[Y] * (k + 0));
                node[i_node[n + 1]++] = (i + 1) + n_nodes[X] * ((j + 1) + n_nodes[Y] * (k + 0));
                node[i_node[n + 1]++] = (i + 0) + n_nodes[X] * ((j + 1) + n_nodes[Y] * (k + 0));
                n += 1;
            }
        }
    }
    mesh->cell.i_node = i_node;
    mesh->cell.node = memory_realloc(node, i_node[mesh->n_cells], sizeof(*node));
}

static void create_entities(Mesh *mesh, const long *n_cells)
{
    mesh->n_entities = 7;

    mesh->entity.name = memory_calloc(mesh->n_entities, sizeof(*mesh->entity.name));
    strcpy(mesh->entity.name[0], "domain");
    strcpy(mesh->entity.name[1], "left");
    strcpy(mesh->entity.name[2], "right");
    strcpy(mesh->entity.name[3], "bottom");
    strcpy(mesh->entity.name[4], "top");
    strcpy(mesh->entity.name[5], "back");
    strcpy(mesh->entity.name[6], "front");

    long *j_cell = memory_calloc(mesh->n_entities + 1, sizeof(*j_cell));
    j_cell[1] = mesh->n_inner_cells;
    j_cell[2] = j_cell[1] + n_cells[Z] * n_cells[Y];
    j_cell[3] = j_cell[2] + n_cells[Z] * n_cells[Y];
    j_cell[4] = j_cell[3] + n_cells[Z] * n_cells[X];
    j_cell[5] = j_cell[4] + n_cells[Z] * n_cells[X];
    j_cell[6] = j_cell[5] + n_cells[Y] * n_cells[X];
    j_cell[7] = j_cell[6] + n_cells[Y] * n_cells[X];
    mesh->entity.j_cell = j_cell;

    mesh->entity.offset = memory_calloc(mesh->n_entities, sizeof(*mesh->entity.offset));
}
