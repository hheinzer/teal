#include <assert.h>
#include <string.h>

#include "connectivity.h"
#include "mesh.h"
#include "teal/array.h"
#include "teal/memory.h"
#include "teal/sync.h"

static void create_nodes(Mesh *mesh, const Vector3d a, const Vector3d b, const Vector3l n_cells,
                         const Vector3l n_nodes);

static void create_cells(Mesh *mesh, const Vector3l n_cells, const Vector3l n_nodes);

static void create_entities(Mesh *mesh, const Vector3l n_cells);

Mesh *mesh_create(const Vector3d a, const Vector3d b, const Vector3l n_cells)
{
    for (long i = 0; i < N_DIMS; ++i) assert(a[i] < b[i] && n_cells[i] > 0);

    memory_sum_setzero();

    Mesh *mesh = memory_calloc(1, sizeof(*mesh));
    if (sync.rank != 0) return mesh;

    Vector3l n_nodes;
    for (long i = 0; i < N_DIMS; ++i) n_nodes[i] = n_cells[i] + 1;

    create_nodes(mesh, a, b, n_cells, n_nodes);
    create_cells(mesh, n_cells, n_nodes);
    create_entities(mesh, n_cells);

    connectivity_cells(mesh);

    return mesh;
}

static void create_nodes(Mesh *mesh, const Vector3d a, const Vector3d b, const Vector3l n_cells,
                         const Vector3l n_nodes)
{
    mesh->n_inner_nodes = array_product(n_nodes, N_DIMS);
    mesh->n_nodes = mesh->n_inner_nodes;

    Vector3d *coord = memory_calloc(mesh->n_nodes, sizeof(*coord));
    for (long n = 0, k = 0; k < n_nodes[Z]; ++k) {
        for (long j = 0; j < n_nodes[Y]; ++j) {
            for (long i = 0; i < n_nodes[X]; ++i) {
                coord[n][X] = a[X] + (b[X] - a[X]) * i / n_cells[X];
                coord[n][Y] = a[Y] + (b[Y] - a[Y]) * j / n_cells[Y];
                coord[n][Z] = a[Z] + (b[Z] - a[Z]) * k / n_cells[Z];
                n += 1;
            }
        }
    }
    mesh->node.coord = coord;
}

static void create_cells(Mesh *mesh, const Vector3l n_cells, const Vector3l n_nodes)
{
    mesh->n_inner_cells = array_product(n_cells, N_DIMS);
    for (long i = 0; i < N_DIMS; ++i)
        for (long j = 0; j < N_DIMS; ++j)
            if (i != j) mesh->n_ghost_cells += n_cells[i] * n_cells[j];
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

static void create_entities(Mesh *mesh, const Vector3l n_cells)
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
