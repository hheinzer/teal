#include <string.h>

#include "find.h"
#include "mesh.h"
#include "reorder.h"
#include "teal/array.h"
#include "teal/memory.h"
#include "teal/sync.h"
#include "teal/utils.h"

static void compute_split(const Mesh *mesh, const Vector3d root, const Vector3d normal, long *split,
                          long E);

static void compute_reorder(const Mesh *mesh, const long *split, long *new2old, long *old2new,
                            long E);

static void split_entities(Mesh *mesh, const long *split, long E, long n_cells);

void mesh_split(Mesh *mesh, const char *entity, const Vector3d root, const Vector3d normal)
{
    if (sync.rank != 0) return;

    const long E = mesh_find_entity(mesh, entity);
    const long n_cells = mesh->entity.j_cell[E + 1] - mesh->entity.j_cell[E];

    smart long *split = memory_calloc(n_cells, sizeof(*split));
    compute_split(mesh, root, normal, split, E);

    smart long *new2old = memory_calloc(mesh->n_cells, sizeof(*new2old));
    smart long *old2new = memory_calloc(mesh->n_cells, sizeof(*old2new));
    compute_reorder(mesh, split, new2old, old2new, E);
    reorder_cells(mesh, old2new, new2old);

    split_entities(mesh, split, E, n_cells);
}

static void compute_split(const Mesh *mesh, const Vector3d root, const Vector3d normal, long *split,
                          long E)
{
    const alias(coord, mesh->node.coord);
    const alias(i_node, mesh->cell.i_node);
    const alias(node, mesh->cell.node);
    const alias(j_cell, mesh->entity.j_cell);
    for (long s = 0, j = j_cell[E]; j < j_cell[E + 1]; ++j) {
        const long n_nodes = i_node[j + 1] - i_node[j];

        Vector3d center = {0};
        for (long i = i_node[j]; i < i_node[j + 1]; ++i)
            for (long d = 0; d < N_DIMS; ++d) center[d] += coord[node[i]][d] / n_nodes;

        Vector3d rc;
        for (long d = 0; d < N_DIMS; ++d) rc[d] = center[d] - root[d];

        split[s++] = (array_dot(rc, normal, N_DIMS) > 0 ? 0 : 1);
    }
}

static void compute_reorder(const Mesh *mesh, const long *split, long *new2old, long *old2new,
                            long E)
{
    const alias(j_cell, mesh->entity.j_cell);
    for (long n = 0, e = 0; e < mesh->n_entities; ++e) {
        if (e != E)
            for (long j = j_cell[e]; j < j_cell[e + 1]; ++j) new2old[n++] = j;
        else
            for (long k = 0; k < 2; ++k)
                for (long s = 0, j = j_cell[e]; j < j_cell[e + 1]; ++j)
                    if (split[s++] == k) new2old[n++] = j;
    }
    for (long i = 0; i < mesh->n_cells; ++i) old2new[new2old[i]] = i;
}

static void split_entities(Mesh *mesh, const long *split, long E, long n_cells)
{
    mesh->n_entities += 1;
    String *name = memory_realloc(mesh->entity.name, mesh->n_entities, sizeof(*name));
    long *j_cell = memory_realloc(mesh->entity.j_cell, mesh->n_entities + 1, sizeof(*j_cell));
    Vector3d *offset = memory_realloc(mesh->entity.offset, mesh->n_entities, sizeof(*offset));
    for (long e = mesh->n_entities - 1; e > E; --e) {
        strcpy(name[e], name[e - 1]);
        j_cell[e + 1] = j_cell[e];
        for (long d = 0; d < N_DIMS; ++d) offset[e][d] = offset[e - 1][d];
    }
    strcpy(name[E + 1], name[E]);
    for (long d = 0; d < N_DIMS; ++d) offset[E + 1][d] = offset[E][d];
    strcat(name[E + 0], "-a");
    strcat(name[E + 1], "-b");
    for (long s = 0; s < n_cells; ++s) j_cell[E + 1] -= split[s];
    mesh->entity.name = name;
    mesh->entity.j_cell = j_cell;
    mesh->entity.offset = offset;
}
