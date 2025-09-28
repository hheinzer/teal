#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "mesh.h"
#include "teal/arena.h"
#include "teal/utils.h"
#include "teal/vector.h"

static long find_entity(const MeshEntities *entities, const char *entity)
{
    for (long i = 0; i < entities->num; i++) {
        if (!strcmp(entities->name[i], entity)) {
            return i;
        }
    }
    return -1;
}

/* Grow entities arrays by one slot. */
static void grow_entities(MeshEntities *entities, long idx)
{
    if (idx < entities->num_inner) {
        entities->num_inner += 1;
    }
    else if (idx < entities->num_inner + entities->num_ghost) {
        entities->num_ghost += 1;
    }
    entities->num += 1;

    string *name = arena_calloc(entities->num, sizeof(*name));
    entities->name = memcpy(name, entities->name, (entities->num - 1) * sizeof(*name));

    long *cell_off = arena_malloc(entities->num + 1, sizeof(*cell_off));
    entities->cell_off = memcpy(cell_off, entities->cell_off, entities->num * sizeof(*cell_off));

    vector *offset = arena_malloc(entities->num, sizeof(*offset));
    entities->offset = memcpy(offset, entities->offset, (entities->num - 1) * sizeof(*offset));
}

/* Build a map for cells [beg,end) that encodes which side of the plane a cell center lies on. */
static long *compute_cell_map(const MeshNodes *nodes, const MeshCells *cells, vector root,
                              vector normal, long beg, long end)
{
    long *map = arena_malloc(end - beg, sizeof(*map));
    for (long num = 0, i = beg; i < end; i++, num++) {
        vector center = {0};
        long num_nodes = cells->node.off[i + 1] - cells->node.off[i];
        for (long j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
            vector coord = nodes->coord[cells->node.idx[j]];
            center = vector_add(center, vector_div(coord, num_nodes));
        }
        map[num] = (vector_dot(vector_sub(center, root), normal) <= 0) ? i : end + i;
    }
    return map;
}

/* Reorder cell-to-node connectivity for [beg,end) according to `map`. */
static void reorder_cells(MeshCells *cells, long beg, long end, const long *map)
{
    Arena save = arena_save();

    typedef struct {
        long map;
        long num;
        long node[MAX_CELL_NODES];
    } Cell;

    Cell *cell = arena_malloc(end - beg, sizeof(*cell));

    for (long num = 0, i = beg; i < end; i++, num++) {
        cell[num].map = map[num];
        cell[num].num = cells->node.off[i + 1] - cells->node.off[i];
        for (long k = 0, j = cells->node.off[i]; j < cells->node.off[i + 1]; j++, k++) {
            cell[num].node[k] = cells->node.idx[j];
        }
    }
    qsort(cell, end - beg, sizeof(*cell), lcmp);

    for (long num = 0, i = beg; i < end; i++, num++) {
        cells->node.off[i + 1] = cells->node.off[i] + cell[num].num;
        for (long k = 0, j = cells->node.off[i]; j < cells->node.off[i + 1]; j++, k++) {
            cells->node.idx[j] = cell[num].node[k];
        }
    }

    arena_load(save);
}

/* Split entity metadata at `idx` into two consecutive entities. */
static void split_entities(MeshEntities *entities, long idx, long beg, long end, const long *map)
{
    for (long i = entities->num - 1; i > idx; i--) {
        strcpy(entities->name[i], entities->name[i - 1]);
        entities->cell_off[i + 1] = entities->cell_off[i];
        entities->offset[i] = entities->offset[i - 1];
    }

    strcat(entities->name[idx], "-a");
    strcat(entities->name[idx + 1], "-b");

    for (long num = 0, i = beg; i < end; i++, num++) {
        entities->cell_off[idx + 1] -= map[num] >= end;
    }

    entities->offset[idx + 1] = entities->offset[idx];
}

void mesh_split(Mesh *mesh, const char *entity, vector root, vector normal)
{
    assert(mesh);

    long idx = find_entity(&mesh->entities, entity);
    assert(idx != -1);

    grow_entities(&mesh->entities, idx);

    Arena save = arena_save();

    long beg = mesh->entities.cell_off[idx];
    long end = mesh->entities.cell_off[idx + 1];
    long *map = compute_cell_map(&mesh->nodes, &mesh->cells, root, normal, beg, end);

    reorder_cells(&mesh->cells, beg, end, map);
    split_entities(&mesh->entities, idx, beg, end, map);

    arena_load(save);
}
