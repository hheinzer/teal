#include <stddef.h>
#include <string.h>

#include "mesh.h"
#include "reorder.h"
#include "teal/arena.h"
#include "teal/assert.h"
#include "teal/vector.h"

static number find_entity(const MeshEntities *entities, const char *entity)
{
    for (number i = 0; i < entities->num; i++) {
        if (!strcmp(entities->name[i], entity)) {
            return i;
        }
    }
    return -1;
}

/* Grow entities arrays by one slot. */
static void grow_entities(MeshEntities *entities, number idx)
{
    if (idx < entities->num_inner) {
        entities->num_inner += 1;
    }
    else if (idx < entities->off_ghost) {
        entities->off_ghost += 1;
    }
    entities->num += 1;

    Name *name = arena_calloc(entities->num, sizeof(*name));
    entities->name = memcpy(name, entities->name, (entities->num - 1) * sizeof(*name));

    number *cell_off = arena_malloc(entities->num + 1, sizeof(*cell_off));
    entities->cell_off = memcpy(cell_off, entities->cell_off, entities->num * sizeof(*cell_off));
    entities->cell_off[entities->num] = entities->cell_off[entities->num - 1];

    vector *translation = arena_malloc(entities->num, sizeof(*translation));
    entities->translation =
        memcpy(translation, entities->translation, (entities->num - 1) * sizeof(*translation));
}

/* Reorder cells based on which side of the splitting plane they are on. */
static number reorder(const MeshNodes *nodes, MeshCells *cells, const MeshEntities *entities,
                      vector root, vector normal, number idx)
{
    Arena save = arena_save();

    number beg = entities->cell_off[idx];
    number end = entities->cell_off[idx + 1];
    number tot = end - beg;
    number *map = arena_malloc(tot, sizeof(*map));

    number off = 0;
    for (number num = 0, i = beg; i < end; i++, num++) {
        vector center = {0};
        number num_nodes = cells->node.off[i + 1] - cells->node.off[i];
        for (number j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
            vector coord = nodes->coord[cells->node.idx[j]];
            vector_inc(&center, vector_div(coord, num_nodes));
        }
        map[num] = vector_dot(vector_sub(center, root), normal) <= 0;
        off += map[num];
    }

    number lhs = 0;
    number rhs = 0;
    for (number i = 0; i < tot; i++) {
        map[i] = map[i] ? lhs++ : (off + rhs++);
    }
    assert(lhs == off);
    assert(rhs == tot - off);

    mesh_reorder_cells(cells, 0, beg, end, map);

    arena_load(save);
    return beg + off;
}

/* Split entity metadata at `idx` into two consecutive entities. */
static void split_entities(MeshEntities *entities, number idx, number cell_off)
{
    for (number i = entities->num - 1; i > idx; i--) {
        strcpy(entities->name[i], entities->name[i - 1]);
        entities->cell_off[i + 1] = entities->cell_off[i];
        entities->translation[i] = entities->translation[i - 1];
    }
    strcat(entities->name[idx], "-a");
    strcat(entities->name[idx + 1], "-b");
    entities->cell_off[idx + 1] = cell_off;
    entities->translation[idx + 1] = entities->translation[idx];
}

void mesh_split(Mesh *mesh, const char *entity, vector root, vector normal)
{
    assert(mesh);

    number idx = find_entity(&mesh->entities, entity);
    assert(idx != -1);

    grow_entities(&mesh->entities, idx);

    number cell_off = reorder(&mesh->nodes, &mesh->cells, &mesh->entities, root, normal, idx);

    split_entities(&mesh->entities, idx, cell_off);
}
