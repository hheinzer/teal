#include "reorder.h"

#include <stdlib.h>

#include "teal/arena.h"
#include "teal/array.h"
#include "teal/assert.h"
#include "teal/utils.h"

static bool is_valid(const number *map, number num)
{
    Arena save = arena_save();

    number *count = arena_calloc(num, sizeof(*count));
    for (number i = 0; i < num; i++) {
        if (0 <= map[i] && map[i] < num) {
            count[map[i]] += 1;
        }
    }

    bool ret = true;
    ret &= array_lmin(count, num) == 1;
    ret &= array_lmax(count, num) == 1;
    ret &= array_lsum(count, num) == num;

    arena_load(save);
    return num > 0 ? ret : true;
}

void mesh_reorder_nodes(MeshNodes *nodes, MeshCells *cells, const number *map)
{
    assert(is_valid(map, nodes->num));
    Arena save = arena_save();

    struct {
        number global;
        vector coord;
    } *node = arena_malloc(nodes->num, sizeof(*node));

    for (number i = 0; i < nodes->num; i++) {
        node[map[i]].global = nodes->global[i];
        node[map[i]].coord = nodes->coord[i];
    }

    for (number i = 0; i < nodes->num; i++) {
        nodes->global[i] = node[i].global;
        nodes->coord[i] = node[i].coord;
    }

    if (cells) {
        for (number i = 0; i < cells->num; i++) {
            for (number j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
                cells->node.idx[j] = map[cells->node.idx[j]];
            }
        }
    }

    arena_load(save);
}

void mesh_reorder_cells(MeshCells *cells, MeshFaces *faces, number beg, number end,
                        const number *map)
{
    assert(is_valid(map, end - beg));
    Arena save = arena_save();

    number tot = end - beg;
    struct {
        number num_nodes;
        number node[MAX_CELL_NODES];
        number num_cells;
        number cell[MAX_CELL_FACES];
    } *cell = arena_malloc(tot, sizeof(*cell));

    for (number num = 0, i = beg; i < end; i++, num++) {
        cell[map[num]].num_nodes = cells->node.off[i + 1] - cells->node.off[i];
        for (number k = 0, j = cells->node.off[i]; j < cells->node.off[i + 1]; j++, k++) {
            cell[map[num]].node[k] = cells->node.idx[j];
        }
        if (cells->cell.off && cells->cell.idx) {
            cell[map[num]].num_cells = cells->cell.off[i + 1] - cells->cell.off[i];
            for (number k = 0, j = cells->cell.off[i]; j < cells->cell.off[i + 1]; j++, k++) {
                cell[map[num]].cell[k] = cells->cell.idx[j];
            }
        }
    }

    for (number num = 0, i = beg; i < end; i++, num++) {
        cells->node.off[i + 1] = cells->node.off[i] + cell[num].num_nodes;
        for (number k = 0, j = cells->node.off[i]; j < cells->node.off[i + 1]; j++, k++) {
            cells->node.idx[j] = cell[num].node[k];
        }
        if (cells->cell.off && cells->cell.idx) {
            cells->cell.off[i + 1] = cells->cell.off[i] + cell[num].num_cells;
            for (number k = 0, j = cells->cell.off[i]; j < cells->cell.off[i + 1]; j++, k++) {
                cells->cell.idx[j] = cell[num].cell[k];
            }
        }
    }

    if (cells->cell.off && cells->cell.idx) {
        for (number i = 0; i < cells->num; i++) {
            for (number j = cells->cell.off[i]; j < cells->cell.off[i + 1]; j++) {
                number idx = cells->cell.idx[j];
                if (beg <= idx && idx < end) {
                    cells->cell.idx[j] = beg + map[idx - beg];
                }
            }
        }
    }

    if (faces) {
        for (number i = 0; i < faces->num; i++) {
            number left = faces->cell[i].left;
            if (beg <= left && left < end) {
                faces->cell[i].left = beg + map[left - beg];
            }
            number right = faces->cell[i].right;
            if (beg <= right && right < end) {
                faces->cell[i].right = beg + map[right - beg];
            }
        }
    }

    arena_load(save);
}

typedef struct {
    number key;
    number num;
    number node[MAX_FACE_NODES];
    Adjacent cell;
} Face;

static int cmp_face(const void *lhs_, const void *rhs_)
{
    const Face *lhs = lhs_;
    const Face *rhs = rhs_;
    return cmp_asc(lhs->key, rhs->key);
}

void mesh_reorder_faces(MeshFaces *faces, const number *key)
{
    Arena save = arena_save();

    Face *face = arena_malloc(faces->num, sizeof(*face));
    for (number i = 0; i < faces->num; i++) {
        face[i].key = key[i];
        face[i].num = faces->node.off[i + 1] - faces->node.off[i];
        for (number k = 0, j = faces->node.off[i]; j < faces->node.off[i + 1]; j++, k++) {
            face[i].node[k] = faces->node.idx[j];
        }
        face[i].cell = faces->cell[i];
    }
    qsort(face, faces->num, sizeof(*face), cmp_face);

    for (number i = 0; i < faces->num; i++) {
        faces->node.off[i + 1] = faces->node.off[i] + face[i].num;
        for (number k = 0, j = faces->node.off[i]; j < faces->node.off[i + 1]; j++, k++) {
            faces->node.idx[j] = face[i].node[k];
        }
        faces->cell[i] = face[i].cell;
    }

    arena_load(save);
}
