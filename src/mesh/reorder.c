#include "reorder.h"

#include <stdlib.h>

#include "teal/arena.h"
#include "teal/array.h"
#include "teal/utils.h"

static bool is_valid(const long *map, long num)
{
    Arena save = arena_save();

    long *cnt = arena_calloc(num, sizeof(*cnt));
    for (long i = 0; i < num; i++) {
        if (0 <= map[i] && map[i] < num) {
            cnt[map[i]] += 1;
        }
    }

    bool ret = true;
    ret &= array_lmin(cnt, num) == 1;
    ret &= array_lmax(cnt, num) == 1;
    ret &= array_lsum(cnt, num) == num;

    arena_load(save);
    return ret;
}

void mesh_reorder_nodes(MeshNodes *nodes, MeshCells *cells, const long *map)
{
    assert(is_valid(map, nodes->num));

    Arena save = arena_save();

    struct {
        long global;
        vector coord;
    } *buf = arena_malloc(nodes->num, sizeof(*buf));

    for (long i = 0; i < nodes->num; i++) {
        buf[map[i]].global = nodes->global[i];
        buf[map[i]].coord = nodes->coord[i];
    }

    for (long i = 0; i < nodes->num; i++) {
        nodes->global[i] = buf[i].global;
        nodes->coord[i] = buf[i].coord;
    }

    if (cells) {
        for (long i = 0; i < cells->num; i++) {
            for (long j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
                cells->node.idx[j] = map[cells->node.idx[j]];
            }
        }
    }

    arena_load(save);
}

void mesh_reorder_cells(MeshCells *cells, MeshFaces *faces, long beg, long end, const long *map)
{
    assert(is_valid(map, end - beg));

    Arena save = arena_save();

    long tot = end - beg;
    struct {
        long num_nodes;
        long node[MAX_CELL_NODES];
        long num_cells;
        long cell[MAX_CELL_FACES];
    } *buf = arena_malloc(tot, sizeof(*buf));

    for (long num = 0, i = beg; i < end; i++, num++) {
        buf[map[num]].num_nodes = cells->node.off[i + 1] - cells->node.off[i];
        for (long k = 0, j = cells->node.off[i]; j < cells->node.off[i + 1]; j++, k++) {
            buf[map[num]].node[k] = cells->node.idx[j];
        }
        if (cells->cell.off && cells->cell.idx) {
            buf[map[num]].num_cells = cells->cell.off[i + 1] - cells->cell.off[i];
            for (long k = 0, j = cells->cell.off[i]; j < cells->cell.off[i + 1]; j++, k++) {
                buf[map[num]].cell[k] = cells->cell.idx[j];
            }
        }
    }

    for (long num = 0, i = beg; i < end; i++, num++) {
        cells->node.off[i + 1] = cells->node.off[i] + buf[num].num_nodes;
        for (long k = 0, j = cells->node.off[i]; j < cells->node.off[i + 1]; j++, k++) {
            cells->node.idx[j] = buf[num].node[k];
        }
        if (cells->cell.off && cells->cell.idx) {
            cells->cell.off[i + 1] = cells->cell.off[i] + buf[num].num_cells;
            for (long k = 0, j = cells->cell.off[i]; j < cells->cell.off[i + 1]; j++, k++) {
                cells->cell.idx[j] = buf[num].cell[k];
            }
        }
    }

    if (cells->cell.off && cells->cell.idx) {
        for (long i = 0; i < cells->num; i++) {
            for (long j = cells->cell.off[i]; j < cells->cell.off[i + 1]; j++) {
                long idx = cells->cell.idx[j];
                if (beg <= idx && idx < end) {
                    cells->cell.idx[j] = beg + map[idx - beg];
                }
            }
        }
    }

    if (faces) {
        for (long i = 0; i < faces->num; i++) {
            long left = faces->cell[i].left;
            if (beg <= left && left < end) {
                faces->cell[i].left = beg + map[left - beg];
            }
            long right = faces->cell[i].right;
            if (beg <= right && right < end) {
                faces->cell[i].right = beg + map[right - beg];
            }
        }
    }

    arena_load(save);
}

void mesh_reorder_faces(MeshFaces *faces, const long *key)
{
    Arena save = arena_save();

    struct {
        long key;
        long num;
        long node[MAX_FACE_NODES];
        FaceCells cell;
    } *buf = arena_malloc(faces->num, sizeof(*buf));

    for (long i = 0; i < faces->num; i++) {
        buf[i].key = key[i];
        buf[i].num = faces->node.off[i + 1] - faces->node.off[i];
        for (long k = 0, j = faces->node.off[i]; j < faces->node.off[i + 1]; j++, k++) {
            buf[i].node[k] = faces->node.idx[j];
        }
        buf[i].cell = faces->cell[i];
    }
    qsort(buf, faces->num, sizeof(*buf), lcmp);

    for (long i = 0; i < faces->num; i++) {
        faces->node.off[i + 1] = faces->node.off[i] + buf[i].num;
        for (long k = 0, j = faces->node.off[i]; j < faces->node.off[i + 1]; j++, k++) {
            faces->node.idx[j] = buf[i].node[k];
        }
        faces->cell[i] = buf[i].cell;
    }

    arena_load(save);
}
