#include "reorder.h"

#include <stdlib.h>

#include "teal/arena.h"
#include "teal/dict.h"
#include "teal/utils.h"

void mesh_reorder_nodes(MeshNodes *nodes, MeshCells *cells, const long *key)
{
    Arena save = arena_save();

    struct {
        long key;
        long idx;
        long global;
        vector coord;
    } *buf = arena_malloc(nodes->num, sizeof(*buf));

    for (long i = 0; i < nodes->num; i++) {
        buf[i].key = key[i];
        buf[i].idx = i;
        buf[i].global = nodes->global[i];
        buf[i].coord = nodes->coord[i];
    }
    qsort(buf, nodes->num, sizeof(*buf), lcmp);

    for (long i = 0; i < nodes->num; i++) {
        nodes->global[i] = buf[i].global;
        nodes->coord[i] = buf[i].coord;
    }

    if (cells) {
        Dict *old2new = dict_create(sizeof(long), sizeof(long));
        for (long i = 0; i < nodes->num; i++) {
            dict_insert(old2new, &buf[i].idx, &i);
        }
        for (long i = 0; i < cells->num; i++) {
            for (long j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
                long *new = dict_lookup(old2new, &cells->node.idx[j]);
                if (new) {
                    cells->node.idx[j] = *new;
                }
            }
        }
    }

    arena_load(save);
}

void mesh_reorder_cells(MeshCells *cells, MeshFaces *faces, long beg, long end, const long *key)
{
    Arena save = arena_save();

    long tot = end - beg;
    struct {
        long key;
        long idx;
        long num_nodes;
        long node[MAX_CELL_NODES];
        long num_cells;
        long cell[MAX_CELL_FACES];
    } *buf = arena_malloc(tot, sizeof(*buf));

    for (long num = 0, i = beg; i < end; i++, num++) {
        buf[num].key = key[num];
        buf[num].idx = i;
        buf[num].num_nodes = cells->node.off[i + 1] - cells->node.off[i];
        for (long k = 0, j = cells->node.off[i]; j < cells->node.off[i + 1]; j++, k++) {
            buf[num].node[k] = cells->node.idx[j];
        }
        if (cells->cell.off && cells->cell.idx) {
            buf[num].num_cells = cells->cell.off[i + 1] - cells->cell.off[i];
            for (long k = 0, j = cells->cell.off[i]; j < cells->cell.off[i + 1]; j++, k++) {
                buf[num].cell[k] = cells->cell.idx[j];
            }
        }
    }
    qsort(buf, tot, sizeof(*buf), lcmp);

    Dict *old2new;
    if ((cells->cell.off && cells->cell.idx) || faces) {
        old2new = dict_create(sizeof(long), sizeof(long));
        for (long num = 0, i = beg; i < end; i++, num++) {
            dict_insert(old2new, &buf[num].idx, &i);
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
                long *new = dict_lookup(old2new, &buf[num].cell[k]);
                cells->cell.idx[j] = new ? *new : buf[num].cell[k];
            }
        }
    }

    if (cells->cell.off && cells->cell.idx) {
        for (long i = 0; i < cells->num; i++) {
            if (i < beg || end <= i) {
                for (long j = cells->cell.off[i]; j < cells->cell.off[i + 1]; j++) {
                    long *new = dict_lookup(old2new, &cells->cell.idx[j]);
                    if (new) {
                        cells->cell.idx[j] = *new;
                    }
                }
            }
        }
    }

    if (faces) {
        for (long i = 0; i < faces->num; i++) {
            long *new = dict_lookup(old2new, &faces->cell[i].left);
            if (new) {
                faces->cell[i].left = *new;
            }
            if ((new = dict_lookup(old2new, &faces->cell[i].right))) {
                faces->cell[i].right = *new;
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
