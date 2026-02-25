#include "grid.h"

#include <assert.h>
#include <string.h>

#include "sync2.h"
#include "teal2.h"
#include "utils2.h"

Grid *grid_init(const char *fname)
{
    assert(fname);
    char *ext = strrchr(fname, '.');
    if (!ext) {
        teal2_error("missing file extension (%s)", fname);
    }
    if (!strcmp(ext, ".msh")) {
        return grid_gmsh(fname);
    }
    teal2_error("unsupported file extension (%s)", ext);
}

void grid_deinit(Grid *grid)
{
    assert(grid);

    teal2_free(grid->nodes.tag);
    teal2_free(grid->nodes.coord);

    teal2_free(grid->cells.node_off);
    teal2_free(grid->cells.node_idx);

    teal2_free(grid->entities.name);
    teal2_free(grid->entities.cell_off);
    teal2_free(grid->entities.periodic);
    teal2_free(grid->entities.rotation);
    teal2_free(grid->entities.translation);
    teal2_free(grid->entities.node_off);
    teal2_free(grid->entities.node_idx);

    teal2_free(grid);
}

typedef struct {
    long tag, idx;
} Map;

static MPI_Datatype datatype_map(void)
{
    MPI_Datatype datatype;
    int len[2] = {1, 1};
    MPI_Aint off[2] = {offsetof(Map, tag), offsetof(Map, idx)};
    MPI_Datatype type[2] = {MPI_LONG, MPI_LONG};
    MPI_Type_create_struct(2, len, off, type, &datatype);
    return sync2_resized(datatype, sizeof(Map));
}

static int compare_map(const void *lhs_, const void *rhs_)
{
    const Map *lhs = lhs_;
    const Map *rhs = rhs_;
    return (lhs->tag > rhs->tag) - (lhs->tag < rhs->tag);
}

long *grid_tag_to_idx(const Grid *grid, const long *node_tag, int num_tags)
{
    int cap = grid->nodes.num;
    sync2_max(&cap, 1, MPI_INT);

    long prefix = grid->nodes.num;
    sync2_prefix(&prefix, 1, MPI_LONG);

    Map *map = teal2_calloc(cap, sizeof(*map));
    for (int i = 0; i < grid->nodes.num; i++) {
        map[i].tag = grid->nodes.tag[i];
        map[i].idx = prefix + i;
    }
    sort(map, grid->nodes.num, sizeof(*map), compare_map);

    int *converted = teal2_calloc(num_tags, sizeof(*converted));
    long *node_idx = teal2_calloc(num_tags, sizeof(*node_idx));

    MPI_Datatype type = datatype_map();
    int num = grid->nodes.num;
    for (int rank = 0; rank < sync2.size; rank++) {
        for (int i = 0; i < num_tags; i++) {
            if (!converted[i]) {
                Map key = {.tag = node_tag[i]};
                Map *val = search(&key, map, num, sizeof(*map), compare_map);
                if (val) {
                    converted[i] = 1;
                    node_idx[i] = val->idx;
                }
            }
        }
        sync2_rotate(map, &num, cap, type, 1);
    }
    assert(num == grid->nodes.num);
    MPI_Type_free(&type);

    for (int i = 0; i < num_tags; i++) {
        if (!converted[i]) {
            teal2_error("node tag (%ld) was not converted to index", node_tag[i]);
        }
    }

    teal2_free(map);
    teal2_free(converted);
    return node_idx;
}
