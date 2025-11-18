#include <stdlib.h>

#include "mesh.h"
#include "reorder.h"
#include "teal/arena.h"
#include "teal/assert.h"
#include "teal/h5io.h"
#include "teal/sync.h"

static void read_nodes(MeshNodes *nodes, hid_t loc)
{
    hid_t group = h5io_group_open("nodes", loc);

    bool root = (sync.rank == 0);

    int tot_nodes;
    h5io_dataset_read("tot", &tot_nodes, root, 1, H5IO_INT, group);
    MPI_Bcast(&tot_nodes, 1, MPI_INT, 0, sync.comm);

    int num_nodes;
    if (h5io_dataset_num("num", group) == sync.size) {
        h5io_dataset_read("num", &num_nodes, 1, 1, H5IO_INT, group);
    }
    else {
        num_nodes = (tot_nodes / sync.size) + (sync.rank < tot_nodes % sync.size);
    }
    nodes->num = num_nodes;
    assert(sync_lsum(nodes->num) == tot_nodes);

    nodes->coord = malloc(nodes->num * sizeof(*nodes->coord));
    assert(nodes->coord);

    h5io_dataset_read("coord", nodes->coord, nodes->num, 3, H5IO_SCALAR, group);

    h5io_group_close(group);
}

static void read_node_graph(Graph *node, int num_cells, hid_t loc)
{
    hid_t group = h5io_group_open("node", loc);

    node->off = malloc((num_cells + 1) * sizeof(*node->off));
    assert(node->off);

    node->off[0] = 0;

    int num_off = num_cells + (sync.rank == 0);
    h5io_dataset_read("off", &node->off[sync.rank != 0], num_off, 1, H5IO_INT, group);

    int offset = 0;
    int dst = (sync.rank + 1 < sync.size) ? (sync.rank + 1) : MPI_PROC_NULL;
    int src = (sync.rank - 1 >= 0) ? (sync.rank - 1) : MPI_PROC_NULL;
    MPI_Sendrecv(&node->off[num_cells], 1, MPI_INT, dst, 0, &offset, 1, MPI_INT, src, 0, sync.comm,
                 MPI_STATUS_IGNORE);

    for (int i = 0; i < num_cells; i++) {
        node->off[i + 1] -= offset;  // localize offsets
    }

    node->idx = malloc(node->off[num_cells] * sizeof(*node->idx));
    assert(node->idx);

    int num_idx = node->off[num_cells];
    h5io_dataset_read("idx", node->idx, num_idx, 1, H5IO_INT, group);

    h5io_group_close(group);
}

static void read_cells(MeshCells *cells, hid_t loc)
{
    hid_t group = h5io_group_open("cells", loc);

    bool root = (sync.rank == 0);

    int tot_cells;
    h5io_dataset_read("tot", &tot_cells, root, 1, H5IO_INT, group);
    MPI_Bcast(&tot_cells, 1, MPI_INT, 0, sync.comm);

    int num_cells;
    if (h5io_dataset_num("num", group) == sync.size) {
        h5io_dataset_read("num", &num_cells, 1, 1, H5IO_INT, group);
    }
    else {
        num_cells = (tot_cells / sync.size) + (sync.rank < tot_cells % sync.size);
    }
    cells->num = num_cells;
    assert(sync_lsum(cells->num) == tot_cells);

    read_node_graph(&cells->node, cells->num, group);

    h5io_group_close(group);
}

static void read_entities(MeshEntities *entities, hid_t loc)
{
    hid_t group = h5io_group_open("entities", loc);

    bool root = (sync.rank == 0);

    h5io_dataset_read("num", &entities->num, root, 1, H5IO_INT, group);
    h5io_dataset_read("num_inner", &entities->num_inner, root, 1, H5IO_INT, group);
    h5io_dataset_read("off_ghost", &entities->off_ghost, root, 1, H5IO_INT, group);

    MPI_Bcast(&entities->num, 1, MPI_INT, 0, sync.comm);
    MPI_Bcast(&entities->num_inner, 1, MPI_INT, 0, sync.comm);
    MPI_Bcast(&entities->off_ghost, 1, MPI_INT, 0, sync.comm);

    entities->name = malloc(entities->num * sizeof(*entities->name));
    assert(entities->name);

    int num = root ? entities->num : 0;
    h5io_dataset_read("name", entities->name, num, sizeof(*entities->name), H5IO_STRING, group);

    MPI_Bcast(entities->name, entities->num * sizeof(*entities->name), MPI_CHAR, 0, sync.comm);

    h5io_group_close(group);
}

static void reorder(MeshCells *cells, MeshEntities *entities, hid_t loc)
{
    Arena save = arena_save();

    int *entity = arena_malloc(cells->num, sizeof(*entity));
    h5io_dataset_read("cells/entity", entity, cells->num, 1, H5IO_INT, loc);

    int *num_cells = arena_calloc(entities->num, sizeof(*num_cells));
    for (int i = 0; i < cells->num; i++) {
        num_cells[entity[i]] += 1;
    }

    int *cell_off = malloc((entities->num + 1) * sizeof(*cell_off));
    assert(cell_off);

    cell_off[0] = 0;
    for (int i = 0; i < entities->num; i++) {
        cell_off[i + 1] = cell_off[i] + num_cells[i];
    }

    int *map = arena_malloc(cells->num, sizeof(*map));
    for (int i = 0; i < entities->num; i++) {
        cell_off[i + 1] -= num_cells[i];
    }
    for (int i = 0; i < cells->num; i++) {
        map[i] = cell_off[entity[i] + 1]++;
    }
    mesh_reorder_cells(cells, 0, 0, cells->num, map);

    arena_load(save);

    entities->cell_off = cell_off;
}

void mesh_read_hdf5(Mesh *mesh, const char *fname)
{
    hid_t file = h5io_file_open(fname);

    read_nodes(&mesh->nodes, file);
    read_cells(&mesh->cells, file);
    read_entities(&mesh->entities, file);

    reorder(&mesh->cells, &mesh->entities, file);

    h5io_file_close(file);
}
