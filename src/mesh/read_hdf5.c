#include <stdlib.h>

#include "mesh.h"
#include "reorder.h"
#include "teal/arena.h"
#include "teal/h5io.h"
#include "teal/sync.h"
#include "teal/utils.h"

static void read_nodes(MeshNodes *nodes, hid_t loc)
{
    hid_t group = h5io_group_open("nodes", loc);

    long tot_nodes;
    h5io_attribute_read(0, "num", &tot_nodes, 1, H5IO_LONG, group);

    nodes->num = (tot_nodes / sync.size) + (sync.rank < tot_nodes % sync.size);
    assert(tot_nodes == sync_lsum(nodes->num));

    nodes->coord = malloc(nodes->num * sizeof(*nodes->coord));
    assert(nodes->coord);

    h5io_dataset_read("coord", nodes->coord, (hsize_t[]){nodes->num, 3}, 2, H5IO_SCALAR, group);

    h5io_group_close(group);
}

static void read_node_graph(MeshGraph *node, long num_cells, hid_t loc)
{
    hid_t group = h5io_group_open("node", loc);

    node->off = malloc((num_cells + 1) * sizeof(*node->off));
    assert(node->off);

    node->off[0] = 0;

    long num = num_cells + (sync.rank == 0);
    h5io_dataset_read("off", &node->off[sync.rank != 0], (hsize_t[]){num}, 1, H5IO_LONG, group);

    long offset = 0;
    int dst = (sync.rank + 1 < sync.size) ? sync.rank + 1 : MPI_PROC_NULL;
    int src = (sync.rank - 1 >= 0) ? sync.rank - 1 : MPI_PROC_NULL;
    MPI_Sendrecv(&node->off[num_cells], 1, MPI_LONG, dst, 0, &offset, 1, MPI_LONG, src, 0,
                 sync.comm, MPI_STATUS_IGNORE);

    for (long i = 0; i < num_cells; i++) {
        node->off[i + 1] -= offset;  // localize offsets
    }

    node->idx = malloc(node->off[num_cells] * sizeof(*node->idx));
    assert(node->idx);

    h5io_dataset_read("idx", node->idx, (hsize_t[]){node->off[num_cells]}, 1, H5IO_LONG, group);

    h5io_group_close(group);
}

static void read_cells(MeshCells *cells, hid_t loc)
{
    hid_t group = h5io_group_open("cells", loc);

    long tot_cells;
    h5io_attribute_read(0, "num", &tot_cells, 1, H5IO_LONG, group);

    cells->num = (tot_cells / sync.size) + (sync.rank < tot_cells % sync.size);
    assert(tot_cells == sync_lsum(cells->num));

    read_node_graph(&cells->node, cells->num, group);

    h5io_group_close(group);
}

static void read_entities(MeshEntities *entities, hid_t loc)
{
    hid_t group = h5io_group_open("entities", loc);

    h5io_attribute_read(0, "num", &entities->num, 1, H5IO_LONG, group);
    h5io_attribute_read(0, "num_inner", &entities->num_inner, 1, H5IO_LONG, group);
    h5io_attribute_read(0, "num_ghost", &entities->num_ghost, 1, H5IO_LONG, group);

    entities->name = arena_malloc(entities->num, sizeof(*entities->name));

    long num_entities = (sync.rank == 0) ? entities->num : 0;
    h5io_dataset_read("name", entities->name, (hsize_t[]){num_entities}, 1, H5IO_STRING, group);

    MPI_Bcast(entities->name, entities->num * sizeof(*entities->name), MPI_CHAR, 0, sync.comm);

    h5io_group_close(group);
}

static void reorder(MeshCells *cells, MeshEntities *entities, hid_t loc)
{
    Arena save = arena_save();

    long *entity = arena_malloc(cells->num, sizeof(*entity));
    h5io_dataset_read("cells/entity", entity, (hsize_t[]){cells->num}, 1, H5IO_LONG, loc);

    long *num_cells = arena_calloc(entities->num, sizeof(*num_cells));
    for (long i = 0; i < cells->num; i++) {
        num_cells[entity[i]] += 1;
    }

    long *cell_off = arena_malloc(entities->num + 1, sizeof(*cell_off));
    cell_off[0] = 0;
    for (long i = 0; i < entities->num; i++) {
        cell_off[i + 1] = cell_off[i] + num_cells[i];
    }

    long *map = arena_malloc(cells->num, sizeof(*map));
    for (long i = 0; i < entities->num; i++) {
        cell_off[i + 1] -= num_cells[i];
    }
    for (long i = 0; i < cells->num; i++) {
        map[i] = cell_off[entity[i] + 1]++;
    }
    mesh_reorder_cells(cells, 0, 0, cells->num, map);

    arena_load(save);

    entities->cell_off = arena_smuggle(cell_off, entities->num + 1, sizeof(*cell_off));
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
