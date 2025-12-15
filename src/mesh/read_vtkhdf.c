#include <assert.h>
#include <stdlib.h>

#include "mesh.h"
#include "reorder.h"
#include "teal/arena.h"
#include "teal/h5io.h"
#include "teal/sync.h"

static void read_point_data(MeshNodes *nodes, hid_t loc)
{
    bool root = (sync.rank == 0);

    long tot_nodes;
    h5io_dataset_read("NumberOfPoints", &tot_nodes, root, 1, H5IO_LONG, loc);
    MPI_Bcast(&tot_nodes, 1, MPI_LONG, 0, sync.comm);

    nodes->num = (tot_nodes / sync.size) + (sync.rank < tot_nodes % sync.size);
    assert(sync_lsum(nodes->num) == tot_nodes);

    nodes->coord = malloc(nodes->num * sizeof(*nodes->coord));
    assert(nodes->coord);

    h5io_dataset_read("Points", nodes->coord, nodes->num, 3, H5IO_SCALAR, loc);
}

static void read_cell_data(MeshCells *cells, hid_t loc)
{
    bool root = (sync.rank == 0);

    long tot_cells;
    h5io_dataset_read("NumberOfCells", &tot_cells, root, 1, H5IO_LONG, loc);
    MPI_Bcast(&tot_cells, 1, MPI_LONG, 0, sync.comm);

    cells->num = (tot_cells / sync.size) + (sync.rank < tot_cells % sync.size);
    assert(sync_lsum(cells->num) == tot_cells);

    cells->node.off = malloc((cells->num + 1) * sizeof(*cells->node.off));
    assert(cells->node.off);

    cells->node.off[0] = 0;

    long num_off = cells->num + root;
    h5io_dataset_read("Offsets", &cells->node.off[!root], num_off, 1, H5IO_LONG, loc);

    long offset = 0;
    long dst = (sync.rank + 1 < sync.size) ? (sync.rank + 1) : MPI_PROC_NULL;
    long src = (sync.rank - 1 >= 0) ? (sync.rank - 1) : MPI_PROC_NULL;
    MPI_Sendrecv(&cells->node.off[cells->num], 1, MPI_LONG, dst, 0, &offset, 1, MPI_LONG, src, 0,
                 sync.comm, MPI_STATUS_IGNORE);

    for (long i = 0; i < cells->num; i++) {
        cells->node.off[i + 1] -= offset;
    }

    long num_idx = cells->node.off[cells->num];
    cells->node.idx = malloc(num_idx * sizeof(*cells->node.idx));
    assert(cells->node.idx);

    h5io_dataset_read("Connectivity", cells->node.idx, num_idx, 1, H5IO_LONG, loc);
}

static void read_entities(MeshEntities *entities, hid_t loc)
{
    hid_t group = h5io_group_open("entities", loc);

    bool root = (sync.rank == 0);

    h5io_dataset_read("num", &entities->num, root, 1, H5IO_LONG, group);
    h5io_dataset_read("num_inner", &entities->num_inner, root, 1, H5IO_LONG, group);
    h5io_dataset_read("off_ghost", &entities->off_ghost, root, 1, H5IO_LONG, group);

    MPI_Bcast(&entities->num, 1, MPI_LONG, 0, sync.comm);
    MPI_Bcast(&entities->num_inner, 1, MPI_LONG, 0, sync.comm);
    MPI_Bcast(&entities->off_ghost, 1, MPI_LONG, 0, sync.comm);

    entities->name = malloc(entities->num * sizeof(*entities->name));
    assert(entities->name);

    long num = root ? entities->num : 0;
    h5io_dataset_read("name", entities->name, num, sizeof(*entities->name), H5IO_STRING, group);

    MPI_Bcast(entities->name, entities->num * sizeof(*entities->name), MPI_CHAR, 0, sync.comm);

    h5io_group_close(group);
}

static void reorder(MeshCells *cells, MeshEntities *entities, hid_t loc)
{
    Arena save = arena_save();

    long *entity = arena_malloc(cells->num, sizeof(*entity));
    h5io_dataset_read("VTKHDF/CellData/entity", entity, cells->num, 1, H5IO_LONG, loc);

    long *num_cells = arena_calloc(entities->num, sizeof(*num_cells));
    for (long i = 0; i < cells->num; i++) {
        num_cells[entity[i]] += 1;
    }

    long *cell_off = malloc((entities->num + 1) * sizeof(*cell_off));
    assert(cell_off);

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

    entities->cell_off = cell_off;
}

void mesh_read_vtkhdf(Mesh *mesh, const char *fname)
{
    hid_t file = h5io_file_open(fname);

    hid_t vtkhdf = h5io_group_open("VTKHDF", file);

    read_point_data(&mesh->nodes, vtkhdf);
    read_cell_data(&mesh->cells, vtkhdf);

    h5io_group_close(vtkhdf);

    read_entities(&mesh->entities, file);
    reorder(&mesh->cells, &mesh->entities, file);

    h5io_file_close(file);
}
