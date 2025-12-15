#include <assert.h>

#include "mesh.h"
#include "teal/arena.h"
#include "teal/h5io.h"
#include "teal/sync.h"
#include "teal/utils.h"

enum {
    VTK_TRIANGLE = 5,
    VTK_QUAD = 9,
    VTK_TETRA = 10,
    VTK_PYRAMID = 14,
    VTK_WEDGE = 13,
    VTK_HEXAHEDRON = 12
};

static void write_point_data(const MeshNodes *nodes, hid_t loc)
{
    bool root = (sync.rank == 0);

    long num_nodes = nodes->num_inner;
    long tot_nodes = sync_lsum(num_nodes);
    h5io_dataset_write("NumberOfPoints", &tot_nodes, root, 1, H5IO_LONG, loc);

    h5io_dataset_write("Points", nodes->coord, num_nodes, 3, H5IO_SCALAR, loc);
}

static void write_cell_data(const MeshNodes *nodes, const MeshCells *cells,
                            const MeshEntities *entities, hid_t loc)
{
    Arena save = arena_save();

    bool root = (sync.rank == 0);

    long num_cells = cells->off_periodic;
    long tot_cells = sync_lsum(num_cells);
    h5io_dataset_write("NumberOfCells", &tot_cells, root, 1, H5IO_LONG, loc);

    long num_idx = cells->node.off[num_cells];
    long tot_idx = sync_lsum(num_idx);
    h5io_dataset_write("NumberOfConnectivityIds", &tot_idx, root, 1, H5IO_LONG, loc);

    long num_off = num_cells + root;
    long *off = arena_malloc(num_off, sizeof(*off));
    long offset = sync_exsum(num_idx);
    for (long i = 0; i < num_off; i++) {
        off[i] = offset + cells->node.off[i + !root];
    }
    h5io_dataset_write("Offsets", off, num_off, 1, H5IO_LONG, loc);

    long *idx = arena_malloc(num_idx, sizeof(*idx));
    for (long i = 0; i < num_idx; i++) {
        idx[i] = nodes->global[cells->node.idx[i]];
    }
    h5io_dataset_write("Connectivity", idx, num_idx, 1, H5IO_LONG, loc);

    unsigned char *type = arena_malloc(num_cells, sizeof(*type));
    for (long i = 0; i < num_cells; i++) {
        long num_nodes = cells->node.off[i + 1] - cells->node.off[i];
        if (i < cells->num_inner) {
            switch (num_nodes) {
                case 4: type[i] = VTK_TETRA; break;
                case 5: type[i] = VTK_PYRAMID; break;
                case 6: type[i] = VTK_WEDGE; break;
                case 8: type[i] = VTK_HEXAHEDRON; break;
                default: error("invalid number of nodes (%ld)", num_nodes);
            }
        }
        else {
            switch (num_nodes) {
                case 3: type[i] = VTK_TRIANGLE; break;
                case 4: type[i] = VTK_QUAD; break;
                default: error("invalid number of nodes (%ld)", num_nodes);
            }
        }
    }
    h5io_dataset_write("Types", type, num_cells, 1, H5T_NATIVE_UCHAR, loc);

    long *entity = arena_malloc(num_cells, sizeof(*entity));
    for (long i = 0; i < entities->num; i++) {
        for (long j = entities->cell_off[i]; j < entities->cell_off[i + 1]; j++) {
            entity[j] = i;
        }
    }

    hid_t cell_data = h5io_group_create("CellData", loc);
    h5io_dataset_write("entity", entity, num_cells, 1, H5IO_LONG, cell_data);
    h5io_group_close(cell_data);

    arena_load(save);
}

static void write_entities(const MeshEntities *entities, hid_t loc)
{
    hid_t group = h5io_group_create("entities", loc);

    bool root = (sync.rank == 0);

    h5io_dataset_write("num", &entities->num, root, 1, H5IO_LONG, group);
    h5io_dataset_write("num_inner", &entities->num_inner, root, 1, H5IO_LONG, group);
    h5io_dataset_write("off_ghost", &entities->off_ghost, root, 1, H5IO_LONG, group);

    long num = root ? entities->num : 0;
    h5io_dataset_write("name", entities->name, num, sizeof(*entities->name), H5IO_STRING, group);

    h5io_group_close(group);
}

void mesh_write(const Mesh *mesh, const char *prefix)
{
    assert(mesh && prefix);

    char fname[128];
    sprintf(fname, "%s_mesh.vtkhdf", prefix);

    hid_t file = h5io_file_create(fname);

    hid_t vtkhdf = h5io_group_create("VTKHDF", file);

    h5io_attribute_write("Version", (long[]){1, 0}, 2, H5IO_LONG, vtkhdf);
    h5io_attribute_write("Type", "UnstructuredGrid", 1, H5IO_STRING, vtkhdf);

    write_point_data(&mesh->nodes, vtkhdf);
    write_cell_data(&mesh->nodes, &mesh->cells, &mesh->entities, vtkhdf);

    h5io_group_close(vtkhdf);

    write_entities(&mesh->entities, file);

    h5io_file_close(file);
}
