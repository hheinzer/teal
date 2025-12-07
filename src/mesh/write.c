#include <assert.h>

#include "mesh.h"
#include "teal/arena.h"
#include "teal/h5io.h"
#include "teal/utils.h"

static void write_point_data(const MeshNodes *nodes, hid_t loc)
{
    Arena save = arena_save();

    long num_nodes = nodes->num;
    h5io_dataset_write("NumberOfPoints", &num_nodes, 1, 1, H5IO_LONG, loc);

    h5io_dataset_write("Points", nodes->coord, num_nodes, 3, H5IO_SCALAR, loc);

    unsigned char *ghost = arena_malloc(num_nodes, sizeof(*ghost));
    for (long i = 0; i < num_nodes; i++) {
        ghost[i] = (i < nodes->num_inner) ? 0 : 1;
    }

    hid_t point_data = h5io_group_create("PointData", loc);
    h5io_dataset_write("GlobalNodeId", nodes->global, num_nodes, 1, H5IO_LONG, point_data);
    h5io_dataset_write("vtkGhostType", ghost, num_nodes, 1, H5T_NATIVE_UCHAR, point_data);
    h5io_group_close(point_data);

    arena_load(save);
}

static void write_cell_data(const MeshCells *cells, hid_t loc)
{
    Arena save = arena_save();

    long num_cells = cells->num_inner;
    long num_off = num_cells + 1;
    long num_idx = cells->node.off[num_cells];
    h5io_dataset_write("NumberOfCells", &num_cells, 1, 1, H5IO_LONG, loc);
    h5io_dataset_write("NumberOfConnectivityIds", &num_idx, 1, 1, H5IO_LONG, loc);

    h5io_dataset_write("Offsets", cells->node.off, num_off, 1, H5IO_LONG, loc);
    h5io_dataset_write("Connectivity", cells->node.idx, num_idx, 1, H5IO_LONG, loc);

    unsigned char *type = arena_malloc(num_cells, sizeof(*type));
    for (long i = 0; i < num_cells; i++) {
        long num_nodes = cells->node.off[i + 1] - cells->node.off[i];
        switch (num_nodes) {
            enum { VTK_TETRA = 10, VTK_PYRAMID = 14, VTK_WEDGE = 13, VTK_HEXAHEDRON = 12 };
            case 4: type[i] = VTK_TETRA; break;
            case 5: type[i] = VTK_PYRAMID; break;
            case 6: type[i] = VTK_WEDGE; break;
            case 8: type[i] = VTK_HEXAHEDRON; break;
            default: error("invalid number of nodes (%ld)", num_nodes);
        }
    }
    h5io_dataset_write("Types", type, num_cells, 1, H5T_NATIVE_UCHAR, loc);

    arena_load(save);
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
    write_cell_data(&mesh->cells, vtkhdf);

    h5io_group_close(vtkhdf);
    h5io_file_close(file);
}
