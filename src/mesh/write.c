#include <assert.h>

#include "mesh.h"
#include "teal/arena.h"
#include "teal/h5io.h"
#include "teal/utils.h"

void write_cell_data(const MeshCells *cells, hid_t loc)
{
    Arena save = arena_save();

    long num_cells = cells->num;
    long num_off = num_cells + 1;
    long num_idx = cells->node.off[num_cells];
    h5io_dataset_write("NumberOfCells", &num_cells, 1, 1, H5IO_LONG, loc);
    h5io_dataset_write("NumberOfConnectivityIds", &num_idx, 1, 1, H5IO_LONG, loc);
    h5io_dataset_write("Offsets", cells->node.off, num_off, 1, H5IO_LONG, loc);
    h5io_dataset_write("Connectivity", cells->node.idx, num_idx, 1, H5IO_LONG, loc);

    unsigned char *type = arena_malloc(num_cells, sizeof(*type));
    for (long i = 0; i < cells->num; i++) {
        long num_nodes = cells->node.off[i + 1] - cells->node.off[i];
        if (i < cells->num_inner || cells->off_periodic <= i) {
            enum { VTK_TETRA = 10, VTK_PYRAMID = 14, VTK_WEDGE = 13, VTK_HEXAHEDRON = 12 };
            switch (num_nodes) {
                case 4: type[i] = VTK_TETRA; break;
                case 5: type[i] = VTK_PYRAMID; break;
                case 6: type[i] = VTK_WEDGE; break;
                case 8: type[i] = VTK_HEXAHEDRON; break;
                default: error("invalid number of nodes (%ld)", num_nodes);
            }
        }
        else {
            enum { VTK_TRIANGLE = 5, VTK_QUAD = 9 };
            switch (num_nodes) {
                case 3: type[i] = VTK_TRIANGLE; break;
                case 4: type[i] = VTK_QUAD; break;
                default: error("invalid number of nodes (%ld)", num_nodes);
            }
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

    h5io_dataset_write("NumberOfPoints", &mesh->nodes.num, 1, 1, H5IO_LONG, vtkhdf);
    h5io_dataset_write("Points", mesh->nodes.coord, mesh->nodes.num, 3, H5IO_SCALAR, vtkhdf);

    write_cell_data(&mesh->cells, vtkhdf);

    h5io_group_close(vtkhdf);
    h5io_file_close(file);
}
