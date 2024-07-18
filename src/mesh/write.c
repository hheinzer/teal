#include <stdio.h>

#include "core/h5io.h"
#include "core/memory.h"
#include "core/sync.h"
#include "core/utils.h"
#include "mesh.h"
#include "teal.h"

void mesh_write(const Mesh *mesh, const char *prefix)
{
    // https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html#vtkhdf-file-format
    char fname[128];
    sprintf(fname, "%s_mesh.hdf", prefix);

    hid_t file = h5_file_create(fname);
    hid_t vtkhdf = h5_group_create("VTKHDF", file);

    const long version[] = {2, 2};
    h5_attribute_write("Version", version, 2, vtkhdf);
    h5_attribute_write("Type", "UnstructuredGrid", 0, vtkhdf);

    const int root = (teal.rank == 0);
    const long n_inner_nodes = mesh->n_inner_nodes;
    const long n_inner_cells = mesh->n_inner_cells;
    const long n_global_points = sync_sum(n_inner_nodes);
    const long n_conns = mesh->cell.i_node[n_inner_cells];
    const long o_conns = sync_exscan_sum(n_conns);
    const long n_global_conns = sync_sum(n_conns);
    const long n_offsets = n_inner_cells + root;
    const long n_global_cells = sync_sum(n_inner_cells);

    h5_dataset_write("NumberOfPoints", &n_global_points, H5DIMS(root), vtkhdf);
    h5_dataset_write("NumberOfConnectivityIds", &n_global_conns, H5DIMS(root), vtkhdf);
    h5_dataset_write("NumberOfCells", &n_global_cells, H5DIMS(root), vtkhdf);

    h5_dataset_write("Points", *mesh->node.coord, H5DIMS(n_inner_nodes, N_DIMS), vtkhdf);

    cleanup long *conn = memory_calloc(n_conns, sizeof(*conn));
    for (long i = 0; i < n_conns; ++i) conn[i] = mesh->node.idx[mesh->cell.node[i]];
    h5_dataset_write("Connectivity", conn, H5DIMS(n_conns), vtkhdf);

    cleanup long *offsets = memory_calloc(n_offsets, sizeof(*offsets));
    for (long i = 0; i < n_offsets; ++i) offsets[i] = mesh->cell.i_node[!root + i] + o_conns;
    h5_dataset_write("Offsets", offsets, H5DIMS(n_offsets), vtkhdf);

    cleanup unsigned char *types = memory_calloc(n_inner_cells, sizeof(*types));
    for (long i = 0; i < n_inner_cells; ++i) {
        const long n_nodes = mesh->cell.i_node[i + 1] - mesh->cell.i_node[i];
        switch (n_nodes) {
            case 4: types[i] = 10; break;  // VTK_TETRA
            case 5: types[i] = 14; break;  // VTK_PYRAMID
            case 6: types[i] = 13; break;  // VTK_WEDGE
            case 8: types[i] = 12; break;  // VTK_HEXAHEDRON
            default: error("unsupported number of cell nodes '%ld'", n_nodes);
        }
    }
    h5_dataset_write("Types", types, H5DIMS(n_inner_cells), vtkhdf);

    hid_t cell = h5_group_create("CellData", vtkhdf);
    cleanup long *buf = memory_calloc(n_inner_cells, sizeof(*buf));

    for (long i = 0; i < n_inner_cells; ++i) buf[i] = teal.rank;
    h5_dataset_write("rank", buf, H5DIMS(n_inner_cells), cell);

    for (long i = 0; i < n_inner_cells; ++i) buf[i] = i;
    h5_dataset_write("index", buf, H5DIMS(n_inner_cells), cell);

    h5_dataset_write("volume", mesh->cell.volume, H5DIMS(n_inner_cells), cell);
    h5_dataset_write("center", *mesh->cell.center, H5DIMS(n_inner_cells, N_DIMS), cell);
    h5_dataset_write("projection", *mesh->cell.projection, H5DIMS(n_inner_cells, N_DIMS), cell);

    h5_group_close(cell);
    h5_group_close(vtkhdf);
    h5_file_close(file);
}
