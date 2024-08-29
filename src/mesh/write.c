#include <stdio.h>
#include <stdlib.h>

#include "mesh.h"
#include "teal/h5io.h"
#include "teal/memory.h"
#include "teal/sync.h"

void mesh_write(const Mesh *mesh, const char *prefix)
{
    // https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html#vtkhdf-file-format
    String fname;
    sprintf(fname, "%s_mesh.hdf", prefix);

    const int root = (sync.rank == 0);
    hid_t file = h5_file_create(fname);
    hid_t vtkhdf = h5_group_create("VTKHDF", file);

    const Vector2l version = {2, 2};
    h5_attribute_write("Version", version, 2, vtkhdf);
    h5_attribute_write("Type", "UnstructuredGrid", 0, vtkhdf);

    const long n_inner_nodes = mesh->n_inner_nodes;
    const long n_inner_cells = mesh->n_inner_cells;
    const long n_global_points = sync_sum(n_inner_nodes);
    const long n_conns = mesh->cell.i_node[n_inner_cells];
    const long o_conns = sync_exsum(n_conns);
    const long n_global_conns = sync_sum(n_conns);
    const long n_offsets = n_inner_cells + root;
    const long n_global_cells = sync_sum(n_inner_cells);

    h5_dataset_write("NumberOfPoints", &n_global_points, h5_dims(root), vtkhdf);
    h5_dataset_write("NumberOfConnectivityIds", &n_global_conns, h5_dims(root), vtkhdf);
    h5_dataset_write("NumberOfCells", &n_global_cells, h5_dims(root), vtkhdf);

    h5_dataset_write("Points", *mesh->node.coord, h5_dims(n_inner_nodes, N_DIMS), vtkhdf);

    smart long *conn = memory_calloc(n_conns, sizeof(*conn));
    for (long i = 0; i < n_conns; ++i) conn[i] = mesh->node.global[mesh->cell.node[i]];
    h5_dataset_write("Connectivity", conn, h5_dims(n_conns), vtkhdf);

    smart long *offsets = memory_calloc(n_offsets, sizeof(*offsets));
    for (long i = 0; i < n_offsets; ++i) offsets[i] = mesh->cell.i_node[!root + i] + o_conns;
    h5_dataset_write("Offsets", offsets, h5_dims(n_offsets), vtkhdf);

    smart unsigned char *types = memory_calloc(n_inner_cells, sizeof(*types));
    for (long i = 0; i < n_inner_cells; ++i) {
        const long n_nodes = mesh->cell.i_node[i + 1] - mesh->cell.i_node[i];
        switch (n_nodes) {
            case 4: types[i] = 10; break;  // VTK_TETRA
            case 5: types[i] = 14; break;  // VTK_PYRAMID
            case 6: types[i] = 13; break;  // VTK_WEDGE
            case 8: types[i] = 12; break;  // VTK_HEXAHEDRON
            default: abort();
        }
    }
    h5_dataset_write("Types", types, h5_dims(n_inner_cells), vtkhdf);

    hid_t cell = h5_group_create("CellData", vtkhdf);
    smart long *buf = memory_calloc(n_inner_cells, sizeof(*buf));

    for (long i = 0; i < n_inner_cells; ++i) buf[i] = sync.rank;
    h5_dataset_write("rank", buf, h5_dims(n_inner_cells), cell);

    for (long i = 0; i < n_inner_cells; ++i) buf[i] = i;
    h5_dataset_write("index", buf, h5_dims(n_inner_cells), cell);

    h5_dataset_write("volume", mesh->cell.volume, h5_dims(n_inner_cells), cell);
    h5_dataset_write("center", *mesh->cell.center, h5_dims(n_inner_cells, N_DIMS), cell);
    h5_dataset_write("projection", *mesh->cell.projection, h5_dims(n_inner_cells, N_DIMS), cell);

    h5_group_close(cell);
    h5_group_close(vtkhdf);
    h5_file_close(file);
}
