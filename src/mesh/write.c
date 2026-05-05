#include <assert.h>
#include <stdio.h>

#include "h5io.h"
#include "mesh.h"
#include "sync.h"
#include "teal.h"

enum VTKCellType {
    VTK_TRIANGLE = 5,
    VTK_QUAD = 9,
    VTK_TETRA = 10,
    VTK_PYRAMID = 14,
    VTK_WEDGE = 13,
    VTK_HEXAHEDRON = 12
};

enum CellGhostTypes {
    DUPLICATECELL = 1,  // present on multiple partitions/ranks
    EXTERIORCELL = 16,  // on the exterior of the data set
    HIDDENCELL = 32,    // only present for connectivity; ignore data values
};

static unsigned char *compute_types(const Mesh *mesh, int num_cells)
{
    unsigned char *type = teal_calloc(num_cells, sizeof(*type));
    for (int i = 0; i < num_cells; i++) {
        int num_nodes = mesh->cells.node.off[i + 1] - mesh->cells.node.off[i];
        if (i < mesh->cells.num_inner || mesh->cells.off_periodic <= i) {
            switch (num_nodes) {
                case 4: type[i] = VTK_TETRA; break;
                case 5: type[i] = VTK_PYRAMID; break;
                case 6: type[i] = VTK_WEDGE; break;
                case 8: type[i] = VTK_HEXAHEDRON; break;
                default: teal_error("invalid number of nodes (%d)", num_nodes);
            }
        }
        else {
            switch (num_nodes) {
                case 3: type[i] = VTK_TRIANGLE; break;
                case 4: type[i] = VTK_QUAD; break;
                default: teal_error("invalid number of nodes (%d)", num_nodes);
            }
        }
    }
    return type;
}

static void write_mesh_data_partitioned(const Mesh *mesh, int num_nodes, int num_cells, hid_t loc)
{
    int num_indices = mesh->cells.node.off[num_cells];

    h5io_dataset_write("NumberOfPoints", &(long){num_nodes}, 1, 1, H5T_NATIVE_LONG, loc);
    h5io_dataset_write("NumberOfCells", &(long){num_cells}, 1, 1, H5T_NATIVE_LONG, loc);
    h5io_dataset_write("NumberOfConnectivityIds", &(long){num_indices}, 1, 1, H5T_NATIVE_LONG, loc);

    unsigned char *type = compute_types(mesh, num_cells);

    h5io_dataset_write("Points", mesh->nodes.coord, num_nodes, 3, H5T_NATIVE_DOUBLE, loc);
    h5io_dataset_write("Offsets", mesh->cells.node.off, num_cells + 1, 1, H5T_NATIVE_INT, loc);
    h5io_dataset_write("Connectivity", mesh->cells.node.idx, num_indices, 1, H5T_NATIVE_INT, loc);
    h5io_dataset_write("Types", type, num_cells, 1, H5T_NATIVE_UCHAR, loc);

    teal_free(type);
}

static int *compute_offsets(const Mesh *mesh, int num_cells, int num_indices)
{
    int prefix = num_indices;
    sync_prefix(&prefix, 1, MPI_INT);

    int *offset = teal_calloc(num_cells + 1, sizeof(*offset));
    for (int i = 0; i < num_cells + 1; i++) {
        offset[i] = prefix + mesh->cells.node.off[i];
    }
    return offset;
}

static int *compute_indices(const Mesh *mesh, int num_indices)
{
    int *index = teal_calloc(num_indices, sizeof(*index));
    for (int i = 0; i < num_indices; i++) {
        assert(mesh->nodes.global[mesh->cells.node.idx[i]] <= INT_MAX);
        index[i] = (int)mesh->nodes.global[mesh->cells.node.idx[i]];
    }
    return index;
}

static void write_mesh_data(const Mesh *mesh, int num_nodes, int num_cells, hid_t loc)
{
    int root = (sync.rank == 0);

    long tot_nodes = num_nodes;
    sync_sum(&tot_nodes, 1, MPI_LONG);

    long tot_cells = num_cells;
    sync_sum(&tot_cells, 1, MPI_LONG);

    int num_indices = mesh->cells.node.off[num_cells];
    long tot_indices = num_indices;
    sync_sum(&tot_indices, 1, MPI_LONG);

    h5io_dataset_write("NumberOfPoints", &tot_nodes, root, 1, H5T_NATIVE_LONG, loc);
    h5io_dataset_write("NumberOfCells", &tot_cells, root, 1, H5T_NATIVE_LONG, loc);
    h5io_dataset_write("NumberOfConnectivityIds", &tot_indices, root, 1, H5T_NATIVE_LONG, loc);

    int *offset = compute_offsets(mesh, num_cells, num_indices);
    int *index = compute_indices(mesh, num_indices);
    unsigned char *type = compute_types(mesh, num_cells);

    h5io_dataset_write("Points", mesh->nodes.coord, num_nodes, 3, H5T_NATIVE_DOUBLE, loc);
    h5io_dataset_write("Offsets", &offset[!root], num_cells + root, 1, H5T_NATIVE_INT, loc);
    h5io_dataset_write("Connectivity", index, num_indices, 1, H5T_NATIVE_INT, loc);
    h5io_dataset_write("Types", type, num_cells, 1, H5T_NATIVE_UCHAR, loc);

    teal_free(offset);
    teal_free(index);
    teal_free(type);
}

static unsigned char *compute_ghosts(const Mesh *mesh, int num_cells)
{
    unsigned char *ghost = teal_calloc(num_cells, sizeof(*ghost));
    for (int i = mesh->cells.num_inner; i < num_cells; i++) {
        if (i < mesh->cells.off_periodic) {
            ghost[i] = EXTERIORCELL | HIDDENCELL;
        }
        else {
            ghost[i] = DUPLICATECELL;
        }
    }
    return ghost;
}

static int *compute_entities(const Mesh *mesh, int num_cells)
{
    int *entity = teal_calloc(num_cells, sizeof(*entity));
    for (int i = 0; i < mesh->entities.num; i++) {
        int beg = mesh->entities.cell_off[i];
        int end = mesh->entities.cell_off[i + 1];
        if (beg >= num_cells) {
            return entity;
        }
        if (end > num_cells) {
            end = num_cells;
        }
        for (int j = beg; j < end; j++) {
            entity[j] = i;
        }
    }
    for (int i = mesh->entities.cell_off[mesh->entities.num]; i < num_cells; i++) {
        entity[i] = mesh->entities.num;
    }
    return entity;
}

static int *compute_ranks(int num_cells)
{
    int *rank = teal_calloc(num_cells, sizeof(*rank));
    for (int i = 0; i < num_cells; i++) {
        rank[i] = sync.rank;
    }
    return rank;
}

static void write_cell_data(const Mesh *mesh, int num_cells, hid_t loc)
{
    unsigned char *ghost = compute_ghosts(mesh, num_cells);
    int *entity = compute_entities(mesh, num_cells);
    int *rank = compute_ranks(num_cells);

    hid_t group = h5io_group_create("CellData", loc);
    h5io_dataset_write("vtkGhostType", ghost, num_cells, 1, H5T_NATIVE_UCHAR, group);
    h5io_dataset_write("entity", entity, num_cells, 1, H5T_NATIVE_INT, group);
    h5io_dataset_write("rank", rank, num_cells, 1, H5T_NATIVE_INT, group);
    h5io_dataset_write("volume", mesh->cells.volume, num_cells, 1, H5T_NATIVE_DOUBLE, group);
    h5io_dataset_write("center", mesh->cells.center, num_cells, 3, H5T_NATIVE_DOUBLE, group);
    h5io_group_close(group);

    teal_free(ghost);
    teal_free(entity);
    teal_free(rank);
}

void mesh_write(const Mesh *mesh, const char *name)
{
    assert(mesh && name);

    String fname;
    sprintf(fname, "%s_mesh.hdf", name);

    hid_t file = h5io_file_create(fname);
    hid_t vtkhdf = h5io_group_create("VTKHDF", file);

    h5io_attribute_write("Version", (int[]){1, 0}, 2, H5T_NATIVE_INT, vtkhdf);
    h5io_attribute_write("Type", "UnstructuredGrid", 16, H5T_C_S1, vtkhdf);

    int num_cells = mesh->cells.off_periodic;
    if (teal.partitioned) {
        int num_nodes = mesh->nodes.num;
        write_mesh_data_partitioned(mesh, num_nodes, num_cells, vtkhdf);
    }
    else {
        int num_nodes = mesh->nodes.num_inner;
        write_mesh_data(mesh, num_nodes, num_cells, vtkhdf);
    }
    write_cell_data(mesh, num_cells, vtkhdf);

    h5io_group_close(vtkhdf);
    h5io_file_close(file);
}
