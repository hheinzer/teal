#include <assert.h>
#include <math.h>
#include <string.h>

#include "mesh.h"
#include "teal/arena.h"
#include "teal/h5io.h"
#include "teal/sync.h"
#include "teal/utils.h"
#include "teal/vector.h"

enum { VTK_TETRA = 10, VTK_PYRAMID = 14, VTK_WEDGE = 13, VTK_HEXAHEDRON = 12 };
enum { DUPLICATECELL = 1, EXTERIORCELL = 16, HIDDENCELL = 32 };

static void write_point_data1(const MeshNodes *nodes, hid_t loc)
{
    bool root = (sync.rank == 0);

    long num_nodes = nodes->num_inner;
    long tot_nodes = sync_lsum(num_nodes);
    h5io_dataset_write("NumberOfPoints", &tot_nodes, root, 1, H5IO_LONG, loc);

    h5io_dataset_write("Points", nodes->coord, num_nodes, 3, H5IO_SCALAR, loc);
}

static void write_cell_data1(const MeshNodes *nodes, const MeshCells *cells, hid_t loc)
{
    Arena save = arena_save();

    bool root = (sync.rank == 0);

    long num_cells = cells->num_inner;
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
        switch (num_nodes) {
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

void mesh_write1(const Mesh *mesh, const char *prefix)
{
    assert(mesh && prefix);

    char fname[128];
    sprintf(fname, "%s_mesh.vtkhdf", prefix);

    hid_t file = h5io_file_create(fname);
    hid_t vtkhdf = h5io_group_create("VTKHDF", file);

    h5io_attribute_write("Version", (long[]){1, 0}, 2, H5IO_LONG, vtkhdf);
    h5io_attribute_write("Type", "UnstructuredGrid", 1, H5IO_STRING, vtkhdf);

    write_point_data1(&mesh->nodes, vtkhdf);
    write_cell_data1(&mesh->nodes, &mesh->cells, vtkhdf);

    h5io_group_close(vtkhdf);
    h5io_file_close(file);
}

static bool is_flat(const MeshNodes *nodes, const MeshCells *cells, long idx)
{
    long num_nodes = cells->node.off[idx + 1] - cells->node.off[idx];
    if (num_nodes == 3) {
        return true;
    }
    if (num_nodes == 4) {
        vector coord[4];
        for (long k = 0, j = cells->node.off[idx]; j < cells->node.off[idx + 1]; j++, k++) {
            coord[k] = nodes->coord[cells->node.idx[j]];
        }
        vector a2b = vector_sub(coord[1], coord[0]);
        vector a2c = vector_sub(coord[2], coord[0]);
        vector a2d = vector_sub(coord[3], coord[0]);
        scalar volume = fabs(vector_dot(a2b, vector_cross(a2c, a2d))) / 6;
        return is_close(volume, 0);
    }
    return false;
}

static void write_point_data2(const MeshNodes *nodes, const MeshCells *cells,
                              const MeshFaces *faces, hid_t loc)
{
    Arena save = arena_save();

    long num_points = nodes->num;
    for (long i = cells->num_inner; i < cells->num; i++) {
        if (is_flat(nodes, cells, i)) {
            long num_nodes = cells->node.off[i + 1] - cells->node.off[i];
            num_points += num_nodes;
        }
    }
    h5io_dataset_write("NumberOfPoints", &num_points, 1, 1, H5IO_LONG, loc);

    vector *point = arena_malloc(num_points, sizeof(*point));
    memcpy(point, nodes->coord, nodes->num * sizeof(*point));
    long num = nodes->num;
    for (long i = faces->num_inner; i < faces->num; i++) {
        long right = faces->cell[i].right;
        if (is_flat(nodes, cells, right)) {
            vector offset = vector_sub(cells->center[right], faces->center[i]);
            for (long j = cells->node.off[right]; j < cells->node.off[right + 1]; j++) {
                vector coord = nodes->coord[cells->node.idx[j]];
                point[num++] = vector_add(coord, vector_mul(2, offset));
            }
        }
    }
    assert(num == num_points);
    h5io_dataset_write("Points", point, num_points, 3, H5IO_SCALAR, loc);

    arena_load(save);
}

static void write_cell_data2(const MeshNodes *nodes, const MeshCells *cells, hid_t loc)
{
    Arena save = arena_save();

    long num_cells = cells->num;
    h5io_dataset_write("NumberOfCells", &num_cells, 1, 1, H5IO_LONG, loc);

    long num_idx = cells->node.off[num_cells];
    for (long i = cells->num_inner; i < cells->num; i++) {
        if (is_flat(nodes, cells, i)) {
            long num_nodes = cells->node.off[i + 1] - cells->node.off[i];
            num_idx += num_nodes;
        }
    }
    h5io_dataset_write("NumberOfConnectivityIds", &num_idx, 1, 1, H5IO_LONG, loc);

    long num_off = num_cells + 1;
    long *off = arena_malloc(num_off, sizeof(*off));
    long *idx = arena_malloc(num_idx, sizeof(*idx));
    off[0] = 0;
    long num = nodes->num;
    for (long i = 0; i < cells->num; i++) {
        off[i + 1] = off[i];
        for (long j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
            idx[off[i + 1]++] = cells->node.idx[j];
        }
        if (is_flat(nodes, cells, i)) {
            for (long j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
                idx[off[i + 1]++] = num++;
            }
        }
    }
    h5io_dataset_write("Offsets", off, num_off, 1, H5IO_LONG, loc);
    h5io_dataset_write("Connectivity", idx, num_idx, 1, H5IO_LONG, loc);

    unsigned char *type = arena_malloc(num_cells, sizeof(*type));
    unsigned char *ghost = arena_malloc(num_cells, sizeof(*ghost));
    for (long i = 0; i < num_cells; i++) {
        long num_nodes = off[i + 1] - off[i];
        switch (num_nodes) {
            case 4: type[i] = VTK_TETRA; break;
            case 5: type[i] = VTK_PYRAMID; break;
            case 6: type[i] = VTK_WEDGE; break;
            case 8: type[i] = VTK_HEXAHEDRON; break;
            default: error("invalid number of nodes (%ld)", num_nodes);
        }
        if (i < cells->num_inner) {
            ghost[i] = 0;
        }
        else if (i < cells->off_periodic) {
            ghost[i] = EXTERIORCELL | HIDDENCELL;
        }
        else {
            ghost[i] = DUPLICATECELL;
        }
    }
    h5io_dataset_write("Types", type, num_cells, 1, H5T_NATIVE_UCHAR, loc);

    hid_t cell_data = h5io_group_create("CellData", loc);
    h5io_dataset_write("vtkGhostType", ghost, num_cells, 1, H5T_NATIVE_UCHAR, cell_data);
    h5io_group_close(cell_data);

    arena_load(save);
}

void mesh_write2(const Mesh *mesh, const char *prefix)
{
    assert(mesh && prefix);

    char fname[128];
    sprintf(fname, "%s_mesh.vtkhdf", prefix);

    hid_t file = h5io_file_create(fname);
    hid_t vtkhdf = h5io_group_create("VTKHDF", file);

    h5io_attribute_write("Version", (long[]){1, 0}, 2, H5IO_LONG, vtkhdf);
    h5io_attribute_write("Type", "UnstructuredGrid", 1, H5IO_STRING, vtkhdf);

    write_point_data2(&mesh->nodes, &mesh->cells, &mesh->faces, vtkhdf);
    write_cell_data2(&mesh->nodes, &mesh->cells, vtkhdf);

    h5io_group_close(vtkhdf);
    h5io_file_close(file);
}
