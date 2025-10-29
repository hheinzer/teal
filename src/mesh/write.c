#include "mesh.h"
#include "teal/arena.h"
#include "teal/assert.h"
#include "teal/h5io.h"
#include "teal/sync.h"

static void write_nodes(const MeshNodes *nodes, hid_t loc)
{
    hid_t group = h5io_group_create("nodes", loc);

    long num_inner = nodes->num_inner;
    long tot_inner = sync_lsum(num_inner);
    h5io_attribute_write(0, "num", &tot_inner, 1, H5IO_LONG, group);

    h5io_dataset_write("coord", nodes->coord, (hsize_t[]){num_inner, 3}, 2, H5IO_SCALAR, group);

    h5io_group_close(group);
}

static void write_node_graph(Graph node, const long *global, long num_cells, hid_t loc)
{
    Arena save = arena_save();

    hid_t group = h5io_group_create("node", loc);

    long num_off = num_cells + (sync.rank == 0);
    long num_idx = node.off[num_cells];

    long *off = arena_malloc(num_off, sizeof(*off));
    long offset = sync_lexsum(num_idx);
    for (long i = 0; i < num_off; i++) {
        off[i] = offset + node.off[i + (sync.rank != 0)];  // globalize offsets
    }
    h5io_dataset_write("off", off, (hsize_t[]){num_off}, 1, H5IO_LONG, group);

    long *idx = arena_malloc(num_idx, sizeof(*idx));
    for (long i = 0; i < num_idx; i++) {
        idx[i] = global[node.idx[i]];  // remap indices
    }
    h5io_dataset_write("idx", idx, (hsize_t[]){num_idx}, 1, H5IO_LONG, group);

    h5io_group_close(group);

    arena_load(save);
}

static void write_cells(const MeshNodes *nodes, const MeshCells *cells,
                        const MeshEntities *entities, hid_t loc)
{
    Arena save = arena_save();

    hid_t group = h5io_group_create("cells", loc);

    long num_inner = cells->num_inner;
    long num_outer = cells->num_ghost + cells->num_periodic;
    long num_cells = num_inner + num_outer;
    long tot_cells = sync_lsum(num_cells);
    h5io_attribute_write(0, "num", &tot_cells, 1, H5IO_LONG, group);

    write_node_graph(cells->node, nodes->global, num_cells, group);

    unsigned char *type = arena_malloc(num_cells, sizeof(*type));
    long *entity = arena_malloc(num_cells, sizeof(*entity));

    long *index = arena_malloc(num_cells, sizeof(*index));
    int *rank = arena_malloc(num_cells, sizeof(*rank));

    long num = 0;
    for (long i = 0; i < entities->num; i++) {
        for (long j = entities->cell_off[i]; j < entities->cell_off[i + 1]; j++) {
            long num_nodes = cells->node.off[j + 1] - cells->node.off[j];
            if (i < entities->num_inner) {
                enum { VTK_TETRA = 10, VTK_PYRAMID = 14, VTK_WEDGE = 13, VTK_HEXAHEDRON = 12 };
                switch (num_nodes) {
                    case 4: type[num] = VTK_TETRA; break;
                    case 5: type[num] = VTK_PYRAMID; break;
                    case 6: type[num] = VTK_WEDGE; break;
                    case 8: type[num] = VTK_HEXAHEDRON; break;
                    default: assert(false);
                }
            }
            else {
                enum { VTK_TRIANGLE = 5, VTK_QUAD = 9 };
                switch (num_nodes) {
                    case 3: type[num] = VTK_TRIANGLE; break;
                    case 4: type[num] = VTK_QUAD; break;
                    default: assert(false);
                }
            }
            entity[num] = i;
            index[num] = j;
            rank[num] = sync.rank;
            num += 1;
        }
    }
    assert(num == num_cells);

    h5io_dataset_write("type", type, (hsize_t[]){num_cells}, 1, H5IO_UCHAR, group);
    h5io_dataset_write("entity", entity, (hsize_t[]){num_cells}, 1, H5IO_LONG, group);

    h5io_dataset_write("index", index, (hsize_t[]){num_cells}, 1, H5IO_LONG, group);
    h5io_dataset_write("rank", rank, (hsize_t[]){num_cells}, 1, H5IO_INT, group);

    h5io_dataset_write("volume", cells->volume, (hsize_t[]){num_cells}, 1, H5IO_SCALAR, group);
    h5io_dataset_write("center", cells->center, (hsize_t[]){num_cells, 3}, 2, H5IO_SCALAR, group);
    h5io_dataset_write("projection", cells->projection, (hsize_t[]){num_cells, 3}, 2, H5IO_SCALAR,
                       group);

    h5io_group_close(group);

    arena_load(save);
}

static void write_entities(const MeshEntities *entities, hid_t loc)
{
    hid_t group = h5io_group_create("entities", loc);

    h5io_attribute_write(0, "num", &entities->num, 1, H5IO_LONG, group);
    h5io_attribute_write(0, "num_inner", &entities->num_inner, 1, H5IO_LONG, group);
    h5io_attribute_write(0, "num_ghost", &entities->num_ghost, 1, H5IO_LONG, group);

    long num_entities = (sync.rank == 0) ? entities->num : 0;
    h5io_dataset_write("name", entities->name, (hsize_t[]){num_entities, sizeof(*entities->name)},
                       2, H5IO_STRING, group);

    h5io_group_close(group);
}

void mesh_write(const Mesh *mesh, const char *prefix)
{
    assert(mesh && prefix);

    char fname[128];
    sprintf(fname, "%s_mesh.h5", prefix);

    hid_t file = h5io_file_create(fname);

    write_nodes(&mesh->nodes, file);
    write_cells(&mesh->nodes, &mesh->cells, &mesh->entities, file);
    write_entities(&mesh->entities, file);

    h5io_file_close(file);
}
