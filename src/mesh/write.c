#include "mesh.h"
#include "teal/arena.h"
#include "teal/assert.h"
#include "teal/h5io.h"
#include "teal/sync.h"

static void write_nodes(const MeshNodes *nodes, hid_t loc)
{
    hid_t group = h5io_group_create("nodes", loc);

    bool root = (sync.rank == 0);

    number num_nodes = nodes->num_inner;
    h5io_dataset_write("num", &num_nodes, 1, 1, H5IO_NUMBER, group);

    number tot_nodes = sync_lsum(num_nodes);
    h5io_dataset_write("tot", &tot_nodes, root, 1, H5IO_NUMBER, group);

    h5io_dataset_write("coord", nodes->coord, num_nodes, 3, H5IO_SCALAR, group);

    h5io_group_close(group);
}

static void write_node_graph(Graph node, const number *global, number num_cells, hid_t loc)
{
    Arena save = arena_save();

    hid_t group = h5io_group_create("node", loc);

    number num_off = num_cells + (sync.rank == 0);
    number num_idx = node.off[num_cells];

    number *off = arena_malloc(num_off, sizeof(*off));
    number offset = sync_lexsum(num_idx);
    for (number i = 0; i < num_off; i++) {
        off[i] = offset + node.off[i + (sync.rank != 0)];  // globalize offsets
    }
    h5io_dataset_write("off", off, num_off, 1, H5IO_NUMBER, group);

    number *idx = arena_malloc(num_idx, sizeof(*idx));
    for (number i = 0; i < num_idx; i++) {
        idx[i] = global[node.idx[i]];  // remap indices
    }
    h5io_dataset_write("idx", idx, num_idx, 1, H5IO_NUMBER, group);

    h5io_group_close(group);

    arena_load(save);
}

static void write_cells(const MeshNodes *nodes, const MeshCells *cells,
                        const MeshEntities *entities, hid_t loc)
{
    Arena save = arena_save();

    hid_t group = h5io_group_create("cells", loc);

    bool root = (sync.rank == 0);

    number num_cells = cells->off_periodic;
    h5io_dataset_write("num", &num_cells, 1, 1, H5IO_NUMBER, group);

    number tot_cells = sync_lsum(num_cells);
    h5io_dataset_write("tot", &tot_cells, root, 1, H5IO_NUMBER, group);

    number num_idx = cells->node.off[num_cells];
    number tot_idx = sync_lsum(num_idx);
    h5io_dataset_write("tot_idx", &tot_idx, root, 1, H5IO_NUMBER, group);

    write_node_graph(cells->node, nodes->global, num_cells, group);

    unsigned char *type = arena_malloc(num_cells, sizeof(*type));
    number *entity = arena_malloc(num_cells, sizeof(*entity));
    number *index = arena_malloc(num_cells, sizeof(*index));
    number *rank = arena_malloc(num_cells, sizeof(*rank));

    number num = 0;
    for (number i = 0; i < entities->num; i++) {
        for (number j = entities->cell_off[i]; j < entities->cell_off[i + 1]; j++) {
            number num_nodes = cells->node.off[j + 1] - cells->node.off[j];
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

    h5io_dataset_write("type", type, num_cells, 1, H5T_NATIVE_UCHAR, group);
    h5io_dataset_write("entity", entity, num_cells, 1, H5IO_NUMBER, group);
    h5io_dataset_write("index", index, num_cells, 1, H5IO_NUMBER, group);
    h5io_dataset_write("rank", rank, num_cells, 1, H5IO_NUMBER, group);

    h5io_dataset_write("volume", cells->volume, num_cells, 1, H5IO_SCALAR, group);
    h5io_dataset_write("center", cells->center, num_cells, 3, H5IO_SCALAR, group);
    h5io_dataset_write("projection", cells->projection, num_cells, 3, H5IO_SCALAR, group);

    h5io_group_close(group);

    arena_load(save);
}

static void write_entities(const MeshEntities *entities, hid_t loc)
{
    hid_t group = h5io_group_create("entities", loc);

    bool root = (sync.rank == 0);

    h5io_dataset_write("num", &entities->num, root, 1, H5IO_NUMBER, group);
    h5io_dataset_write("num_inner", &entities->num_inner, root, 1, H5IO_NUMBER, group);
    h5io_dataset_write("off_ghost", &entities->off_ghost, root, 1, H5IO_NUMBER, group);

    number num = root ? entities->num : 0;
    h5io_dataset_write("name", entities->name, num, sizeof(*entities->name), H5IO_STRING, group);

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
