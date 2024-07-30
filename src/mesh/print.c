#include <stdio.h>

#include "core/array.h"
#include "core/sync.h"
#include "core/utils.h"
#include "mesh.h"
#include "teal.h"

static double compute_size(const Mesh *mesh);

void mesh_print(const Mesh *mesh)
{
    if (teal.quiet) return;

    const long n_inner_nodes = sync_sum(mesh->n_inner_nodes);
    const long n_inner_cells = sync_sum(mesh->n_inner_cells);
    const long n_inner_faces = sync_sum(mesh->n_inner_faces);
    const long n_bound_faces = sync_sum(mesh->n_bound_faces);
    const long n_sync_faces =
        sync_sum(mesh->n_faces - mesh->n_inner_faces - mesh->n_bound_faces) / 2;

    const double min_vol = sync_min(array_min(mesh->cell.volume, mesh->n_inner_cells));
    const double max_vol = sync_max(array_max(mesh->cell.volume, mesh->n_inner_cells));
    const double tot_vol = mesh->volume;

    double size = sync_sum(compute_size(mesh));
    const char mod = sizefmt(&size);

    if (teal.rank == 0) {
        printf("Mesh summary:\n");
        printf(" | " KEYFMT ": %ld\n", "number of inner nodes", n_inner_nodes);
        printf(" | " KEYFMT ": %ld\n", "number of inner cells", n_inner_cells);
        printf(" | " KEYFMT ": %ld\n", "number of inner faces", n_inner_faces);
        printf(" | " KEYFMT ": %ld\n", "number of boundary faces", n_bound_faces);
        printf(" | " KEYFMT ": %ld\n", "number of sync faces", n_sync_faces);
        printf(" | " KEYFMT ": %g / %g / %g\n", "min/max/tot volume", min_vol, max_vol, tot_vol);
        printf(" | " KEYFMT ": %g %cB\n", "memory size", size, mod);
    }
}

static double compute_size(const Mesh *mesh)
{
    double size = sizeof(*mesh);

    size += mesh->n_nodes * sizeof(*mesh->node.idx);
    size += mesh->n_nodes * sizeof(*mesh->node.coord);

    size += (mesh->n_cells + 1) * sizeof(*mesh->cell.i_node);
    size += mesh->cell.i_node[mesh->n_cells] * sizeof(*mesh->cell.node);
    size += (mesh->n_cells + 1) * sizeof(*mesh->cell.i_cell);
    size += mesh->cell.i_cell[mesh->n_cells] * sizeof(*mesh->cell.cell);
    size += mesh->n_inner_cells * sizeof(*mesh->cell.volume);
    size += mesh->n_cells * sizeof(*mesh->cell.center);
    size += mesh->n_inner_cells * sizeof(*mesh->cell.projection);
    size += mesh->cell.i_cell[mesh->n_inner_cells] * sizeof(*mesh->cell.reconstruction);

    size += (mesh->n_faces + 1) * sizeof(*mesh->face.i_node);
    size += mesh->face.i_node[mesh->n_faces] * sizeof(*mesh->face.node);
    size += mesh->n_faces * sizeof(*mesh->face.cell);
    size += mesh->n_faces * sizeof(*mesh->face.area);
    size += mesh->n_faces * sizeof(*mesh->face.center);
    size += mesh->n_faces * sizeof(*mesh->face.normal);
    size += mesh->n_faces * sizeof(*mesh->face.gradient_weight);
    size += mesh->n_faces * sizeof(*mesh->face.reconstruction);
    size += mesh->n_faces * sizeof(*mesh->face.gradient_correction);

    size += mesh->n_entities * sizeof(*mesh->entity.name);
    size += (mesh->n_entities + 1) * sizeof(*mesh->entity.j_cell);
    size += (mesh->n_entities + 1) * sizeof(*mesh->entity.j_face);
    size += mesh->n_entities * sizeof(*mesh->entity.offset);

    size += (teal.size + 1) * sizeof(*mesh->sync.j_recv);
    size += (teal.size + 1) * sizeof(*mesh->sync.i_send);
    size += mesh->sync.i_send[teal.size] * sizeof(*mesh->sync.send);

    return size;
}
