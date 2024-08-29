#include "teal/print.h"

#include <stdio.h>

#include "mesh.h"
#include "teal/array.h"
#include "teal/memory.h"
#include "teal/option.h"
#include "teal/sync.h"

void mesh_print(const Mesh *mesh)
{
    if (option.quiet) return;

    const long n_inner_nodes = sync_sum(mesh->n_inner_nodes);
    const long n_inner_cells = sync_sum(mesh->n_inner_cells);
    const long n_inner_faces = sync_sum(mesh->n_inner_faces);
    const long n_ghost_faces = sync_sum(mesh->n_ghost_faces);
    const long n_sync_faces =
        sync_sum(mesh->n_faces - mesh->n_inner_faces - mesh->n_ghost_faces) / 2;

    const double min_vol = sync_min(array_min(mesh->cell.volume, mesh->n_inner_cells));
    const double max_vol = sync_max(array_max(mesh->cell.volume, mesh->n_inner_cells));
    const double sum_vol = mesh->volume;

    const double size = sync_sum(memory_sum_get());

    if (sync.rank == 0) {
        printf("Mesh summary:\n");

        print_key("number of inner nodes", "%ld", n_inner_nodes);
        print_key("number of inner cells", "%ld", n_inner_cells);
        print_key("number of inner faces", "%ld", n_inner_faces);
        print_key("number of ghost faces", "%ld", n_ghost_faces);
        print_key("number of sync faces", "%ld", n_sync_faces);

        print_key("min/max/sum cell volume", "%g / %g / %g", min_vol, max_vol, sum_vol);

        print_size(size);
    }
}
