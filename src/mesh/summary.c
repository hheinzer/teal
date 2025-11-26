#include <assert.h>

#include "mesh.h"
#include "teal/array.h"
#include "teal/sync.h"
#include "teal/utils.h"
#include "teal/vector.h"

/* Count cells received from ranks other than self. */
static long count_exchange_cells(const MeshNeighbors *neighbors)
{
    long count = 0;
    for (long i = 0; i < neighbors->num; i++) {
        if (neighbors->rank[i] != sync.rank) {
            count += neighbors->recv_off[i + 1] - neighbors->recv_off[i];
        }
    }
    return count;
}

void mesh_summary(const Mesh *mesh)
{
    assert(mesh);

    long inner_nodes = sync_lsum(mesh->nodes.num_inner);
    long inner_cells = sync_lsum(mesh->cells.num_inner);
    long inner_faces = sync_lsum(mesh->faces.num_inner);
    long ghost_faces = sync_lsum(mesh->faces.off_ghost - mesh->faces.num_inner);
    long outer_faces = sync_lsum(mesh->faces.num - mesh->faces.off_ghost);

    long periodic_cells = sync_lsum(mesh->cells.off_periodic - mesh->cells.off_ghost);
    long neighbor_cells = sync_lsum(mesh->cells.num - mesh->cells.off_periodic);
    long exchange_cells = sync_lsum(count_exchange_cells(&mesh->neighbors));

    vector min_coord = sync_vector_min(vector_min(mesh->nodes.coord, mesh->nodes.num_inner));
    vector max_coord = sync_vector_max(vector_max(mesh->nodes.coord, mesh->nodes.num_inner));

    scalar min_volume = sync_fmin(array_fmin(mesh->cells.volume, mesh->cells.num_inner));
    scalar max_volume = sync_fmax(array_fmax(mesh->cells.volume, mesh->cells.num_inner));
    scalar sum_volume = sync_fsum(array_fsum(mesh->cells.volume, mesh->cells.num_inner));

    println("Mesh summary");
    println("\t number of inner nodes    : %ld", inner_nodes);
    println("\t number of inner cells    : %ld", inner_cells);
    println("\t number of inner faces    : %ld", inner_faces);
    println("\t number of ghost faces    : %ld", ghost_faces);
    println("\t number of outer faces    : %ld", outer_faces);
    println("\t number of periodic cells : %ld", periodic_cells);
    println("\t number of neighbor cells : %ld", neighbor_cells);
    println("\t number of exchange cells : %ld", exchange_cells);
    println("\t min/max coord            : %g %g %g / %g %g %g", min_coord.x, min_coord.y,
            min_coord.z, max_coord.x, max_coord.y, max_coord.z);
    println("\t min/max volume           : %g / %g / %g", min_volume, max_volume, sum_volume);
}
