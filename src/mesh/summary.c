#include <assert.h>

#include "mesh.h"
#include "teal/sync.h"
#include "teal/utils.h"

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

    println("Mesh summary");
    println("\t number of inner nodes    : %ld", inner_nodes);
    println("\t number of inner cells    : %ld", inner_cells);
    println("\t number of inner faces    : %ld", inner_faces);
    println("\t number of ghost faces    : %ld", ghost_faces);
    println("\t number of outer faces    : %ld", outer_faces);
    println("\t number of periodic cells : %ld", periodic_cells);
    println("\t number of neighbor cells : %ld", neighbor_cells);
    println("\t number of exchange cells : %ld", exchange_cells);
}
