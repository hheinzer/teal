#include <assert.h>

#include "mesh.h"
#include "teal/sync.h"
#include "teal/utils.h"

static int count_exchange_cells(const MeshNeighbors *neighbors)
{
    int count = 0;
    for (int i = 0; i < neighbors->num; i++) {
        if (neighbors->rank[i] != sync.rank) {
            count += neighbors->recv_off[i + 1] - neighbors->recv_off[i];
        }
    }
    return count;
}

void mesh_summary(const Mesh *mesh)
{
    assert(mesh);

    int inner_nodes = sync_lsum(mesh->nodes.num_inner);
    int inner_cells = sync_lsum(mesh->cells.num_inner);
    int inner_faces = sync_lsum(mesh->faces.num_inner);
    int ghost_faces = sync_lsum(mesh->faces.off_ghost - mesh->faces.num_inner);
    int outer_faces = sync_lsum(mesh->faces.num - mesh->faces.off_ghost);

    int periodic_cells = sync_lsum(mesh->cells.off_periodic - mesh->cells.off_ghost);
    int neighbor_cells = sync_lsum(mesh->cells.num - mesh->cells.off_periodic);
    int exchange_cells = sync_lsum(count_exchange_cells(&mesh->neighbors));

    println("Mesh summary");
    println("\t number of inner nodes    : %d", inner_nodes);
    println("\t number of inner cells    : %d", inner_cells);
    println("\t number of inner faces    : %d", inner_faces);
    println("\t number of ghost faces    : %d", ghost_faces);
    println("\t number of outer faces    : %d", outer_faces);
    println("\t number of periodic cells : %d", periodic_cells);
    println("\t number of neighbor cells : %d", neighbor_cells);
    println("\t number of exchange cells : %d", exchange_cells);
}
