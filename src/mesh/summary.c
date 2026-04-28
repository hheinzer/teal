#include <assert.h>

#include "mesh.h"
#include "sync.h"
#include "teal.h"

void mesh_summary(const Mesh *mesh)
{
    assert(mesh);

    long inner_nodes = mesh->nodes.num_inner;
    sync_sum(&inner_nodes, 1, MPI_LONG);

    long inner_cells = mesh->cells.num_inner;
    sync_sum(&inner_cells, 1, MPI_LONG);

    long sync_cells = 0;
    for (int i = 0; i < mesh->neighbors.num; i++) {
        if (mesh->neighbors.rank[i] != sync.rank) {
            sync_cells += mesh->neighbors.recv_off[i + 1] - mesh->neighbors.recv_off[i];
        }
    }
    sync_sum(&sync_cells, 1, MPI_LONG);

    teal_print("Mesh summary");
    teal_print("\t number of inner nodes : %ld", inner_nodes);
    teal_print("\t number of inner cells : %ld", inner_cells);
    teal_print("\t number of sync cells  : %ld", sync_cells);
}
