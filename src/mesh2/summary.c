#include <assert.h>

#include "mesh2.h"
#include "sync2.h"
#include "teal2.h"

void mesh2_summary(const Mesh *mesh)
{
    assert(mesh);

    long inner_nodes = mesh->nodes.num_inner;
    sync2_sum(&inner_nodes, 1, MPI_LONG);

    long inner_cells = mesh->cells.num_inner;
    sync2_sum(&inner_cells, 1, MPI_LONG);

    long sync_cells = 0;
    for (int i = 0; i < mesh->neighbors.num; i++) {
        if (mesh->neighbors.rank[i] != sync2.rank) {
            sync_cells += mesh->neighbors.recv_off[i + 1] - mesh->neighbors.recv_off[i];
        }
    }
    sync2_sum(&sync_cells, 1, MPI_LONG);

    teal2_print("Mesh summary");
    teal2_print("\t number of inner nodes : %ld", inner_nodes);
    teal2_print("\t number of inner cells : %ld", inner_cells);
    teal2_print("\t number of sync cells  : %ld", sync_cells);
}
