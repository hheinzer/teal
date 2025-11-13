#include "mesh.h"
#include "teal/assert.h"
#include "teal/sync.h"
#include "teal/utils.h"

static number count_exchange_cells(const MeshNeighbors *neighbors)
{
    number count = 0;
    for (number i = 0; i < neighbors->num; i++) {
        if (neighbors->rank[i] != sync.rank) {
            count += neighbors->recv_off[i + 1] - neighbors->recv_off[i];
        }
    }
    return count;
}

void mesh_summary(const Mesh *mesh)
{
    assert(mesh);

    number inner_nodes = sync_lsum(mesh->nodes.num_inner);
    number inner_cells = sync_lsum(mesh->cells.num_inner);
    number inner_faces = sync_lsum(mesh->faces.num_inner);
    number ghost_faces = sync_lsum(mesh->faces.off_ghost - mesh->faces.num_inner);
    number outer_faces = sync_lsum(mesh->faces.num - mesh->faces.off_ghost);

    number periodic_cells = sync_lsum(mesh->cells.off_periodic - mesh->cells.off_ghost);
    number neighbor_cells = sync_lsum(mesh->cells.num - mesh->cells.off_periodic);
    number exchange_cells = sync_lsum(count_exchange_cells(&mesh->neighbors));

    println("Mesh summary");
    println("\t number of inner nodes    : %td", inner_nodes);
    println("\t number of inner cells    : %td", inner_cells);
    println("\t number of inner faces    : %td", inner_faces);
    println("\t number of ghost faces    : %td", ghost_faces);
    println("\t number of outer faces    : %td", outer_faces);
    println("\t number of periodic cells : %td", periodic_cells);
    println("\t number of neighbor cells : %td", neighbor_cells);
    println("\t number of exchange cells : %td", exchange_cells);
}
