#include "mesh.h"
#include "teal/array.h"
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
    long inner_nodes = sync_lsum(mesh->nodes.num_inner);
    long inner_cells = sync_lsum(mesh->cells.num_inner);
    long inner_faces = sync_lsum(mesh->faces.num_inner);
    long ghost_faces = sync_lsum(mesh->faces.num_ghost);
    long outer_faces = sync_lsum(mesh->faces.num - mesh->faces.num_inner - mesh->faces.num_ghost);

    long periodic_cells = sync_lsum(mesh->cells.num_periodic);
    long neighbor_cells = sync_lsum(mesh->cells.num - mesh->cells.num_inner -
                                    mesh->cells.num_ghost - mesh->cells.num_periodic);
    long exchange_cells = count_exchange_cells(&mesh->neighbors);

    vector min_coord = sync_vmin(array_vmin(mesh->nodes.coord, mesh->nodes.num_inner));
    vector max_coord = sync_vmax(array_vmax(mesh->nodes.coord, mesh->nodes.num_inner));

    double min_volume = sync_fmin(array_fmin(mesh->cells.volume, mesh->cells.num_inner));
    double max_volume = sync_fmax(array_fmax(mesh->cells.volume, mesh->cells.num_inner));
    double sum_volume = sync_fsum(array_fsum(mesh->cells.volume, mesh->cells.num_inner));

    print("Mesh summary:\n");
    print("\t number of inner nodes:    %ld\n", inner_nodes);
    print("\t number of inner cells:    %ld\n", inner_cells);
    print("\t number of inner faces:    %ld\n", inner_faces);
    print("\t number of ghost faces:    %ld\n", ghost_faces);
    print("\t number of outer faces:    %ld\n", outer_faces);
    print("\t number of periodic cells: %ld\n", periodic_cells);
    print("\t number of neighbor cells: %ld\n", neighbor_cells);
    print("\t number of exchange cells: %ld\n", exchange_cells);
    print("\t min/max coord:            %g %g %g / %g %g %g\n", min_coord.x, min_coord.y,
          min_coord.z, max_coord.x, max_coord.y, max_coord.z);
    print("\t min/max/sum volume:       %g / %g / %g\n", min_volume, max_volume, sum_volume);
}
