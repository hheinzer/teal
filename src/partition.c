#include "partition.h"

#include <assert.h>
#include <metis.h>
#include <mpi.h>
#include <stdlib.h>

#include "array.h"
#include "dict.h"
#include "memory.h"
#include "mesh.h"

static void compute_cell_partitioning(const Mesh *mesh, long *part);
static void reorder_cells(Mesh *mesh, long *part);
static void compute_cell_map(const Mesh *mesh, const idx_t *part, long (*cell_map)[2],
                             long *n_inner_cells, long *n_ghost_cells, long *n_cells);
static void compute_node_partitioning(const Mesh *mesh, const long *cpart, long *npart);
static void compute_node_map(Mesh *mesh, const long *part, const long (*cell_map)[2],
                             long (*node_map)[2], long *n_inner_nodes, long *n_nodes);
static void compute_mesh_node_map(Mesh *mesh, const long *part, const long (*node_map)[2],
                                  const long n_nodes);
static void partition_nodes(Mesh *mesh, const long (*node_map)[2]);
static void partition_cells(Mesh *mesh, const long (*node_map)[2], const long (*cell_map)[2]);
static void partition_entities(Mesh *mesh, const long (*cell_map)[2]);
static void compute_sync(Mesh *mesh, const long (*cell_map)[2]);
static int mapcmp(const void *a, const void *b);

void partition(Mesh *mesh, bool reorder)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &mesh->rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mesh->size);
    assert(mesh->size <= mesh->n_inner_cells && "more ranks than cells is not supported");

    cleanup long *cpart = memory_calloc(mesh->n_cells, sizeof(*cpart));
    compute_cell_partitioning(mesh, cpart);

    if (reorder) reorder_cells(mesh, cpart);

    cleanup long(*cell_map)[2] = memory_calloc(mesh->n_cells, sizeof(*cell_map));
    long n_inner_cells, n_ghost_cells, n_cells;
    compute_cell_map(mesh, cpart, cell_map, &n_inner_cells, &n_ghost_cells, &n_cells);

    cleanup long *npart = memory_calloc(mesh->n_nodes, sizeof(*npart));
    compute_node_partitioning(mesh, cpart, npart);

    mesh->n_inner_cells = n_inner_cells;
    mesh->n_ghost_cells = n_ghost_cells;
    mesh->n_cells = n_cells;

    cleanup long(*node_map)[2] = memory_calloc(mesh->n_nodes, sizeof(*node_map));
    long n_inner_nodes, n_nodes;
    compute_node_map(mesh, npart, cell_map, node_map, &n_inner_nodes, &n_nodes);
    compute_mesh_node_map(mesh, npart, node_map, n_nodes);

    mesh->n_inner_nodes = n_inner_nodes;
    mesh->n_nodes = n_nodes;

    partition_nodes(mesh, node_map);
    partition_cells(mesh, node_map, cell_map);
    partition_entities(mesh, cell_map);

    compute_sync(mesh, cell_map);
}

static void compute_cell_partitioning(const Mesh *mesh, long *part)
{
    if (mesh->size == 1) {
        for (long i = 0; i < mesh->n_inner_cells; ++i) part[i] = mesh->rank;
    }
    else if (mesh->size == mesh->n_inner_cells) {
        for (long i = 0; i < mesh->n_inner_cells; ++i) part[i] = i;
    }
    else {
        idx_t n_cells = mesh->n_cells;
        idx_t options[METIS_NOPTIONS], objval;
        METIS_SetDefaultOptions(options);
        options[METIS_OPTION_NUMBERING] = 0;
        METIS_PartGraphKway(&n_cells, (idx_t[]){1}, mesh->cell.i_cell, mesh->cell.cell, 0, 0, 0,
                            (idx_t[]){mesh->size}, 0, 0, options, &objval, part);
    }
    assert(array_contains(part, mesh->n_inner_cells, mesh->rank) && "partitioning error");

    // correct partitioning of ghost cells
    for (long j = mesh->n_inner_cells; j < mesh->n_cells; ++j)
        part[j] = part[mesh->cell.cell[mesh->cell.i_cell[j]]];
}

static void reorder_cells(Mesh *mesh, long *part)
{
    // mark ghost cells as their own partition, separate from all other partitions
    for (long i = mesh->n_inner_cells; i < mesh->n_cells; ++i) part[i] += mesh->size;

    // compute reordering (only for inner cells)
    cleanup long *old2new = memory_calloc(mesh->n_cells, sizeof(*old2new));
    cleanup long *new2old = memory_calloc(mesh->n_cells, sizeof(*new2old));
    METIS_CacheFriendlyReordering(mesh->n_cells, mesh->cell.i_cell, mesh->cell.cell, part, old2new);
    for (long i = mesh->n_inner_cells; i < mesh->n_cells; ++i) old2new[i] = i;
    for (long i = 0; i < mesh->n_cells; ++i) new2old[old2new[i]] = i;

    // undo marking of ghost cells
    for (long i = mesh->n_inner_cells; i < mesh->n_cells; ++i) part[i] -= mesh->size;

    // reorder cell to node
    long *i_node = memory_calloc(mesh->n_cells + 1, sizeof(*i_node));
    long *node = memory_calloc(mesh->cell.i_node[mesh->n_cells], sizeof(*node));
    for (long jnew = 0; jnew < mesh->n_cells; ++jnew) {
        const long jold = new2old[jnew];
        i_node[jnew + 1] = i_node[jnew];
        for (long i = mesh->cell.i_node[jold]; i < mesh->cell.i_node[jold + 1]; ++i)
            node[i_node[jnew + 1]++] = mesh->cell.node[i];
    }
    free(mesh->cell.i_node);
    free(mesh->cell.node);
    mesh->cell.i_node = i_node;
    mesh->cell.node = node;

    // reorder cell to cell
    long *i_cell = memory_calloc(mesh->n_cells + 1, sizeof(*i_cell));
    long *cell = memory_calloc(mesh->cell.i_cell[mesh->n_cells], sizeof(*cell));
    for (long jnew = 0; jnew < mesh->n_cells; ++jnew) {
        const long jold = new2old[jnew];
        i_cell[jnew + 1] = i_cell[jnew];
        for (long i = mesh->cell.i_cell[jold]; i < mesh->cell.i_cell[jold + 1]; ++i)
            cell[i_cell[jnew + 1]++] = old2new[mesh->cell.cell[i]];
    }
    free(mesh->cell.i_cell);
    free(mesh->cell.cell);
    mesh->cell.i_cell = i_cell;
    mesh->cell.cell = cell;

    // reorder part
    cleanup long *new_part = memory_calloc(mesh->n_cells, sizeof(*new_part));
    for (long i = 0; i < mesh->n_cells; ++i) new_part[i] = part[new2old[i]];
    for (long i = 0; i < mesh->n_cells; ++i) part[i] = new_part[i];
}

static void compute_cell_map(const Mesh *mesh, const long *part, long (*cell_map)[2],
                             long *n_inner_cells, long *n_ghost_cells, long *n_cells)
{
    // inner cells
    long n = 0;
    for (long i = 0; i < mesh->n_inner_cells; ++i) {
        if (part[i] == mesh->rank) {
            cell_map[n][0] = i;
            cell_map[n][1] = mesh->rank;
            n += 1;
        }
    }
    *n_inner_cells = n;

    // ghost cells
    for (long e = 0; e < mesh->n_entities; ++e) {
        if (mesh->entity.j_cell[e] < mesh->n_inner_cells) continue;
        for (long j = mesh->entity.j_cell[e]; j < mesh->entity.j_cell[e + 1]; ++j) {
            if (part[j] == mesh->rank) {
                cell_map[n][0] = j;
                cell_map[n][1] = mesh->rank;
                n += 1;
            }
        }
    }
    *n_ghost_cells = n - *n_inner_cells;

    // mpi cells
    cleanup long *seen = memory_calloc(mesh->n_inner_cells, sizeof(*seen));
    for (long jl = 0; jl < *n_inner_cells; ++jl) {
        const long jg = cell_map[jl][0];
        for (long i = mesh->cell.i_cell[jg]; i < mesh->cell.i_cell[jg + 1]; ++i) {
            if (mesh->cell.cell[i] >= mesh->n_inner_cells) continue;
            if (part[mesh->cell.cell[i]] == mesh->rank) continue;
            if (seen[mesh->cell.cell[i]]++) continue;
            cell_map[n][0] = mesh->cell.cell[i];
            cell_map[n][1] = part[mesh->cell.cell[i]];
            n += 1;
        }
    }
    *n_cells = n;

    // sort mpi cells, first by rank, then by global index
    const long n_mpi_cells = *n_cells - *n_inner_cells - *n_ghost_cells;
    qsort(&cell_map[*n_cells - n_mpi_cells], n_mpi_cells, sizeof(*cell_map), mapcmp);
}

static void compute_node_partitioning(const Mesh *mesh, const long *cpart, long *npart)
{
    // count how often each node is part of each partition
    cleanup long(*count)[mesh->size] = memory_calloc(mesh->n_nodes, sizeof(*count));
    for (long j = 0; j < mesh->n_inner_cells; ++j)
        for (long i = mesh->cell.i_node[j]; i < mesh->cell.i_node[j + 1]; ++i)
            count[mesh->cell.node[i]][cpart[j]] += 1;

    // node is part of the rank where it appears most often
    for (long i = 0; i < mesh->n_nodes; ++i) npart[i] = array_argmax(count[i], mesh->size);
}

static void compute_node_map(Mesh *mesh, const long *part, const long (*cell_map)[2],
                             long (*node_map)[2], long *n_inner_nodes, long *n_nodes)
{
    // inner nodes
    long n = 0;
    cleanup long *seen = memory_calloc(mesh->n_nodes, sizeof(*seen));
    for (long i = 0; i < mesh->n_nodes; ++i) {
        if (part[i] == mesh->rank) {
            if (seen[i]++) continue;
            node_map[n][0] = i;
            node_map[n][1] = part[i];
            n += 1;
        }
    }
    *n_inner_nodes = n;

    // mpi nodes
    for (long jl = 0; jl < mesh->n_inner_cells; ++jl) {
        const long jg = cell_map[jl][0];
        for (long i = mesh->cell.i_node[jg]; i < mesh->cell.i_node[jg + 1]; ++i) {
            if (seen[mesh->cell.node[i]]++) continue;
            node_map[n][0] = mesh->cell.node[i];
            node_map[n][1] = part[mesh->cell.node[i]];
            n += 1;
        }
    }
    *n_nodes = n;

    // sort mpi nodes, first by rank, then by global index
    const long n_mpi_nodes = *n_nodes - *n_inner_nodes;
    qsort(node_map[*n_nodes - n_mpi_nodes], n_mpi_nodes, sizeof(*node_map), mapcmp);
}

static void compute_mesh_node_map(Mesh *mesh, const long *part, const long (*node_map)[2],
                                  const long n_nodes)
{
    // count number of nodes per partition and compute global to local node indices
    cleanup long *count = memory_calloc(mesh->size, sizeof(*count));
    cleanup long *g2l = memory_calloc(mesh->n_nodes, sizeof(*g2l));
    for (long i = 0; i < mesh->n_nodes; ++i) g2l[i] = count[part[i]]++;

    // compute partition node offsets
    cleanup long *offset = memory_calloc(mesh->size, sizeof(*offset));
    for (long i = 0; i < mesh->size - 1; ++i) offset[i + 1] = offset[i] + count[i];

    // compute local to global node map
    mesh->node.map = memory_calloc(n_nodes, sizeof(*mesh->node.map));
    for (long il = 0; il < n_nodes; ++il) {
        const long ig = node_map[il][0];
        const long rank = node_map[il][1];
        mesh->node.map[il] = g2l[ig] + offset[rank];
    }
}

static void partition_nodes(Mesh *mesh, const long (*node_map)[2])
{
    double(*x)[N_DIMS] = memory_calloc(mesh->n_nodes, sizeof(*x));
    for (long jl = 0; jl < mesh->n_nodes; ++jl) {
        const long jg = node_map[jl][0];
        for (long d = 0; d < N_DIMS; ++d) x[jl][d] = mesh->node.x[jg][d];
    }
    free(mesh->node.x);
    mesh->node.x = x;
}

static void partition_cells(Mesh *mesh, const long (*node_map)[2], const long (*cell_map)[2])
{
    // compute global to local node map
    fcleanup(dict_free) Dict g2ln = dict_create(mesh->n_nodes);
    for (long i = 0; i < mesh->n_nodes; ++i) dict_insert(&g2ln, &node_map[i][0], 1, &i, 1);

    // compute new cell to node map
    long *i_node = memory_calloc(mesh->n_cells + 1, sizeof(*i_node));
    long *node = memory_calloc(mesh->n_cells * MAX_CELL_NODES, sizeof(*node));
    for (long jl = 0; jl < mesh->n_cells; ++jl) {
        const long jg = cell_map[jl][0];
        const long n_nodes = mesh->cell.i_node[jg + 1] - mesh->cell.i_node[jg];
        i_node[jl + 1] = i_node[jl];
        for (long i = 0; i < n_nodes; ++i) {
            long *n;
            if (dict_lookup(&g2ln, &mesh->cell.node[mesh->cell.i_node[jg] + i], 1, &n))
                node[i_node[jl + 1]++] = *n;
        }
    }
    free(mesh->cell.i_node);
    free(mesh->cell.node);
    mesh->cell.i_node = i_node;
    mesh->cell.node = memory_realloc(node, i_node[mesh->n_cells], sizeof(*node));

    // compute global to local cell map
    fcleanup(dict_free) Dict g2lc = dict_create(mesh->n_cells);
    for (long i = 0; i < mesh->n_cells; ++i) dict_insert(&g2lc, &cell_map[i][0], 1, &i, 1);

    // compute new cell to cell map
    long *i_cell = memory_calloc(mesh->n_cells + 1, sizeof(*i_cell));
    long *cell = memory_calloc(mesh->n_cells * MAX_CELL_FACES, sizeof(*cell));
    for (long jl = 0; jl < mesh->n_cells; ++jl) {
        const long jg = cell_map[jl][0];
        const long n_cells = mesh->cell.i_cell[jg + 1] - mesh->cell.i_cell[jg];
        i_cell[jl + 1] = i_cell[jl];
        for (long i = 0; i < n_cells; ++i) {
            long *c;
            if (!dict_lookup(&g2lc, &mesh->cell.cell[mesh->cell.i_cell[jg] + i], 1, &c)) continue;
            if (jl >= mesh->n_inner_cells + mesh->n_ghost_cells && *c >= mesh->n_inner_cells)
                continue;  // filter out connectivity between mpi cells
            cell[i_cell[jl + 1]++] = *c;
        }
    }
    free(mesh->cell.i_cell);
    free(mesh->cell.cell);
    mesh->cell.i_cell = i_cell;
    mesh->cell.cell = memory_realloc(cell, i_cell[mesh->n_cells], sizeof(*cell));
}

static void partition_entities(Mesh *mesh, const long (*cell_map)[2])
{
    // compute entity cell offsets
    cleanup long *n_cells = memory_calloc(mesh->n_entities, sizeof(*n_cells));
    for (long i = 0; i < mesh->n_cells; ++i) {
        if (cell_map[i][1] != mesh->rank) continue;
        for (long e = 0; e < mesh->n_entities; ++e) {
            if (mesh->entity.j_cell[e] <= cell_map[i][0] &&
                cell_map[i][0] < mesh->entity.j_cell[e + 1])
                n_cells[e] += 1;
        }
    }

    // update entity cells
    for (long e = 0; e < mesh->n_entities; ++e)
        mesh->entity.j_cell[e + 1] = mesh->entity.j_cell[e] + n_cells[e];
}

static void compute_sync(Mesh *mesh, const long (*cell_map)[2])
{
    // compute receive map and find all cells that need to be sent to other ranks
    mesh->sync.i_recv = memory_calloc(mesh->size + 1, sizeof(*mesh->sync.i_recv));
    const long n_mpi_cells = mesh->n_cells - mesh->n_inner_cells - mesh->n_ghost_cells;
    fcleanup(dict_free) Dict sync_dict = dict_create(mesh->size * n_mpi_cells);
    for (long j = mesh->n_inner_cells + mesh->n_ghost_cells; j < mesh->n_cells; ++j) {
        const long rank = cell_map[j][1];
        mesh->sync.i_recv[rank + 1] += 1;
        for (long i = mesh->cell.i_cell[j]; i < mesh->cell.i_cell[j + 1]; ++i)
            dict_insert(&sync_dict, (long[]){rank, mesh->cell.cell[i]}, 2, 0, 0);
    }
    mesh->sync.i_recv[0] = mesh->n_inner_cells + mesh->n_ghost_cells;
    for (long i = 0; i < mesh->size; ++i) mesh->sync.i_recv[i + 1] += mesh->sync.i_recv[i];
    cleanup DictItem *sync = dict_serialize(&sync_dict);

    // compute send map
    mesh->sync.i_send = memory_calloc(mesh->size + 1, sizeof(*mesh->sync.i_send));
    mesh->sync.send = memory_calloc(sync_dict.n_items, sizeof(*mesh->sync.send));
    for (long i = 0; i < sync_dict.n_items; ++i) {
        const long rank = sync[i].key[0];
        const long cell = sync[i].key[1];
        mesh->sync.i_send[rank + 1] += 1;
        mesh->sync.send[i] = cell;
    }
    for (long i = 0; i < mesh->size; ++i) mesh->sync.i_send[i + 1] += mesh->sync.i_send[i];
}

static int mapcmp(const void *a, const void *b)
{
    const long *map_a = a;
    const long *map_b = b;

    // sort by rank
    if (map_a[1] < map_b[1]) return -1;
    if (map_a[1] > map_b[1]) return 1;

    // sort by global index
    if (map_a[0] < map_b[0]) return -1;
    if (map_a[0] > map_b[0]) return 1;

    return 0;
}
