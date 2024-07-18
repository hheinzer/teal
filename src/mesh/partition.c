#include "partition.h"

#include <metis.h>
#include <stdlib.h>
#include <string.h>

#include "core/array.h"
#include "core/dict.h"
#include "core/memory.h"
#include "core/utils.h"
#include "reorder.h"
#include "sendrecv.h"
#include "teal.h"

static void compute_cell_partition(Mesh *mesh, long *cell_part);
static void compute_node_partition(const Mesh *mesh, const long *cell_part, long *node_part);
static void compute_cell_map(const Mesh *mesh, const long *cell_part, long (*cell_map)[2],
                             Mesh *part, long rank);
static void compute_node_map(const Mesh *mesh, const long *node_part, const long (*cell_map)[2],
                             long (*node_map)[2], Mesh *part, long rank);
static void compute_node_idx(const Mesh *mesh, const long *node_part, const long (*node_map)[2],
                             Mesh *part);
static void partition_nodes(const Mesh *mesh, const long (*node_map)[2], Mesh *part);
static void partition_cells(const Mesh *mesh, const long (*node_map)[2], const long (*cell_map)[2],
                            Mesh *part);
static void partition_entities(const Mesh *mesh, const long *cell_part, Mesh *part, long rank);
static void compute_sync(Mesh *part, const long (*cell_map)[2]);
static int mapcmp(const void *a, const void *b);

void partition(Mesh *mesh)
{
    if (teal.rank == 0) {
        if (teal.size > mesh->n_inner_cells)
            error("number of ranks '%d' cannot be larger than number of inner cells '%ld'",
                  teal.size, mesh->n_inner_cells);

        cleanup long *cell_part = memory_calloc(mesh->n_cells, sizeof(*cell_part));
        cleanup long *node_part = memory_calloc(mesh->n_nodes, sizeof(*node_part));
        compute_cell_partition(mesh, cell_part);
        compute_node_partition(mesh, cell_part, node_part);

        Mesh part = {0};
        for (long rank = teal.size - 1; rank >= 0; --rank) {
            mesh_free(&part);

            cleanup long(*cell_map)[2] = memory_calloc(mesh->n_cells, sizeof(*cell_map));
            cleanup long(*node_map)[2] = memory_calloc(mesh->n_nodes, sizeof(*node_map));
            compute_cell_map(mesh, cell_part, cell_map, &part, rank);
            compute_node_map(mesh, node_part, cell_map, node_map, &part, rank);
            compute_node_idx(mesh, node_part, node_map, &part);

            partition_nodes(mesh, node_map, &part);
            partition_cells(mesh, node_map, cell_map, &part);
            partition_entities(mesh, cell_part, &part, rank);
            compute_sync(&part, cell_map);

            if (rank > 0) send_mesh(&part, rank);
        }

        mesh_free(mesh);
        *mesh = part;
    }
    else
        recv_mesh(mesh, 0);
}

static void compute_cell_partition(Mesh *mesh, long *cell_part)
{
    const ALIAS(i_cell, mesh->cell.i_cell);
    const ALIAS(cell, mesh->cell.cell);
    idx_t nvtxs = mesh->n_inner_cells;
    idx_t ncon = 1;
    cleanup idx_t *xadj = memory_calloc(nvtxs + 1, sizeof(*xadj));
    cleanup idx_t *adjncy = memory_calloc(nvtxs * MAX_CELL_FACES, sizeof(*adjncy));
    for (long j = 0; j < nvtxs; ++j) {
        xadj[j + 1] = xadj[j];
        for (long i = i_cell[j]; i < i_cell[j + 1]; ++i) {
            if (cell[i] >= mesh->n_inner_cells) continue;  // remove non-inner connections
            adjncy[xadj[j + 1]++] = cell[i];
        }
    }
    idx_t nparts = teal.size;
    idx_t edgecut;
    if (nparts > 8)
        METIS_PartGraphKway(&nvtxs, &ncon, xadj, adjncy, 0, 0, 0, &nparts, 0, 0, 0, &edgecut,
                            cell_part);
    else if (nparts > 1)
        METIS_PartGraphRecursive(&nvtxs, &ncon, xadj, adjncy, 0, 0, 0, &nparts, 0, 0, 0, &edgecut,
                                 cell_part);

    for (long rank = 0; rank < teal.size; ++rank)
        if (!array_contains(cell_part, mesh->n_inner_cells, rank))
            error("mesh contains no cells for rank '%ld'", rank);

    for (long j = mesh->n_inner_cells; j < mesh->n_cells; ++j)
        cell_part[j] = cell_part[cell[i_cell[j]]];  // assign ghost cells to partition of inner cell

    cleanup long *old2new = memory_calloc(mesh->n_cells, sizeof(*old2new));
    cleanup long *new2old = memory_calloc(mesh->n_cells, sizeof(*new2old));
    METIS_CacheFriendlyReordering(nvtxs, xadj, adjncy, cell_part, old2new);
    for (long j = mesh->n_inner_cells; j < mesh->n_cells; ++j) old2new[j] = j;
    for (long j = 0; j < mesh->n_cells; ++j) new2old[old2new[j]] = j;
    reorder_cells(mesh, old2new, new2old);

    cleanup long *part = memory_calloc(mesh->n_cells, sizeof(*part));
    for (long j = 0; j < mesh->n_cells; ++j) part[j] = cell_part[new2old[j]];
    memory_copy(cell_part, part, mesh->n_cells, sizeof(*part));
}

static void compute_node_partition(const Mesh *mesh, const long *cell_part, long *node_part)
{
    const ALIAS(i_node, mesh->cell.i_node);
    const ALIAS(node, mesh->cell.node);
    cleanup long(*count)[teal.size] = memory_calloc(mesh->n_nodes, sizeof(*count));
    for (long j = 0; j < mesh->n_inner_cells; ++j)
        for (long i = i_node[j]; i < i_node[j + 1]; ++i) count[node[i]][cell_part[j]] += 1;
    for (long i = 0; i < mesh->n_nodes; ++i) node_part[i] = array_argmax(count[i], teal.size);
}

static void compute_cell_map(const Mesh *mesh, const long *cell_part, long (*cell_map)[2],
                             Mesh *part, long rank)
{
    long n = 0;
    for (long j = 0; j < mesh->n_inner_cells; ++j) {
        if (cell_part[j] == rank) {
            cell_map[n][0] = j;
            cell_map[n][1] = cell_part[j];
            n += 1;
        }
    }
    part->n_inner_cells = n;

    for (long j = mesh->n_inner_cells; j < mesh->n_inner_cells + mesh->n_ghost_cells; ++j) {
        if (cell_part[j] == rank) {
            cell_map[n][0] = j;
            cell_map[n][1] = cell_part[j];
            n += 1;
        }
    }
    part->n_ghost_cells = n - part->n_inner_cells;

    const ALIAS(i_cell, mesh->cell.i_cell);
    const ALIAS(cell, mesh->cell.cell);
    cleanup long *seen = memory_calloc(mesh->n_inner_cells, sizeof(*seen));
    for (long jl = 0; jl < part->n_inner_cells; ++jl) {
        const long jg = cell_map[jl][0];
        for (long i = i_cell[jg]; i < i_cell[jg + 1]; ++i) {
            if (cell[i] >= mesh->n_inner_cells) continue;  // skip ghost cells
            if (cell_part[cell[i]] == rank) continue;      // skip local cells
            if (seen[cell[i]]++) continue;                 // skip duplicates
            cell_map[n][0] = cell[i];
            cell_map[n][1] = cell_part[cell[i]];
            n += 1;
        }
    }
    part->n_sync_cells = n - part->n_inner_cells - part->n_ghost_cells;
    part->n_cells = n;
    qsort(cell_map[n - part->n_sync_cells], part->n_sync_cells, sizeof(*cell_map), mapcmp);
}

static void compute_node_map(const Mesh *mesh, const long *node_part, const long (*cell_map)[2],
                             long (*node_map)[2], Mesh *part, long rank)
{
    long n = 0;
    cleanup long *seen = memory_calloc(mesh->n_nodes, sizeof(*seen));
    for (long i = 0; i < mesh->n_nodes; ++i) {
        if (node_part[i] == rank) {
            if (seen[i]++) continue;
            node_map[n][0] = i;
            node_map[n][1] = node_part[i];
            n += 1;
        }
    }
    part->n_inner_nodes = n;

    const ALIAS(i_node, mesh->cell.i_node);
    const ALIAS(node, mesh->cell.node);
    for (long jl = 0; jl < mesh->n_inner_cells; ++jl) {
        const long jg = cell_map[jl][0];
        for (long i = i_node[jg]; i < i_node[jg + 1]; ++i) {
            if (seen[node[i]]++) continue;
            node_map[n][0] = node[i];
            node_map[n][1] = node_part[node[i]];
            n += 1;
        }
    }
    part->n_sync_nodes = n - part->n_inner_nodes;
    part->n_nodes = n;
    qsort(node_map[n - part->n_sync_nodes], part->n_sync_nodes, sizeof(*node_map), mapcmp);
}

static void compute_node_idx(const Mesh *mesh, const long *node_part, const long (*node_map)[2],
                             Mesh *part)
{
    long *idx = memory_calloc(part->n_nodes, sizeof(*idx));
    cleanup long *count = memory_calloc(teal.size, sizeof(*count));
    cleanup long *g2l = memory_calloc(mesh->n_nodes, sizeof(*g2l));
    cleanup long *offset = memory_calloc(teal.size, sizeof(*offset));
    for (long i = 0; i < mesh->n_nodes; ++i) g2l[i] = count[node_part[i]]++;
    for (long i = 0; i < teal.size - 1; ++i) offset[i + 1] = offset[i] + count[i];
    for (long il = 0; il < part->n_nodes; ++il) {
        const long ig = node_map[il][0];
        const long rank = node_map[il][1];
        idx[il] = g2l[ig] + offset[rank];
    }
    part->node.idx = idx;
}

static void partition_nodes(const Mesh *mesh, const long (*node_map)[2], Mesh *part)
{
    double(*coord)[N_DIMS] = memory_calloc(part->n_nodes, sizeof(*coord));
    for (long jl = 0; jl < part->n_nodes; ++jl) {
        const long jg = node_map[jl][0];
        for (long d = 0; d < N_DIMS; ++d) coord[jl][d] = mesh->node.coord[jg][d];
    }
    part->node.coord = coord;
}

static void partition_cells(const Mesh *mesh, const long (*node_map)[2], const long (*cell_map)[2],
                            Mesh *part)
{
    cleanup long *node_idx = memory_calloc(mesh->n_nodes, sizeof(*node_idx));
    cleanup long *cell_idx = memory_calloc(mesh->n_cells, sizeof(*cell_idx));
    for (long i = 0; i < mesh->n_nodes; ++i) node_idx[i] = -1;
    for (long i = 0; i < mesh->n_cells; ++i) cell_idx[i] = -1;
    for (long i = 0; i < part->n_nodes; ++i) node_idx[node_map[i][0]] = i;
    for (long i = 0; i < part->n_cells; ++i) cell_idx[cell_map[i][0]] = i;

    long *i_node = memory_calloc(part->n_cells + 1, sizeof(*i_node));
    long *i_cell = memory_calloc(part->n_cells + 1, sizeof(*i_cell));
    long *node = memory_calloc(part->n_cells * MAX_CELL_NODES, sizeof(*node));
    long *cell = memory_calloc(part->n_cells * MAX_CELL_FACES, sizeof(*cell));
    for (long jl = 0; jl < part->n_cells; ++jl) {
        const long jg = cell_map[jl][0];

        i_node[jl + 1] = i_node[jl];
        for (long i = mesh->cell.i_node[jg]; i < mesh->cell.i_node[jg + 1]; ++i) {
            const long idx = node_idx[mesh->cell.node[i]];
            if (idx >= 0) node[i_node[jl + 1]++] = idx;
        }

        i_cell[jl + 1] = i_cell[jl];
        for (long i = mesh->cell.i_cell[jg]; i < mesh->cell.i_cell[jg + 1]; ++i) {
            const long idx = cell_idx[mesh->cell.cell[i]];
            if (idx < 0) continue;
            if (jl >= part->n_inner_cells + part->n_ghost_cells && idx >= part->n_inner_cells)
                continue;  // filter out connectivity between sync cells
            cell[i_cell[jl + 1]++] = idx;
        }
    }
    part->cell.i_node = i_node;
    part->cell.i_cell = i_cell;
    part->cell.node = memory_realloc(node, i_node[part->n_cells], sizeof(*node));
    part->cell.cell = memory_realloc(cell, i_cell[part->n_cells], sizeof(*cell));
}

static void partition_entities(const Mesh *mesh, const long *cell_part, Mesh *part, long rank)
{
    char(*name)[NAMELEN] = memory_calloc(mesh->n_entities, sizeof(*name));
    long *j_cell = memory_calloc(mesh->n_entities + 1, sizeof(*j_cell));
    double(*offset)[N_DIMS] = memory_calloc(mesh->n_entities, sizeof(*offset));
    for (long e = 0; e < mesh->n_entities; ++e) {
        strcpy(name[e], mesh->entity.name[e]);

        j_cell[e + 1] = j_cell[e];
        for (long j = mesh->entity.j_cell[e]; j < mesh->entity.j_cell[e + 1]; ++j)
            if (cell_part[j] == rank) j_cell[e + 1] += 1;

        for (long d = 0; d < N_DIMS; ++d) offset[e][d] = mesh->entity.offset[e][d];
    }
    part->n_entities = mesh->n_entities;
    part->entity.name = name;
    part->entity.j_cell = j_cell;
    part->entity.offset = offset;
}

static void compute_sync(Mesh *part, const long (*cell_map)[2])
{
    long *j_recv = memory_calloc(teal.size + 1, sizeof(*j_recv));
    fcleanup(dict_free) Dict sync = dict_create(teal.size * part->n_sync_cells);
    for (long j = part->n_inner_cells + part->n_ghost_cells; j < part->n_cells; ++j) {
        const long rank = cell_map[j][1];
        j_recv[rank + 1] += 1;
        for (long i = part->cell.i_cell[j]; i < part->cell.i_cell[j + 1]; ++i)
            dict_insert(&sync, (long[]){rank, part->cell.cell[i]}, 2, 0, 0);
    }
    j_recv[0] = part->n_inner_cells + part->n_ghost_cells;
    for (long i = 0; i < teal.size; ++i) j_recv[i + 1] += j_recv[i];
    part->sync.j_recv = j_recv;

    long *i_send = memory_calloc(teal.size + 1, sizeof(*i_send));
    long *send = memory_calloc(sync.n_items, sizeof(*send));
    cleanup DictItem *item = dict_serialize_by_key(&sync);
    for (long i = 0; i < sync.n_items; ++i) {
        const long rank = item[i].key[0];
        const long cell = item[i].key[1];
        i_send[rank + 1] += 1;
        send[i] = cell;
    }
    for (long i = 0; i < teal.size; ++i) i_send[i + 1] += i_send[i];
    part->sync.i_send = i_send;
    part->sync.send = send;
}

static int mapcmp(const void *a, const void *b)
{
    const long *aa = a;
    const long *bb = b;
    if (aa[1] < bb[1]) return -1;
    if (aa[1] > bb[1]) return 1;
    if (aa[0] < bb[0]) return -1;
    if (aa[0] > bb[0]) return 1;
    return 0;
}
