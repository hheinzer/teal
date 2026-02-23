#include "dual.h"

#include <assert.h>
#include <limits.h>

#include "sync2.h"
#include "teal2.h"
#include "utils2.h"

Dual *dual_init(const Grid *grid)
{
    Dual *dual = teal2_calloc(1, sizeof(*dual));

    dual->dist = teal2_calloc(sync2.size + 1, sizeof(*dual->dist));
    sync2_offsets(&(long){grid->cells.off_periodic}, dual->dist, 1, MPI_LONG);

    long numflag = 0;
    long ncommon = 3;
    int ret = ParMETIS_V3_Mesh2Dual(dual->dist, grid->cells.node_off, grid->cells.node_idx,
                                    &numflag, &ncommon, &dual->xadj, &dual->adjncy, &sync2.comm);
    if (ret != METIS_OK) {
        teal2_error("ParMETIS_V3_Mesh2Dual failed (%d)", ret);
    }

    return dual;
}

void dual_deinit(Dual *dual)
{
    teal2_free(dual->dist);
    teal2_free(dual->xadj);
    teal2_free(dual->adjncy);
    teal2_free(dual->part);
    teal2_free(dual);
}

static void collect_outer(long *part, const Dual *dual, const Grid *grid)
{
    int num_cells = grid->cells.num - grid->cells.num_inner;
    long *adjncy = teal2_calloc(num_cells, sizeof(*adjncy));
    int num = 0;
    for (int i = grid->cells.num_inner; i < grid->cells.num; i++) {
        assert(dual->xadj[i + 1] - dual->xadj[i] == 1);
        adjncy[num++] = dual->adjncy[dual->xadj[i]];
    }
    assert(num == num_cells);

    long *recv = teal2_calloc(num_cells, sizeof(*recv));
    sync2_collect(part, recv, adjncy, grid->cells.off_periodic, num_cells, MPI_LONG);

    num = 0;
    for (int i = grid->cells.num_inner; i < grid->cells.num; i++) {
        part[i] = recv[num++];
    }
    assert(num == num_cells);

    teal2_free(adjncy);
    teal2_free(recv);
}

void dual_partition(Dual *dual, const Grid *grid)
{
    long *vwgt = teal2_calloc(grid->cells.off_periodic, sizeof(*vwgt));
    for (int i = 0; i < grid->cells.num_inner; i++) {
        vwgt[i] = 1;
    }

    long wgtflag = 2;
    long numflag = 0;
    long ncon = 1;
    long nparts = sync2.size;

    assert(ncon <= INT_MAX / nparts);
    double *tpwgts = teal2_calloc((int)(ncon * nparts), sizeof(*tpwgts));
    for (int i = 0; i < (int)(ncon * nparts); i++) {
        tpwgts[i] = 1 / (double)nparts;
    }

    double *ubvec = teal2_calloc((int)ncon, sizeof(*ubvec));
    for (int i = 0; i < ncon; i++) {
        ubvec[i] = 1.05;
    }

    long *part = teal2_calloc(grid->cells.off_periodic, sizeof(*part));

    long options[3] = {0};
    long edgecut;
    int ret =
        ParMETIS_V3_PartKway(dual->dist, dual->xadj, dual->adjncy, vwgt, 0, &wgtflag, &numflag,
                             &ncon, &nparts, tpwgts, ubvec, options, &edgecut, part, &sync2.comm);
    if (ret != METIS_OK) {
        teal2_error("ParMETIS_V3_PartKway failed (%d)", ret);
    }

    if (sync2.rank == 0) {
        teal2_verbose("ParMETIS_V3_PartKway edgecut = %ld", edgecut);
    }

    long *num_inner = teal2_calloc(sync2.size, sizeof(*num_inner));
    for (int i = 0; i < grid->cells.num_inner; i++) {
        num_inner[part[i]] += 1;
    }
    sync2_sum(num_inner, sync2.size, MPI_LONG);

    if (num_inner[sync2.rank] <= 0 || INT_MAX < num_inner[sync2.rank]) {
        teal2_error("invalid partition size (%ld)", num_inner[sync2.rank]);
    }

    // this should not be necessary since outer cells have a weight of 0, and should therefore be
    // assigned to same partition as their connected inner cells; do it anyways to be safe
    collect_outer(part, dual, grid);

    dual->part = part;

    teal2_free(vwgt);
    teal2_free(tpwgts);
    teal2_free(ubvec);
    teal2_free(num_inner);
}

typedef struct {
    int entity;
    long node;
    long peer;
} Map;

static MPI_Datatype datatype_map(void)
{
    MPI_Datatype datatype;
    int len[3] = {1, 1, 1};
    MPI_Aint off[3] = {offsetof(Map, entity), offsetof(Map, node), offsetof(Map, peer)};
    MPI_Datatype type[3] = {MPI_INT, MPI_LONG, MPI_LONG};
    MPI_Type_create_struct(3, len, off, type, &datatype);
    return sync2_resized(datatype, sizeof(Map));
}

static int compare_map(const void *lhs_, const void *rhs_)
{
    const Map *lhs = lhs_;
    const Map *rhs = rhs_;
    if (lhs->entity != rhs->entity) {
        return (lhs->entity < rhs->entity) ? -1 : +1;
    }
    return (lhs->node > rhs->node) - (lhs->node < rhs->node);
}

static Map *compute_node_map(const Grid *grid, int *num_nodes)
{
    long max_nodes = grid->cells.node_off[grid->cells.off_periodic] -
                     grid->cells.node_off[grid->cells.off_boundary];

    assert(max_nodes <= INT_MAX);
    Map *map = teal2_calloc((int)max_nodes, sizeof(*map));
    int num = 0;
    for (int i = grid->entities.off_boundary; i < grid->entities.num; i++) {
        for (int j = grid->entities.cell_off[i]; j < grid->entities.cell_off[i + 1]; j++) {
            for (long k = grid->cells.node_off[j]; k < grid->cells.node_off[j + 1]; k++) {
                map[num].entity = i;
                map[num].node = grid->cells.node_idx[k];
                map[num].peer = -1;
                num += 1;
            }
        }
    }
    *num_nodes = unique(map, num, sizeof(*map), compare_map);

    assert(grid->entities.node_off[grid->entities.num] <= INT_MAX);
    int num_periodic = (int)grid->entities.node_off[grid->entities.num];
    int cap = num_periodic;
    sync2_max(&cap, 1, MPI_INT);

    Map *periodic = teal2_calloc(cap, sizeof(*periodic));
    num = 0;
    for (int i = grid->entities.off_boundary; i < grid->entities.num; i++) {
        assert(grid->entities.periodic[i] >= grid->entities.off_boundary);
        long *off_node = &grid->entities.node_off[i];
        long *off_peer = &grid->entities.node_off[grid->entities.periodic[i]];
        assert(off_node[1] - off_node[0] == off_peer[1] - off_peer[0]);
        for (int j = 0; j < off_node[1] - off_node[0]; j++) {
            periodic[num].entity = i;
            periodic[num].node = grid->entities.node_idx[off_node[0] + j];
            periodic[num].peer = grid->entities.node_idx[off_peer[0] + j];
            num += 1;
        }
    }
    assert(num == num_periodic);
    sort(periodic, num, sizeof(*periodic), compare_map);

    MPI_Datatype type = datatype_map();
    for (int rank = 0; rank < sync2.size; rank++) {
        for (int i = 0; i < *num_nodes; i++) {
            if (map[i].peer == -1) {
                Map key = map[i];
                Map *val = search(&key, periodic, num, sizeof(*periodic), compare_map);
                if (val) {
                    map[i].peer = val->peer;
                }
            }
        }
        sync2_rotate(periodic, &num, cap, type);
    }
    assert(num == num_periodic);
    MPI_Type_free(&type);

    for (int i = 0; i < *num_nodes; i++) {
        if (map[i].peer == -1) {
            teal2_error("missing periodic peer node (%d, %" PRIi64 ")", map[i].entity, map[i].node);
        }
    }

    teal2_free(periodic);
    return map;
}

typedef struct {
    int entity;
    int num_nodes;
    long node[MAX_CELL_NODES];
    long inner;
    int idx;
} Cell;

static MPI_Datatype datatype_cell(void)
{
    MPI_Datatype datatype;
    int len[5] = {1, 1, MAX_CELL_NODES, 1, 1};
    MPI_Aint off[5] = {offsetof(Cell, entity), offsetof(Cell, num_nodes), offsetof(Cell, node),
                       offsetof(Cell, inner), offsetof(Cell, idx)};
    MPI_Datatype type[5] = {MPI_INT, MPI_INT, MPI_LONG, MPI_LONG, MPI_INT};
    MPI_Type_create_struct(5, len, off, type, &datatype);
    return sync2_resized(datatype, sizeof(Cell));
}

static int compare_cell(const void *lhs_, const void *rhs_)
{
    const Cell *lhs = lhs_;
    const Cell *rhs = rhs_;
    if (lhs->entity != rhs->entity) {
        return (lhs->entity < rhs->entity) ? -1 : +1;
    }
    if (lhs->num_nodes != rhs->num_nodes) {
        return (lhs->num_nodes < rhs->num_nodes) ? -1 : +1;
    }
    for (int i = 0; i < lhs->num_nodes; i++) {
        if (lhs->node[i] != rhs->node[i]) {
            return (lhs->node[i] < rhs->node[i]) ? -1 : +1;
        }
    }
    return 0;
}

static long *collect_edges(const Dual *dual, const Grid *grid, int num_edges)
{
    int num_nodes;
    Map *map = compute_node_map(grid, &num_nodes);

    Cell *cell = teal2_calloc(num_edges, sizeof(*cell));
    int num = 0;
    for (int i = grid->entities.off_boundary; i < grid->entities.num; i++) {
        for (int j = grid->entities.cell_off[i]; j < grid->entities.cell_off[i + 1]; j++) {
            cell[num].entity = i;
            cell[num].num_nodes = 0;
            for (long k = grid->cells.node_off[j]; k < grid->cells.node_off[j + 1]; k++) {
                assert(cell[num].num_nodes < MAX_CELL_NODES);
                cell[num].node[cell[num].num_nodes++] = grid->cells.node_idx[k];
            }
            sort(cell[num].node, cell[num].num_nodes, sizeof(*cell[num].node), compare_idx);
            assert(dual->xadj[j + 1] - dual->xadj[j] == 1);
            cell[num].inner = dual->adjncy[dual->xadj[j]];
            cell[num].idx = num;
            num += 1;
        }
    }
    assert(num == num_edges);
    sort(cell, num_edges, sizeof(*cell), compare_cell);

    int cap = num_edges;
    sync2_max(&cap, 1, MPI_INT);

    Cell *peer = teal2_calloc(cap, sizeof(*peer));
    for (int i = 0; i < num_edges; i++) {
        peer[i].entity = grid->entities.periodic[cell[i].entity];
        peer[i].num_nodes = cell[i].num_nodes;
        for (int j = 0; j < cell[i].num_nodes; j++) {
            Map key = {.entity = cell[i].entity, .node = cell[i].node[j]};
            Map *val = search(&key, map, num_nodes, sizeof(*map), compare_map);
            assert(val);
            peer[i].node[j] = val->peer;
        }
        sort(peer[i].node, peer[i].num_nodes, sizeof(*peer[i].node), compare_idx);
        peer[i].inner = -1;
    }

    MPI_Datatype type = datatype_cell();
    for (int rank = 0; rank < sync2.size; rank++) {
        for (int i = 0; i < num; i++) {
            if (peer[i].inner == -1) {
                Cell key = peer[i];
                Cell *val = search(&key, cell, num_edges, sizeof(*cell), compare_cell);
                if (val) {
                    peer[i].inner = val->inner;
                }
            }
        }
        sync2_rotate(peer, &num, cap, type);
    }
    assert(num == num_edges);
    MPI_Type_free(&type);

    long *edge = teal2_calloc(num_edges, sizeof(*edge));
    for (int i = 0; i < num_edges; i++) {
        if (peer[i].inner == -1) {
            teal2_error("missing periodic peer cell (%d, %ld)", peer[i].entity, cell[i].inner);
        }
        edge[cell[i].idx] = peer[i].inner;
    }

    teal2_free(map);
    teal2_free(cell);
    teal2_free(peer);
    return edge;
}

void dual_periodic(Dual *dual, const Grid *grid)
{
    int num_edges = grid->cells.off_periodic - grid->cells.off_boundary;
    long *edge = collect_edges(dual, grid, num_edges);

    assert(dual->xadj[grid->cells.off_periodic] <= INT_MAX - num_edges);
    int num_adjncys = (int)dual->xadj[grid->cells.off_periodic] + num_edges;

    long *xadj = teal2_calloc(grid->cells.off_periodic + 1, sizeof(*xadj));
    long *adjncy = teal2_calloc(num_adjncys, sizeof(*adjncy));

    int num = 0;
    for (int i = 0; i < grid->cells.off_periodic; i++) {
        xadj[i + 1] = xadj[i];
        for (long j = dual->xadj[i]; j < dual->xadj[i + 1]; j++) {
            adjncy[xadj[i + 1]++] = dual->adjncy[j];
        }
        if (i >= grid->cells.off_boundary) {
            adjncy[xadj[i + 1]++] = edge[num++];
        }
    }
    assert(num == num_edges && xadj[grid->cells.off_periodic] == num_adjncys);

    teal2_free(dual->xadj);
    teal2_free(dual->adjncy);

    dual->xadj = xadj;
    dual->adjncy = adjncy;

    teal2_free(edge);
}

int compare_idx(const void *lhs_, const void *rhs_)
{
    const long *lhs = lhs_;
    const long *rhs = rhs_;
    return (*lhs > *rhs) - (*lhs < *rhs);
}
