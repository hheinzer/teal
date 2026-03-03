#include <assert.h>
#include <limits.h>

#include "dual.h"
#include "grid.h"
#include "mesh2.h"
#include "sync2.h"
#include "teal2.h"
#include "utils2.h"

static long *collect_adjncys(const Dual *dual, const long *part, const Grid *grid, int *num_adjncys)
{
    int *num_send = teal2_calloc(sync2.size, sizeof(*num_send));
    for (int i = 0; i < grid->cells.off_periodic; i++) {
        long num_cells = dual->xadj[i + 1] - dual->xadj[i];
        assert(num_cells <= INT_MAX);
        num_send[part[i]] += (int)num_cells;
    }

    int *off_send = teal2_calloc(sync2.size + 1, sizeof(*off_send));
    for (int i = 0; i < sync2.size; i++) {
        off_send[i + 1] = off_send[i] + num_send[i];
    }

    int *cur_send = teal2_calloc(sync2.size, sizeof(*cur_send));
    copy(cur_send, off_send, sync2.size, sizeof(*cur_send));

    int tot_send = off_send[sync2.size];
    long *send = teal2_calloc(tot_send, sizeof(*send));
    for (int i = 0; i < grid->cells.off_periodic; i++) {
        for (long j = dual->xadj[i]; j < dual->xadj[i + 1]; j++) {
            send[cur_send[part[i]]++] = dual->adjncy[j];
        }
    }
    for (int i = 0; i < sync2.size; i++) {
        assert(cur_send[i] == off_send[i + 1]);
    }

    int *num_recv = teal2_calloc(sync2.size, sizeof(*num_recv));
    MPI_Alltoall(num_send, 1, MPI_INT, num_recv, 1, MPI_INT, sync2.comm);

    int *off_recv = teal2_calloc(sync2.size + 1, sizeof(*off_recv));
    for (int i = 0; i < sync2.size; i++) {
        off_recv[i + 1] = off_recv[i] + num_recv[i];
    }

    int tot_recv = off_recv[sync2.size];
    long *recv = teal2_calloc(tot_recv, sizeof(*recv));
    MPI_Alltoallv(send, num_send, off_send, MPI_LONG, recv, num_recv, off_recv, MPI_LONG,
                  sync2.comm);

    *num_adjncys = unique(recv, tot_recv, sizeof(*recv), compare_long);

    teal2_free(num_send);
    teal2_free(off_send);
    teal2_free(cur_send);
    teal2_free(send);
    teal2_free(num_recv);
    teal2_free(off_recv);
    return recv;
}

static long *collect_parts(const long *adjncy, const Dual *dual, const long *part, int num_adjncys)
{
    long num_send = dual->dist[sync2.rank + 1] - dual->dist[sync2.rank];
    assert(0 <= num_send && num_send <= INT_MAX);

    long *recv = teal2_calloc(num_adjncys, sizeof(*recv));
    sync2_collect(part, recv, adjncy, (int)num_send, num_adjncys, MPI_LONG, 1);

    return recv;
}

static void decompose_comm(const Dual *dual, const long *part_, const Grid *grid)
{
    int num_adjncys;
    long *adjncy = collect_adjncys(dual, part_, grid, &num_adjncys);
    long *part = collect_parts(adjncy, dual, part_, num_adjncys);

    int *count = teal2_calloc(sync2.size, sizeof(*count));
    for (int i = 0; i < num_adjncys; i++) {
        if (part[i] != sync2.rank) {
            count[part[i]] += 1;
        }
    }

    int degree = 0;
    for (int i = 0; i < sync2.size; i++) {
        if (count[i] > 0) {
            degree += 1;
        }
    }

    int *rank = teal2_calloc(degree, sizeof(*rank));
    int *weight = teal2_calloc(degree, sizeof(*weight));

    int num = 0;
    for (int i = 0; i < sync2.size; i++) {
        if (count[i] > 0) {
            rank[num] = i;
            weight[num] = count[i];
            num += 1;
        }
    }
    assert(num == degree);

    int reorder = 1;  // WARN: I have not found a cluster yet which actually does this
    MPI_Comm comm;
    int ret = MPI_Dist_graph_create_adjacent(sync2.comm, degree, rank, weight, degree, rank, weight,
                                             MPI_INFO_NULL, reorder, &comm);
    if (ret != MPI_SUCCESS) {
        teal2_error("MPI_Dist_graph_create_adjacent failed (%d)", ret);
    }
    sync2_reinit(comm);

    teal2_free(adjncy);
    teal2_free(part);
    teal2_free(count);
    teal2_free(rank);
    teal2_free(weight);
}

typedef struct {
    int entity;
    int rank;
    int idx;
    int num_nodes;
    long node[MAX_CELL_NODES];
} Cell;

static MPI_Datatype datatype_cell(void)
{
    MPI_Datatype datatype;
    int len[5] = {1, 1, 1, 1, MAX_CELL_NODES};
    MPI_Aint off[5] = {offsetof(Cell, entity), offsetof(Cell, rank), offsetof(Cell, idx),
                       offsetof(Cell, num_nodes), offsetof(Cell, node)};
    MPI_Datatype type[5] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_LONG};
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
    return (lhs->rank > rhs->rank) - (lhs->rank < rhs->rank);
}

static void redistribute_cells(Grid *grid, const long *part)
{
    int *num_send = teal2_calloc(sync2.size, sizeof(*num_send));
    for (int i = 0; i < grid->cells.off_periodic; i++) {
        num_send[part[i]] += 1;
    }

    int *off_send = teal2_calloc(sync2.size + 1, sizeof(*off_send));
    for (int i = 0; i < sync2.size; i++) {
        off_send[i + 1] = off_send[i] + num_send[i];
    }

    int *cur_send = teal2_calloc(sync2.size, sizeof(*cur_send));
    copy(cur_send, off_send, sync2.size, sizeof(*cur_send));

    Cell *send = teal2_calloc(grid->cells.off_periodic, sizeof(*send));
    for (int i = 0; i < grid->entities.num; i++) {
        for (int j = grid->entities.cell_off[i]; j < grid->entities.cell_off[i + 1]; j++) {
            Cell *cell = &send[cur_send[part[j]]++];
            cell->entity = i;
            cell->num_nodes = 0;
            for (long k = grid->cells.node_off[j]; k < grid->cells.node_off[j + 1]; k++) {
                assert(cell->num_nodes < MAX_CELL_NODES);
                cell->node[cell->num_nodes++] = grid->cells.node_idx[k];
            }
        }
    }
    for (int i = 0; i < sync2.size; i++) {
        assert(cur_send[i] == off_send[i + 1]);
    }

    int *num_recv = teal2_calloc(sync2.size, sizeof(*num_recv));
    MPI_Alltoall(num_send, 1, MPI_INT, num_recv, 1, MPI_INT, sync2.comm);

    int *off_recv = teal2_calloc(sync2.size + 1, sizeof(*off_recv));
    for (int i = 0; i < sync2.size; i++) {
        off_recv[i + 1] = off_recv[i] + num_recv[i];
    }

    int tot_recv = off_recv[sync2.size];
    Cell *recv = teal2_calloc(tot_recv, sizeof(*recv));

    MPI_Datatype type = datatype_cell();
    MPI_Alltoallv(send, num_send, off_send, type, recv, num_recv, off_recv, type, sync2.comm);
    MPI_Type_free(&type);

    sort(recv, tot_recv, sizeof(*recv), compare_cell);

    int *cell_off = teal2_calloc(grid->entities.num + 1, sizeof(*cell_off));
    long *node_off = teal2_calloc(tot_recv + 1, sizeof(*node_off));
    long *node_idx = teal2_calloc(tot_recv * MAX_CELL_NODES, sizeof(*node_idx));

    for (int i = 0; i < tot_recv; i++) {
        cell_off[recv[i].entity + 1] += 1;
        node_off[i + 1] = node_off[i];
        for (int j = 0; j < recv[i].num_nodes; j++) {
            assert(node_off[i + 1] - node_off[i] < MAX_CELL_NODES);
            node_idx[node_off[i + 1]++] = recv[i].node[j];
        }
    }
    for (int i = 0; i < grid->entities.num; i++) {
        cell_off[i + 1] += cell_off[i];
    }

    teal2_free(grid->cells.node_off);
    teal2_free(grid->cells.node_idx);
    teal2_free(grid->entities.cell_off);

    grid->cells.num = tot_recv;
    grid->cells.num_inner = cell_off[grid->entities.num_inner];
    grid->cells.off_boundary = cell_off[grid->entities.off_boundary];
    grid->cells.off_periodic = cell_off[grid->entities.num];
    assert(grid->cells.off_periodic == grid->cells.num);

    assert(node_off[tot_recv] <= INT_MAX);
    int num_indices = (int)node_off[tot_recv];

    grid->cells.node_off = node_off;
    grid->cells.node_idx = teal2_realloc(node_idx, num_indices, sizeof(*node_idx));

    grid->entities.cell_off = cell_off;

    teal2_free(num_send);
    teal2_free(off_send);
    teal2_free(cur_send);
    teal2_free(send);
    teal2_free(num_recv);
    teal2_free(off_recv);
    teal2_free(recv);
}

static void redistribute_nodes(Grid *grid, int dst)
{
    int *rank = teal2_calloc(sync2.size, sizeof(*rank));
    sync2_gather(&dst, rank, 1, MPI_INT);

    int src = -1;
    for (int i = 0; i < sync2.size; i++) {
        if (rank[i] == sync2.rank) {
            src = i;
            break;
        }
    }
    assert(src != -1);

    teal2_free(rank);

    if (dst == src) {
        return;
    }

    int num_nodes;
    sync2_exchange(&grid->nodes.num, &num_nodes, 1, 1, dst, src, MPI_INT, 1);

    Vector *coord = teal2_calloc(num_nodes, sizeof(*coord));
    sync2_exchange(grid->nodes.coord, coord, grid->nodes.num, num_nodes, dst, src, MPI_DOUBLE, 3);

    teal2_free(grid->nodes.coord);

    grid->nodes.num = num_nodes;
    grid->nodes.coord = coord;
}

static void partition_cells(Grid *grid)
{
    Dual *dual = dual_init(grid);
    long *part = dual_partition(dual, grid);

    if (grid->entities.num > grid->entities.off_boundary) {
        dual_periodic(dual, grid);
    }

    int rank = sync2.rank;
    decompose_comm(dual, part, grid);

    redistribute_cells(grid, part);
    redistribute_nodes(grid, rank);

    teal2_free(part);
    dual_deinit(dual);
}

static void reorder_periodic_cells(Grid *grid, const Cell *recv, int num_periodic)
{
    int num = 0;
    for (int i = grid->entities.off_boundary; i < grid->entities.num; i++) {
        for (int j = grid->entities.cell_off[i]; j < grid->entities.cell_off[i + 1]; j++) {
            grid->cells.node_off[j + 1] = grid->cells.node_off[j];
            for (int k = 0; k < recv[num].num_nodes; k++) {
                assert(grid->cells.node_off[j + 1] - grid->cells.node_off[j] < MAX_CELL_NODES);
                grid->cells.node_idx[grid->cells.node_off[j + 1]++] = recv[num].node[k];
            }
            num += 1;
        }
    }
    assert(num == num_periodic);
}

static void create_periodic_graph(Mesh2 *mesh, const Grid *grid, const Cell *recv, int num_periodic)
{
    int num_neighbors = 1;
    for (int i = 1; i < num_periodic; i++) {
        if (compare_cell(&recv[i], &recv[i - 1]) != 0) {
            num_neighbors += 1;
        }
    }

    int (*tag)[2] = teal2_calloc(num_neighbors, sizeof(*tag));
    int *rank = teal2_calloc(num_neighbors, sizeof(*rank));
    int *recv_off = teal2_calloc(num_neighbors + 1, sizeof(*recv_off));
    int *send_off = teal2_calloc(num_neighbors + 1, sizeof(*send_off));
    int *send_idx = teal2_calloc(num_periodic, sizeof(*send_idx));

    tag[0][0] = recv[0].entity;
    tag[0][1] = grid->entities.periodic[recv[0].entity];
    rank[0] = recv[0].rank;
    recv_off[0] = grid->cells.off_boundary;

    int num = 1;
    for (int i = 1; i < num_periodic; i++) {
        if (compare_cell(&recv[i], &recv[i - 1]) != 0) {
            tag[num][0] = recv[i].entity;
            tag[num][1] = grid->entities.periodic[recv[i].entity];
            rank[num] = recv[i].rank;
            recv_off[num] = recv_off[0] + i;
            send_off[num] = i;
            num += 1;
        }
    }
    assert(num == num_neighbors);

    recv_off[num] = grid->cells.off_periodic;
    send_off[num] = num_periodic;

    int *recv_idx = teal2_calloc(num_periodic, sizeof(*recv_idx));
    for (int i = 0; i < num_periodic; i++) {
        recv_idx[i] = recv[i].idx;
    }

    MPI_Request *req_recv = teal2_calloc(num_neighbors, sizeof(*req_recv));
    MPI_Request *req_send = teal2_calloc(num_neighbors, sizeof(*req_send));
    for (int i = 0; i < num_neighbors; i++) {
        int off = send_off[i];
        int count = recv_off[i + 1] - recv_off[i];
        MPI_Irecv(&send_idx[off], count, MPI_INT, rank[i], tag[i][0], sync2.comm, &req_recv[i]);
        MPI_Isend(&recv_idx[off], count, MPI_INT, rank[i], tag[i][1], sync2.comm, &req_send[i]);
    }
    MPI_Waitall(num_neighbors, req_recv, MPI_STATUSES_IGNORE);
    MPI_Waitall(num_neighbors, req_send, MPI_STATUSES_IGNORE);

    mesh->neighbors.num = num_neighbors;
    mesh->neighbors.tag = tag;
    mesh->neighbors.rank = rank;
    mesh->neighbors.recv_off = recv_off;
    mesh->neighbors.send.off = send_off;
    mesh->neighbors.send.idx = send_idx;

    teal2_free(recv_idx);
    teal2_free(req_recv);
    teal2_free(req_send);
}

static void collect_periodic(Mesh2 *mesh, Grid *grid)
{
    Dual *dual = dual_init(grid);
    dual_periodic(dual, grid);

    if (grid->cells.off_periodic == grid->cells.off_boundary) {
        goto out;
    }

    int num_periodic = grid->cells.off_periodic - grid->cells.off_boundary;
    Cell *recv = teal2_calloc(num_periodic, sizeof(*recv));

    int num = 0;
    for (int i = grid->entities.off_boundary; i < grid->entities.num; i++) {
        for (int j = grid->entities.cell_off[i]; j < grid->entities.cell_off[i + 1]; j++) {
            assert(dual->xadj[j + 1] - dual->xadj[j] == 2);
            long adjncy = dual->adjncy[dual->xadj[j] + 1];
            int rank = digitize(&adjncy, dual->dist, sync2.size, sizeof(*dual->dist), compare_long);
            assert(0 <= rank && rank < sync2.size);
            long idx = adjncy - dual->dist[rank];
            assert(idx <= INT_MAX);
            recv[num].entity = i;
            recv[num].rank = rank;
            recv[num].idx = (int)idx;
            recv[num].num_nodes = 0;
            for (long k = grid->cells.node_off[j]; k < grid->cells.node_off[j + 1]; k++) {
                assert(recv[num].num_nodes < MAX_CELL_NODES);
                recv[num].node[recv[num].num_nodes++] = grid->cells.node_idx[k];
            }
            num += 1;
        }
    }
    assert(num == num_periodic);

    sort(recv, num_periodic, sizeof(*recv), compare_cell);

    reorder_periodic_cells(grid, recv, num_periodic);
    create_periodic_graph(mesh, grid, recv, num_periodic);

    teal2_free(recv);
out:
    dual_deinit(dual);
}

static void append_neighbor_cells(Grid *grid, const Cell *recv, int tot_recv)
{
    assert(grid->cells.node_off[grid->cells.num_inner] <= INT_MAX);
    int num_inner = (int)grid->cells.node_off[grid->cells.num_inner];
    long *inner = teal2_calloc(num_inner, sizeof(*inner));
    copy(inner, grid->cells.node_idx, num_inner, sizeof(*inner));
    unique(inner, num_inner, sizeof(*inner), compare_long);

    int num_cells = grid->cells.off_periodic + tot_recv;
    int max_indices = num_cells * MAX_CELL_NODES;

    long *node_off = teal2_realloc(grid->cells.node_off, num_cells + 1, sizeof(*node_off));
    long *node_idx = teal2_realloc(grid->cells.node_idx, max_indices, sizeof(*node_idx));

    int num = grid->cells.off_periodic;
    for (int i = 0; i < tot_recv; i++) {
        node_off[num + 1] = node_off[num];
        for (int j = 0; j < recv[i].num_nodes; j++) {
            if (search(&recv[i].node[j], inner, num_inner, sizeof(*inner), compare_long)) {
                assert(node_off[num + 1] - node_off[num] < MAX_CELL_NODES);
                node_idx[node_off[num + 1]++] = recv[i].node[j];
            }
        }
        num += 1;
    }
    assert(num == num_cells);

    assert(node_off[num_cells] <= INT_MAX);
    int num_indices = (int)node_off[num_cells];

    grid->cells.num = num_cells;
    grid->cells.node_off = node_off;
    grid->cells.node_idx = teal2_realloc(node_idx, num_indices, sizeof(*node_idx));

    teal2_free(inner);
}

static void append_neighbor_graph(Mesh2 *mesh, const int *num_recv, const int *num_send,
                                  const int *off_send, const int *idx_send, const Grid *grid)
{
    int num_neighbors = mesh->neighbors.num;
    for (int i = 0; i < sync2.size; i++) {
        if (num_recv[i] || num_send[i]) {
            num_neighbors += 1;
        }
    }

    int num_periodic = grid->cells.off_periodic - grid->cells.off_boundary;
    int num_indices = off_send[sync2.size] + num_periodic;

    int (*tag)[2] = teal2_realloc(mesh->neighbors.tag, num_neighbors, sizeof(*tag));
    int *rank = teal2_realloc(mesh->neighbors.rank, num_neighbors, sizeof(*rank));
    int *recv_off = teal2_realloc(mesh->neighbors.recv_off, num_neighbors + 1, sizeof(*recv_off));
    int *send_off = teal2_realloc(mesh->neighbors.send.off, num_neighbors + 1, sizeof(*send_off));
    int *send_idx = teal2_realloc(mesh->neighbors.send.idx, num_indices, sizeof(*send_idx));

    recv_off[mesh->neighbors.num] = grid->cells.off_periodic;
    send_off[mesh->neighbors.num] = num_periodic;

    int num = mesh->neighbors.num;
    for (int i = 0; i < sync2.size; i++) {
        if (num_recv[i] || num_send[i]) {
            tag[num][0] = grid->entities.num;
            tag[num][1] = grid->entities.num;
            rank[num] = i;
            recv_off[num + 1] = recv_off[num] + num_recv[i];
            send_off[num + 1] = send_off[num];
            for (int j = off_send[i]; j < off_send[i + 1]; j++) {
                send_idx[send_off[num + 1]++] = idx_send[j];
            }
            num += 1;
        }
    }
    assert(num == num_neighbors);

    mesh->neighbors.num = num_neighbors;
    mesh->neighbors.tag = tag;
    mesh->neighbors.rank = rank;
    mesh->neighbors.recv_off = recv_off;
    mesh->neighbors.send.off = send_off;
    mesh->neighbors.send.idx = send_idx;
}

static void create_neighbors(Mesh2 *mesh, Grid *grid)
{
    if (grid->entities.num > grid->entities.off_boundary) {
        collect_periodic(mesh, grid);
    }

    Dual *dual = dual_init(grid);

    assert(dual->xadj[grid->cells.num_inner] <= INT_MAX);
    int num_adjncys = (int)dual->xadj[grid->cells.num_inner];

    long *adjncy = teal2_calloc(num_adjncys, sizeof(*adjncy));
    copy(adjncy, dual->adjncy, num_adjncys, sizeof(*adjncy));
    num_adjncys = unique(adjncy, num_adjncys, sizeof(*adjncy), compare_long);

    int *num_inner = teal2_calloc(sync2.size, sizeof(*num_inner));
    sync2_gather(&grid->cells.num_inner, num_inner, 1, MPI_INT);

    int num = 0;
    while (num < num_adjncys) {
        int rank =
            digitize(&adjncy[num], dual->dist, sync2.size, sizeof(*dual->dist), compare_long);
        long beg_inner = dual->dist[rank];
        long end_inner = dual->dist[rank] + num_inner[rank];
        if (rank == sync2.rank || adjncy[num] < beg_inner || end_inner <= adjncy[num]) {
            adjncy[num] = adjncy[--num_adjncys];
        }
        else {
            num += 1;
        }
    }
    sort(adjncy, num_adjncys, sizeof(*adjncy), compare_long);

    int *rank = teal2_calloc(num_adjncys, sizeof(*rank));
    for (int i = 0; i < num_adjncys; i++) {
        rank[i] = digitize(&adjncy[i], dual->dist, sync2.size, sizeof(*dual->dist), compare_long);
        assert(0 <= rank[i] && rank[i] < sync2.size);
    }

    int *num_recv = teal2_calloc(sync2.size, sizeof(*num_recv));
    for (int i = 0; i < num_adjncys; i++) {
        num_recv[rank[i]] += 1;
    }

    int *off_recv = teal2_calloc(sync2.size + 1, sizeof(*off_recv));
    for (int i = 0; i < sync2.size; i++) {
        off_recv[i + 1] = off_recv[i] + num_recv[i];
    }

    int *cur_recv = teal2_calloc(sync2.size, sizeof(*cur_recv));
    copy(cur_recv, off_recv, sync2.size, sizeof(*cur_recv));

    int *idx_recv = teal2_calloc(num_adjncys, sizeof(*idx_recv));
    for (int i = 0; i < num_adjncys; i++) {
        long idx_local = adjncy[i] - dual->dist[rank[i]];
        assert(idx_local <= INT_MAX);
        idx_recv[cur_recv[rank[i]]++] = (int)idx_local;
    }
    for (int i = 0; i < sync2.size; i++) {
        assert(cur_recv[i] == off_recv[i + 1]);
    }

    int *num_send = teal2_calloc(sync2.size, sizeof(*num_send));
    MPI_Alltoall(num_recv, 1, MPI_INT, num_send, 1, MPI_INT, sync2.comm);

    int *off_send = teal2_calloc(sync2.size + 1, sizeof(*off_send));
    for (int i = 0; i < sync2.size; i++) {
        off_send[i + 1] = off_send[i] + num_send[i];
    }

    int tot_send = off_send[sync2.size];
    int *idx_send = teal2_calloc(tot_send, sizeof(*idx_send));
    MPI_Alltoallv(idx_recv, num_recv, off_recv, MPI_INT, idx_send, num_send, off_send, MPI_INT,
                  sync2.comm);

    Cell *send = teal2_calloc(tot_send, sizeof(*send));
    for (int i = 0; i < tot_send; i++) {
        int idx = idx_send[i];
        send[i].num_nodes = 0;
        for (long j = grid->cells.node_off[idx]; j < grid->cells.node_off[idx + 1]; j++) {
            assert(send[i].num_nodes < MAX_CELL_NODES);
            send[i].node[send[i].num_nodes++] = grid->cells.node_idx[j];
        }
    }

    int tot_recv = off_recv[sync2.size];
    Cell *recv = teal2_calloc(tot_recv, sizeof(*recv));

    MPI_Datatype type = datatype_cell();
    MPI_Alltoallv(send, num_send, off_send, type, recv, num_recv, off_recv, type, sync2.comm);
    MPI_Type_free(&type);

    append_neighbor_cells(grid, recv, tot_recv);
    append_neighbor_graph(mesh, num_recv, num_send, off_send, idx_send, grid);

    teal2_free(adjncy);
    teal2_free(num_inner);
    teal2_free(rank);
    teal2_free(num_recv);
    teal2_free(off_recv);
    teal2_free(cur_recv);
    teal2_free(idx_recv);
    teal2_free(num_send);
    teal2_free(off_send);
    teal2_free(idx_send);
    teal2_free(send);
    teal2_free(recv);
    dual_deinit(dual);
}

typedef struct {
    int rank;
    long old;
    long new;
} Node;

static MPI_Datatype datatype_node(void)
{
    MPI_Datatype datatype;
    int len[3] = {1, 1, 1};
    MPI_Aint off[3] = {offsetof(Node, rank), offsetof(Node, old), offsetof(Node, new)};
    MPI_Datatype type[3] = {MPI_INT, MPI_LONG, MPI_LONG};
    MPI_Type_create_struct(3, len, off, type, &datatype);
    return sync2_resized(datatype, sizeof(Node));
}

static int compare_node(const void *lhs_, const void *rhs_)
{
    const Node *lhs = lhs_;
    const Node *rhs = rhs_;
    return (lhs->rank > rhs->rank) - (lhs->rank < rhs->rank);
}

static void determine_owners(Node *node, const long *node_idx, int num_nodes, int cap)
{
    MPI_Datatype type = datatype_node();
    int num = num_nodes;
    for (int rank = 0; rank < sync2.size; rank++) {
        for (int i = 0; i < num; i++) {
            long key = node[i].old;
            if (search(&key, node_idx, num_nodes, sizeof(*node_idx), compare_long)) {
                if (sync2.rank < node[i].rank) {
                    node[i].rank = sync2.rank;
                }
            }
        }
        sync2_rotate(node, &num, cap, type, 1);
    }
    assert(num == num_nodes);
    MPI_Type_free(&type);
}

typedef struct {
    long old, new, idx;
} Map;

static int compare_map(const void *lhs_, const void *rhs_)
{
    const Map *lhs = lhs_;
    const Map *rhs = rhs_;
    return (lhs->old > rhs->old) - (lhs->old < rhs->old);
}

static void propagate_indices(Node *node, const Map *map, int num_nodes, int num_inner, int cap)
{
    MPI_Datatype type = datatype_node();
    int num = num_nodes;
    for (int rank = 0; rank < sync2.size; rank++) {
        for (int i = 0; i < num; i++) {
            if (node[i].new != -1) {
                continue;
            }
            Map key = {.old = node[i].old};
            Map *val = search(&key, map, num_inner, sizeof(*map), compare_map);
            if (val) {
                node[i].new = val->new;
            }
        }
        sync2_rotate(node, &num, cap, type, 1);
    }
    assert(num == num_nodes);
    MPI_Type_free(&type);
}

static long *collect_tags(const Node *node, int num_nodes)
{
    long *tag = teal2_calloc(num_nodes, sizeof(*tag));
    for (int i = 0; i < num_nodes; i++) {
        if (node[i].new == -1) {
            teal2_error("global node index was not converted (%ld)", node[i].old);
        }
        tag[i] = node[i].new;
    }
    return tag;
}

static Vector *collect_coords(const Grid *grid, const Map *map, const long *node_idx, int num_nodes)
{
    Vector *coord = teal2_calloc(num_nodes, sizeof(*coord));
    sync2_collect(grid->nodes.coord, coord, node_idx, grid->nodes.num, num_nodes, MPI_DOUBLE, 3);

    Vector *sorted = teal2_calloc(num_nodes, sizeof(*grid->nodes.coord));
    for (int i = 0; i < num_nodes; i++) {
        Map key = {.old = node_idx[i]};
        Map *val = search(&key, map, num_nodes, sizeof(*map), compare_map);
        if (!val) {
            teal2_error("node index not found (%ld)", node_idx[i]);
        }
        sorted[val->idx] = coord[i];
    }

    teal2_free(coord);
    return sorted;
}

static void remap_cells_node_graph(Grid *grid, const Map *map, int num_nodes)
{
    for (int i = 0; i < grid->cells.num; i++) {
        for (long j = grid->cells.node_off[i]; j < grid->cells.node_off[i + 1]; j++) {
            Map key = {.old = grid->cells.node_idx[j]};
            Map *val = search(&key, map, num_nodes, sizeof(*map), compare_map);
            if (!val) {
                teal2_error("node index not found (%ld)", grid->cells.node_idx[j]);
            }
            grid->cells.node_idx[j] = val->idx;
        }
    }
}

static void remap_entities_node_graph(Grid *grid, const Map *map, int num_nodes)
{
    if (grid->entities.num == grid->entities.off_boundary) {
        return;
    }

    assert(grid->entities.node_off[grid->entities.num] <= INT_MAX);
    int num = (int)grid->entities.node_off[grid->entities.num];
    int cap = num;
    sync2_max(&cap, 1, MPI_INT);

    long *node_idx = teal2_calloc(cap, sizeof(*node_idx));
    for (int i = 0; i < num; i++) {
        node_idx[i] = -(grid->entities.node_idx[i] + 1);
    }

    for (int rank = 0; rank < sync2.size; rank++) {
        for (int i = 0; i < num; i++) {
            if (node_idx[i] < 0) {
                Map key = {.old = -(node_idx[i] + 1)};
                Map *val = search(&key, map, num_nodes, sizeof(*map), compare_map);
                if (val) {
                    node_idx[i] = val->new;
                }
            }
        }
        sync2_rotate(node_idx, &num, cap, MPI_LONG, 1);
    }
    assert(num == grid->entities.node_off[grid->entities.num]);

    for (int i = 0; i < num; i++) {
        if (node_idx[i] < 0) {
            teal2_error("periodic node index was not converted (%ld)", node_idx[i]);
        }
    }

    copy(grid->entities.node_idx, node_idx, num, sizeof(*node_idx));

    teal2_free(node_idx);
}

static int partition_nodes(Grid *grid)
{
    assert(grid->cells.node_off[grid->cells.num] <= INT_MAX);
    int max_nodes = (int)grid->cells.node_off[grid->cells.num];

    long *node_idx = teal2_calloc(max_nodes, sizeof(*node_idx));
    int num = 0;
    for (int i = 0; i < grid->cells.num; i++) {
        for (long j = grid->cells.node_off[i]; j < grid->cells.node_off[i + 1]; j++) {
            node_idx[num++] = grid->cells.node_idx[j];
        }
    }
    int num_nodes = unique(node_idx, num, sizeof(*node_idx), compare_long);

    int cap = num_nodes;
    sync2_max(&cap, 1, MPI_INT);

    Node *node = teal2_calloc(cap, sizeof(*node));
    for (int i = 0; i < num_nodes; i++) {
        node[i].rank = sync2.rank;
        node[i].old = node_idx[i];
        node[i].new = -1;
    }
    determine_owners(node, node_idx, num_nodes, cap);

    int num_inner = 0;
    for (int i = 0; i < num_nodes; i++) {
        if (node[i].rank == sync2.rank) {
            node[i].rank = -node[i].rank;
            num_inner += 1;
        }
    }
    sort(node, num_nodes, sizeof(*node), compare_node);

    long prefix = num_inner;
    sync2_prefix(&prefix, 1, MPI_LONG);

    Map *map = teal2_calloc(num_nodes, sizeof(*map));
    num = 0;
    for (int i = 0; i < num_nodes; i++) {
        if (node[i].rank == -sync2.rank) {
            node[i].new = prefix + i;
            map[num].old = node[i].old;
            map[num].new = node[i].new;
            num += 1;
        }
    }
    assert(num == num_inner);

    sort(map, num_inner, sizeof(*map), compare_map);
    propagate_indices(node, map, num_nodes, num_inner, cap);

    for (int i = 0; i < num_nodes; i++) {
        map[i].old = node[i].old;
        map[i].new = node[i].new;
        map[i].idx = i;
    }
    sort(map, num_nodes, sizeof(*map), compare_map);

    long *tag = collect_tags(node, num_nodes);
    Vector *coord = collect_coords(grid, map, node_idx, num_nodes);

    teal2_free(grid->nodes.tag);
    teal2_free(grid->nodes.coord);

    grid->nodes.num = num_nodes;
    grid->nodes.tag = tag;
    grid->nodes.coord = coord;

    remap_cells_node_graph(grid, map, num_nodes);
    remap_entities_node_graph(grid, map, num_nodes);

    teal2_free(node_idx);
    teal2_free(node);
    teal2_free(map);
    return num_inner;
}

static void create_nodes(Mesh2 *mesh, Grid *grid, int num_inner)
{
    long *global = teal2_calloc(grid->nodes.num, sizeof(*global));
    for (int i = 0; i < grid->nodes.num; i++) {
        global[i] = grid->nodes.tag[i];
    }

    mesh->nodes.num = grid->nodes.num;
    mesh->nodes.num_inner = num_inner;
    mesh->nodes.global = global;
    mesh->nodes.coord = grid->nodes.coord;

    grid->nodes.coord = 0;
}

static void create_cells(Mesh2 *mesh, const Grid *grid)
{
    int *node_off = teal2_calloc(grid->cells.num + 1, sizeof(*node_off));
    for (int i = 0; i < grid->cells.num + 1; i++) {
        assert(grid->cells.node_off[i] <= INT_MAX);
        node_off[i] = (int)grid->cells.node_off[i];
    }

    int *node_idx = teal2_calloc(node_off[grid->cells.num], sizeof(*node_idx));
    for (int i = 0; i < node_off[grid->cells.num]; i++) {
        assert(grid->cells.node_idx[i] <= INT_MAX);
        node_idx[i] = (int)grid->cells.node_idx[i];
    }

    mesh->cells.num = grid->cells.num;
    mesh->cells.num_inner = grid->cells.num_inner;
    mesh->cells.off_boundary = grid->cells.off_boundary;
    mesh->cells.off_periodic = grid->cells.off_periodic;
    mesh->cells.node.off = node_off;
    mesh->cells.node.idx = node_idx;
}

static void create_entities(Mesh2 *mesh, Grid *grid)
{
    mesh->entities.num = grid->entities.num;
    mesh->entities.num_inner = grid->entities.num_inner;
    mesh->entities.off_boundary = grid->entities.off_boundary;
    mesh->entities.name = grid->entities.name;
    mesh->entities.cell_off = grid->entities.cell_off;
    mesh->entities.rotation = grid->entities.rotation;
    mesh->entities.translation = grid->entities.translation;

    grid->entities.name = 0;
    grid->entities.cell_off = 0;
    grid->entities.rotation = 0;
    grid->entities.translation = 0;
}

Mesh2 *mesh2_read(const char *fname)
{
    assert(fname);

    Grid *grid = grid_init(fname);
    partition_cells(grid);

    Mesh2 *mesh = teal2_calloc(1, sizeof(*mesh));
    create_neighbors(mesh, grid);

    int num_inner = partition_nodes(grid);
    create_nodes(mesh, grid, num_inner);
    create_cells(mesh, grid);
    create_entities(mesh, grid);

    grid_deinit(grid);
    return mesh;
}
