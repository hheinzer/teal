#include "read.h"

#include <parmetis.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mesh.h"
#include "reorder.h"
#include "teal/arena.h"
#include "teal/array.h"
#include "teal/assert.h"
#include "teal/kdtree.h"
#include "teal/option.h"
#include "teal/sync.h"
#include "teal/utils.h"
#include "teal/vector.h"

/* Dispatch to reader based on file extension. */
static void read_file(Mesh *mesh, const char *fname)
{
    char *ext = strrchr(fname, '.');
    if (ext && !strcmp(ext, ".h5")) {
        mesh_read_hdf5(mesh, fname);
    }
    else if (ext && !strcmp(ext, ".msh")) {
        mesh_read_gmsh(mesh, fname);
    }
    else {
        assert(false);
    }
}

static long count_inner_cells(const MeshEntities *entities)
{
    long count = 0;
    for (long i = 0; i < entities->num_inner; i++) {
        count += entities->cell_off[i + 1] - entities->cell_off[i];
    }
    return count;
}

/* Gather coordinates for requested global node ids to their owning ranks. */
static void collect_coords(const MeshNodes *nodes, const long *global, vector *coord, long num)
{
    Arena save = arena_save();

    long *num_nodes = arena_malloc(sync.size + 1, sizeof(*num_nodes));
    num_nodes[0] = 0;
    MPI_Allgather(&nodes->num, 1, MPI_LONG, &num_nodes[1], 1, MPI_LONG, sync.comm);
    for (long i = 0; i < sync.size; i++) {
        num_nodes[i + 1] += num_nodes[i];
    }

    int *num_recv = arena_calloc(sync.size, sizeof(*num_recv));
    for (long i = 0; i < num; i++) {
        long rank = array_ldigitize(&num_nodes[1], global[i], sync.size);
        num_recv[rank] += 1;
    }

    int *off_recv = arena_malloc(sync.size + 1, sizeof(*off_recv));
    off_recv[0] = 0;
    for (long i = 0; i < sync.size; i++) {
        off_recv[i + 1] = off_recv[i] + num_recv[i];
    }

    long *idx_recv = arena_malloc(num, sizeof(*idx_recv));
    for (long i = 0; i < sync.size; i++) {
        off_recv[i + 1] -= num_recv[i];
    }
    for (long i = 0; i < num; i++) {
        long rank = array_ldigitize(&num_nodes[1], global[i], sync.size);
        idx_recv[off_recv[rank + 1]++] = global[i] - num_nodes[rank];
    }
    for (long i = 0; i < sync.size; i++) {
        assert(off_recv[i + 1] - off_recv[i] == num_recv[i]);
    }

    int *num_send = arena_malloc(sync.size, sizeof(*num_send));
    MPI_Alltoall(num_recv, 1, MPI_INT, num_send, 1, MPI_INT, sync.comm);

    int *off_send = arena_malloc(sync.size + 1, sizeof(*off_send));
    off_send[0] = 0;
    for (long i = 0; i < sync.size; i++) {
        off_send[i + 1] = off_send[i] + num_send[i];
    }

    long tot_send = off_send[sync.size];
    long *idx_send = arena_malloc(tot_send, sizeof(*idx_send));
    MPI_Alltoallv(idx_recv, num_recv, off_recv, MPI_LONG, idx_send, num_send, off_send, MPI_LONG,
                  sync.comm);

    vector *coord_send = arena_malloc(tot_send, sizeof(*coord_send));
    for (long i = 0; i < tot_send; i++) {
        coord_send[i] = nodes->coord[idx_send[i]];
    }

    MPI_Datatype type;
    MPI_Type_contiguous(sizeof(vector), MPI_BYTE, &type);
    MPI_Type_commit(&type);

    vector *coord_recv = arena_malloc(num, sizeof(*coord_recv));
    MPI_Alltoallv(coord_send, num_send, off_send, type, coord_recv, num_recv, off_recv, type,
                  sync.comm);
    MPI_Type_free(&type);

    for (long i = 0; i < sync.size; i++) {
        off_recv[i + 1] -= num_recv[i];
    }
    for (long i = 0; i < num; i++) {
        long rank = array_ldigitize(&num_nodes[1], global[i], sync.size);
        coord[i] = coord_recv[off_recv[rank + 1]++];
    }
    for (long i = 0; i < sync.size; i++) {
        assert(off_recv[i + 1] - off_recv[i] == num_recv[i]);
    }

    arena_load(save);
}

typedef enum { LHS, RHS } Side;

typedef struct {
    Side side;
    long cell;
    long peer;
} Link;

/* Discover periodic partners by ring-rotating shifted centers. */
static void collect_links(const vector *mean, Kdtree *center2link)
{
    Arena save = arena_save();

    typedef struct {
        vector coord;
        long global;
    } Center;

    long cap = sync_lmax(center2link->num);
    Center *center = arena_calloc(cap, sizeof(*center));

    long num = 0;
    for (KdtreeItem *item = center2link->beg; item; item = item->next) {
        vector coord = item->key;
        Link *link = item->val;
        if (link->side == LHS) {
            center[num].coord = vector_add(coord, vector_sub(mean[RHS], mean[LHS]));
        }
        else {
            center[num].coord = vector_add(coord, vector_sub(mean[LHS], mean[RHS]));
        }
        center[num].global = link->peer;
        num += 1;
    }
    assert(num == center2link->num);

    MPI_Datatype type;
    MPI_Type_contiguous(sizeof(Center), MPI_BYTE, &type);
    MPI_Type_commit(&type);

    int dst = (sync.rank + 1) % sync.size;
    int src = (sync.rank - 1 + sync.size) % sync.size;
    int tag = sync_tag();
    for (long step = 0; step < sync.size; step++) {
        for (long i = 0; i < cap; i++) {
            if (center[i].global == -1) {
                Link *link = kdtree_lookup(center2link, center[i].coord);
                if (link) {
                    center[i].global = link->cell;
                }
            }
        }
        MPI_Sendrecv_replace(center, cap, type, dst, tag, src, tag, sync.comm, MPI_STATUS_IGNORE);
    }
    MPI_Type_free(&type);

    num = 0;
    for (KdtreeItem *item = center2link->beg; item; item = item->next) {
        Link *link = item->val;
        link->peer = center[num].global;
        assert(link->peer != -1);
        num += 1;
    }
    assert(num == center2link->num);

    arena_load(save);
}

typedef struct {
    long cell;
    long peer;
} Edge;

static int cmp_edge(const void *lhs_, const void *rhs_)
{
    const Edge *lhs = lhs_;
    const Edge *rhs = rhs_;
    return cmp_asc(lhs->cell, rhs->cell);
}

/* For each periodic link, request the peer cell's global id from its owner. */
static Edge *collect_edges(const MeshCells *cells, const Kdtree *center2link, long *num_edges)
{
    Arena save = arena_save();

    long tot_send = center2link->num;
    Edge *send = arena_malloc(tot_send, sizeof(*send));
    long num = 0;
    for (KdtreeItem *item = center2link->beg; item; item = item->next) {
        Link *link = item->val;
        send[num].cell = link->cell;
        send[num].peer = link->peer;
        num += 1;
    }
    assert(num == tot_send);
    qsort(send, tot_send, sizeof(*send), cmp_edge);

    long *num_cells = arena_malloc(sync.size + 1, sizeof(*num_cells));
    num_cells[0] = 0;
    MPI_Allgather(&cells->num, 1, MPI_LONG, &num_cells[1], 1, MPI_LONG, sync.comm);
    for (long i = 0; i < sync.size; i++) {
        num_cells[i + 1] += num_cells[i];
    }

    int *num_send = arena_calloc(sync.size, sizeof(*num_send));
    for (long i = 0; i < tot_send; i++) {
        long rank = array_ldigitize(&num_cells[1], send[i].cell, sync.size);
        num_send[rank] += 1;
    }

    int *off_send = arena_malloc(sync.size + 1, sizeof(*off_send));
    off_send[0] = 0;
    for (long i = 0; i < sync.size; i++) {
        off_send[i + 1] = off_send[i] + num_send[i];
    }

    int *num_recv = arena_malloc(sync.size, sizeof(*num_recv));
    MPI_Alltoall(num_send, 1, MPI_INT, num_recv, 1, MPI_INT, sync.comm);

    int *off_recv = arena_malloc(sync.size + 1, sizeof(*off_recv));
    off_recv[0] = 0;
    for (long i = 0; i < sync.size; i++) {
        off_recv[i + 1] = off_recv[i] + num_recv[i];
    }

    MPI_Datatype type;
    MPI_Type_contiguous(sizeof(Edge), MPI_BYTE, &type);
    MPI_Type_commit(&type);

    long tot_recv = off_recv[sync.size];
    Edge *recv = arena_malloc(tot_recv, sizeof(*recv));
    MPI_Alltoallv(send, num_send, off_send, type, recv, num_recv, off_recv, type, sync.comm);
    MPI_Type_free(&type);

    long off_cells = sync_lexsum(cells->num);
    for (long i = 0; i < tot_recv; i++) {
        recv[i].cell -= off_cells;
    }
    qsort(recv, tot_recv, sizeof(*recv), cmp_edge);

    arena_load(save);

    *num_edges = tot_recv;
    return arena_smuggle(recv, *num_edges, sizeof(*recv));
}

typedef struct {
    idx_t *xadj, *adjncy;  // CSR graph (ParMETIS)
} Dual;

/* Augment the dual graph with periodic edges between entity sets `lhs` and `rhs`. */
static void connect_periodic(const MeshNodes *nodes, const MeshCells *cells,
                             const MeshEntities *entities, Dual *dual, long lhs, long rhs)
{
    Arena save = arena_save();

    long num_lhs = entities->cell_off[lhs + 1] - entities->cell_off[lhs];
    long num_rhs = entities->cell_off[rhs + 1] - entities->cell_off[rhs];

    long *global = arena_malloc((num_lhs + num_rhs) * MAX_CELL_NODES, sizeof(*global));
    long num_cells[2] = {0};
    long num = 0;
    for (long i = 0; i < entities->num; i++) {
        if (i == lhs || i == rhs) {
            for (long j = entities->cell_off[i]; j < entities->cell_off[i + 1]; j++) {
                for (long k = cells->node.off[j]; k < cells->node.off[j + 1]; k++) {
                    global[num++] = cells->node.idx[k];
                }
                Side side = (i == lhs) ? LHS : RHS;
                num_cells[side] += 1;
            }
        }
    }
    array_lunique(global, &num);

    vector *coord = arena_malloc(num, sizeof(*coord));
    collect_coords(nodes, global, coord, num);

    num_cells[LHS] = sync_lsum(num_cells[LHS]);
    num_cells[RHS] = sync_lsum(num_cells[RHS]);
    assert(num_cells[LHS] == num_cells[RHS]);

    Kdtree *center2link = kdtree_create(sizeof(Link));
    vector mean[2] = {0};
    long off_cells = sync_lexsum(cells->num);
    for (long i = 0; i < entities->num; i++) {
        if (i == lhs || i == rhs) {
            for (long j = entities->cell_off[i]; j < entities->cell_off[i + 1]; j++) {
                vector center = {0};
                long num_nodes = cells->node.off[j + 1] - cells->node.off[j];
                for (long k = cells->node.off[j]; k < cells->node.off[j + 1]; k++) {
                    long key = cells->node.idx[k];
                    long *val = bsearch(&key, global, num, sizeof(*global), cmp_long);
                    assert(val);
                    center = vector_add(center, vector_div(coord[val - global], num_nodes));
                }
                Side side = (i == lhs) ? LHS : RHS;
                long cell = j + off_cells;
                kdtree_insert(center2link, center, &(Link){side, cell, -1});
                mean[side] = vector_add(mean[side], center);
            }
        }
    }

    mean[LHS] = vector_div(sync_vsum(mean[LHS]), num_cells[LHS]);
    mean[RHS] = vector_div(sync_vsum(mean[RHS]), num_cells[RHS]);
    collect_links(mean, center2link);

    long num_edges;
    Edge *edge = collect_edges(cells, center2link, &num_edges);

    idx_t *xadj = malloc((cells->num + 1) * sizeof(*xadj));
    assert(xadj);

    idx_t *adjncy = malloc((dual->xadj[cells->num] + num_edges) * sizeof(*adjncy));
    assert(adjncy);

    xadj[0] = 0;

    for (long i = 0; i < cells->num; i++) {
        xadj[i + 1] = xadj[i];
        for (long j = dual->xadj[i]; j < dual->xadj[i + 1]; j++) {
            adjncy[xadj[i + 1]++] = dual->adjncy[j];
        }
        Edge key = {.cell = i};
        Edge *val = bsearch(&key, edge, num_edges, sizeof(*edge), cmp_edge);
        if (val) {
            adjncy[xadj[i + 1]++] = val->peer;
        }
    }
    free(dual->xadj);
    free(dual->adjncy);

    arena_load(save);

    dual->xadj = xadj;
    dual->adjncy = adjncy;
}

static int rotate_at_char(char *dst, const char *src, char sep)
{
    char *pos = strchr(src, sep);
    return pos ? sprintf(dst, "%s%c%.*s", pos + 1, sep, (int)(pos - src), src) : 0;
}

/* Build the dual graph (ParMETIS) and add edges for periodic entity pairs. */
static Dual connect_cells(const MeshNodes *nodes, const MeshCells *cells,
                          const MeshEntities *entities)
{
    Arena save = arena_save();

    idx_t *elmdist = arena_malloc(sync.size + 1, sizeof(*elmdist));
    elmdist[0] = 0;
    MPI_Allgather(&(idx_t){cells->num}, 1, IDX_T, &elmdist[1], 1, IDX_T, sync.comm);
    for (long i = 0; i < sync.size; i++) {
        elmdist[i + 1] += elmdist[i];
    }

    idx_t *eptr = arena_malloc(cells->num + 1, sizeof(*eptr));
    idx_t *eind = arena_malloc(cells->node.off[cells->num], sizeof(*eind));
    eptr[0] = 0;
    for (long i = 0; i < cells->num; i++) {
        eptr[i + 1] = eptr[i];
        for (long j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
            eind[eptr[i + 1]++] = cells->node.idx[j];
        }
    }

    idx_t numflag = 0;
    idx_t ncommon = 3;
    Dual dual;
    int ret = ParMETIS_V3_Mesh2Dual(elmdist, eptr, eind, &numflag, &ncommon, &dual.xadj,
                                    &dual.adjncy, &sync.comm);
    assert(ret == METIS_OK);

    for (long lhs = 0; lhs < entities->num; lhs++) {
        Name name;
        if (rotate_at_char(name, entities->name[lhs], ':') > 0) {
            for (long rhs = lhs + 1; rhs < entities->num; rhs++) {
                if (!strcmp(entities->name[rhs], name)) {
                    connect_periodic(nodes, cells, entities, &dual, lhs, rhs);
                }
            }
        }
    }

    arena_load(save);
    return dual;
}

/* Set partitions of outer cells to their adjacent inner cell's partition. */
static void collect_outer_parts(const MeshCells *cells, const Dual *dual, idx_t *part)
{
    Arena save = arena_save();

    long *num_cells = arena_malloc(sync.size + 1, sizeof(*num_cells));
    num_cells[0] = 0;
    MPI_Allgather(&cells->num, 1, MPI_LONG, &num_cells[1], 1, MPI_LONG, sync.comm);
    for (long i = 0; i < sync.size; i++) {
        num_cells[i + 1] += num_cells[i];
    }

    int *num_recv = arena_calloc(sync.size, sizeof(*num_recv));
    for (long i = cells->num_inner; i < cells->num; i++) {
        long inner = dual->adjncy[dual->xadj[i]];
        long rank = array_ldigitize(&num_cells[1], inner, sync.size);
        num_recv[rank] += 1;
    }

    int *off_recv = arena_malloc(sync.size + 1, sizeof(*off_recv));
    off_recv[0] = 0;
    for (long i = 0; i < sync.size; i++) {
        off_recv[i + 1] = off_recv[i] + num_recv[i];
    }

    long tot_recv = off_recv[sync.size];
    long *idx_recv = arena_malloc(tot_recv, sizeof(*idx_recv));
    for (long i = 0; i < sync.size; i++) {
        off_recv[i + 1] -= num_recv[i];
    }
    for (long i = cells->num_inner; i < cells->num; i++) {
        long inner = dual->adjncy[dual->xadj[i]];
        long rank = array_ldigitize(&num_cells[1], inner, sync.size);
        idx_recv[off_recv[rank + 1]++] = inner - num_cells[rank];
    }
    for (long i = 0; i < sync.size; i++) {
        assert(off_recv[i + 1] - off_recv[i] == num_recv[i]);
    }

    int *num_send = arena_malloc(sync.size, sizeof(*num_send));
    MPI_Alltoall(num_recv, 1, MPI_INT, num_send, 1, MPI_INT, sync.comm);

    int *off_send = arena_malloc(sync.size + 1, sizeof(*off_send));
    off_send[0] = 0;
    for (long i = 0; i < sync.size; i++) {
        off_send[i + 1] = off_send[i] + num_send[i];
    }

    long tot_send = off_send[sync.size];
    long *idx_send = arena_malloc(tot_send, sizeof(*idx_send));
    MPI_Alltoallv(idx_recv, num_recv, off_recv, MPI_LONG, idx_send, num_send, off_send, MPI_LONG,
                  sync.comm);

    idx_t *part_send = arena_malloc(tot_send, sizeof(*part_send));
    for (long i = 0; i < tot_send; i++) {
        part_send[i] = part[idx_send[i]];
    }

    idx_t *part_recv = arena_malloc(tot_recv, sizeof(*part_recv));
    MPI_Alltoallv(part_send, num_send, off_send, IDX_T, part_recv, num_recv, off_recv, IDX_T,
                  sync.comm);

    for (long i = 0; i < sync.size; i++) {
        off_recv[i + 1] -= num_recv[i];
    }
    for (long i = cells->num_inner; i < cells->num; i++) {
        long inner = dual->adjncy[dual->xadj[i]];
        long rank = array_ldigitize(&num_cells[1], inner, sync.size);
        part[i] = part_recv[off_recv[rank + 1]++];
    }
    for (long i = 0; i < sync.size; i++) {
        assert(off_recv[i + 1] - off_recv[i] == num_recv[i]);
    }

    arena_load(save);
}

/* Partition the dual graph (k-way) and refine with ParMETIS. */
static void compute_partitioning(const MeshCells *cells, const Dual *dual, idx_t *part)
{
    Arena save = arena_save();

    idx_t *vtxdist = arena_malloc(sync.size + 1, sizeof(*vtxdist));
    vtxdist[0] = 0;
    MPI_Allgather(&(idx_t){cells->num}, 1, IDX_T, &vtxdist[1], 1, IDX_T, sync.comm);
    for (long i = 0; i < sync.size; i++) {
        vtxdist[i + 1] += vtxdist[i];
    }

    idx_t wgtflag = 0;
    idx_t numflag = 0;

    idx_t ncon = 1;
    idx_t nparts = sync.size;
    real_t *tpwgts = arena_malloc(ncon * nparts, sizeof(*tpwgts));
    for (long i = 0; i < ncon * nparts; i++) {
        tpwgts[i] = 1.0 / nparts;
    }

    real_t *ubvec = arena_malloc(ncon, sizeof(*ubvec));
    for (long i = 0; i < ncon; i++) {
        ubvec[i] = 1.05;  // NOLINT(readability-magic-numbers)
    }

    idx_t options[4] = {0};
    idx_t edgecut;
    int ret =
        ParMETIS_V3_PartKway(vtxdist, dual->xadj, dual->adjncy, 0, 0, &wgtflag, &numflag, &ncon,
                             &nparts, tpwgts, ubvec, options, &edgecut, part, &sync.comm);
    assert(ret == METIS_OK);

    if (option.num_refines > 0) {
        verbose("Partition refinement:");
        verbose("\t %4s %12s %10s %s", "iter", "edgecut", "delta", "note");
        verbose("\t %4d %12ld %10d %s", 0, edgecut, 0, "initial");

        idx_t last_edgecut = edgecut;
        idx_t best_edgecut = edgecut;
        idx_t *best_part = arena_memdup(part, cells->num, sizeof(*best_part));
        for (long i = 0; i < option.num_refines; i++) {
            ret = ParMETIS_V3_RefineKway(vtxdist, dual->xadj, dual->adjncy, 0, 0, &wgtflag,
                                         &numflag, &ncon, &nparts, tpwgts, ubvec, options, &edgecut,
                                         part, &sync.comm);
            assert(ret == METIS_OK);

            const char *note = "same";
            if (edgecut < best_edgecut) {
                best_edgecut = edgecut;
                memcpy(best_part, part, cells->num * sizeof(*part));
                note = "improved";
            }
            else if (edgecut > best_edgecut) {
                edgecut = best_edgecut;
                memcpy(part, best_part, cells->num * sizeof(*part));
                note = "restored";
            }

            verbose("\t %4ld %12ld %+10ld %s", i + 1, edgecut, edgecut - last_edgecut, note);
            last_edgecut = edgecut;
        }
    }

    long *num_cells = arena_calloc(sync.size, sizeof(*num_cells));
    for (long i = 0; i < cells->num_inner; i++) {
        num_cells[part[i]] += 1;
    }
    MPI_Allreduce(MPI_IN_PLACE, num_cells, sync.size, MPI_LONG, MPI_SUM, sync.comm);
    assert(num_cells[sync.rank] > 0);

    collect_outer_parts(cells, dual, part);

    arena_load(save);
}

/* Gather remote adjacency ids, grouped by destination partition. */
static long *collect_adjncys(const MeshCells *cells, const Dual *dual, const idx_t *part,
                             long *num_adjncy)
{
    Arena save = arena_save();

    long tot_send = 0;
    int *num_send = arena_calloc(sync.size, sizeof(*num_send));
    for (long i = 0; i < cells->num; i++) {
        long num_cells = dual->xadj[i + 1] - dual->xadj[i];
        tot_send += num_cells;
        num_send[part[i]] += num_cells;
    }

    int *off_send = arena_malloc(sync.size + 1, sizeof(*off_send));
    off_send[0] = 0;
    for (long i = 0; i < sync.size; i++) {
        off_send[i + 1] = off_send[i] + num_send[i];
    }

    long *idx_send = arena_malloc(tot_send, sizeof(*idx_send));
    for (long i = 0; i < sync.size; i++) {
        off_send[i + 1] -= num_send[i];
    }
    for (long i = 0; i < cells->num; i++) {
        for (long j = dual->xadj[i]; j < dual->xadj[i + 1]; j++) {
            idx_send[off_send[part[i] + 1]++] = dual->adjncy[j];
        }
    }
    for (long i = 0; i < sync.size; i++) {
        assert(off_send[i + 1] - off_send[i] == num_send[i]);
    }

    int *num_recv = arena_malloc(sync.size, sizeof(*num_recv));
    MPI_Alltoall(num_send, 1, MPI_INT, num_recv, 1, MPI_INT, sync.comm);

    int *off_recv = arena_malloc(sync.size + 1, sizeof(*off_recv));
    off_recv[0] = 0;
    for (long i = 0; i < sync.size; i++) {
        off_recv[i + 1] = off_recv[i] + num_recv[i];
    }

    long tot_recv = off_recv[sync.size];
    long *idx_recv = arena_malloc(tot_recv, sizeof(*idx_recv));
    MPI_Alltoallv(idx_send, num_send, off_send, MPI_LONG, idx_recv, num_recv, off_recv, MPI_LONG,
                  sync.comm);

    array_lunique(idx_recv, &tot_recv);

    arena_load(save);

    *num_adjncy = tot_recv;
    return arena_smuggle(idx_recv, *num_adjncy, sizeof(*idx_recv));
}

/* Resolve partitions for unique remote adjacency ids. */
static void collect_parts(const MeshCells *cells, const idx_t *part, const long *adjncy,
                          long *adjncy_part, long num_adjncy)
{
    Arena save = arena_save();

    long *num_cells = arena_malloc(sync.size + 1, sizeof(*num_cells));
    num_cells[0] = 0;
    MPI_Allgather(&cells->num, 1, MPI_LONG, &num_cells[1], 1, MPI_LONG, sync.comm);
    for (long i = 0; i < sync.size; i++) {
        num_cells[i + 1] += num_cells[i];
    }

    int *num_recv = arena_calloc(sync.size, sizeof(*num_recv));
    for (long i = 0; i < num_adjncy; i++) {
        long rank = array_ldigitize(&num_cells[1], adjncy[i], sync.size);
        num_recv[rank] += 1;
    }

    int *off_recv = arena_malloc(sync.size + 1, sizeof(*off_recv));
    off_recv[0] = 0;
    for (long i = 0; i < sync.size; i++) {
        off_recv[i + 1] = off_recv[i] + num_recv[i];
    }

    long *idx_recv = arena_malloc(num_adjncy, sizeof(*idx_recv));
    for (long i = 0; i < sync.size; i++) {
        off_recv[i + 1] -= num_recv[i];
    }
    for (long i = 0; i < num_adjncy; i++) {
        long rank = array_ldigitize(&num_cells[1], adjncy[i], sync.size);
        idx_recv[off_recv[rank + 1]++] = adjncy[i] - num_cells[rank];
    }
    for (long i = 0; i < sync.size; i++) {
        assert(off_recv[i + 1] - off_recv[i] == num_recv[i]);
    }

    int *num_send = arena_malloc(sync.size, sizeof(*num_send));
    MPI_Alltoall(num_recv, 1, MPI_INT, num_send, 1, MPI_INT, sync.comm);

    int *off_send = arena_malloc(sync.size + 1, sizeof(*off_send));
    off_send[0] = 0;
    for (long i = 0; i < sync.size; i++) {
        off_send[i + 1] = off_send[i] + num_send[i];
    }

    long tot_send = off_send[sync.size];
    long *idx_send = arena_malloc(tot_send, sizeof(*idx_send));
    MPI_Alltoallv(idx_recv, num_recv, off_recv, MPI_LONG, idx_send, num_send, off_send, MPI_LONG,
                  sync.comm);

    long *part_send = arena_malloc(tot_send, sizeof(*part_send));
    for (long i = 0; i < tot_send; i++) {
        part_send[i] = part[idx_send[i]];
    }

    long *part_recv = arena_malloc(num_adjncy, sizeof(*part_recv));
    MPI_Alltoallv(part_send, num_send, off_send, MPI_LONG, part_recv, num_recv, off_recv, MPI_LONG,
                  sync.comm);

    for (long i = 0; i < sync.size; i++) {
        off_recv[i + 1] -= num_recv[i];
    }
    for (long i = 0; i < num_adjncy; i++) {
        long rank = array_ldigitize(&num_cells[1], adjncy[i], sync.size);
        adjncy_part[i] = part_recv[off_recv[rank + 1]++];
    }
    for (long i = 0; i < sync.size; i++) {
        assert(off_recv[i + 1] - off_recv[i] == num_recv[i]);
    }

    arena_load(save);
}

typedef struct {
    long num;
    long node[MAX_CELL_NODES];
} Neighbor;

/* Fetch node lists for remote inner neighbors. */
static Neighbor *collect_neighbor_cells(const MeshCells *cells, const long *adjncy,
                                        const long *adjncy_part, long num_adjncy,
                                        long *num_neighbors)
{
    Arena save = arena_save();

    long *num_cells = arena_malloc(sync.size + 1, sizeof(*num_cells));
    num_cells[0] = 0;
    MPI_Allgather(&cells->num, 1, MPI_LONG, &num_cells[1], 1, MPI_LONG, sync.comm);
    for (long i = 0; i < sync.size; i++) {
        num_cells[i + 1] += num_cells[i];
    }

    long *num_inner = arena_malloc(sync.size, sizeof(*num_inner));
    MPI_Allgather(&cells->num_inner, 1, MPI_LONG, num_inner, 1, MPI_LONG, sync.comm);

    int *num_recv = arena_calloc(sync.size, sizeof(*num_recv));
    for (long i = 0; i < num_adjncy; i++) {
        if (adjncy_part[i] != sync.rank) {
            long rank = array_ldigitize(&num_cells[1], adjncy[i], sync.size);
            if (adjncy[i] < num_cells[rank] + num_inner[rank]) {
                num_recv[rank] += 1;
            }
        }
    }

    int *off_recv = arena_malloc(sync.size + 1, sizeof(*off_recv));
    off_recv[0] = 0;
    for (long i = 0; i < sync.size; i++) {
        off_recv[i + 1] = off_recv[i] + num_recv[i];
    }

    long tot_recv = off_recv[sync.size];
    long *idx_recv = arena_malloc(tot_recv, sizeof(*idx_recv));
    for (long i = 0; i < sync.size; i++) {
        off_recv[i + 1] -= num_recv[i];
    }
    for (long i = 0; i < num_adjncy; i++) {
        if (adjncy_part[i] != sync.rank) {
            long rank = array_ldigitize(&num_cells[1], adjncy[i], sync.size);
            if (adjncy[i] < num_cells[rank] + num_inner[rank]) {
                idx_recv[off_recv[rank + 1]++] = adjncy[i] - num_cells[rank];
            }
        }
    }
    for (long i = 0; i < sync.size; i++) {
        assert(off_recv[i + 1] - off_recv[i] == num_recv[i]);
    }

    int *num_send = arena_malloc(sync.size, sizeof(*num_send));
    MPI_Alltoall(num_recv, 1, MPI_INT, num_send, 1, MPI_INT, sync.comm);

    int *off_send = arena_malloc(sync.size + 1, sizeof(*off_send));
    off_send[0] = 0;
    for (long i = 0; i < sync.size; i++) {
        off_send[i + 1] = off_send[i] + num_send[i];
    }

    long tot_send = off_send[sync.size];
    long *idx_send = arena_malloc(tot_send, sizeof(*idx_send));
    MPI_Alltoallv(idx_recv, num_recv, off_recv, MPI_LONG, idx_send, num_send, off_send, MPI_LONG,
                  sync.comm);

    Neighbor *send = arena_malloc(tot_send, sizeof(*send));
    for (long i = 0; i < tot_send; i++) {
        long idx = idx_send[i];
        send[i].num = cells->node.off[idx + 1] - cells->node.off[idx];
        for (long k = 0, j = cells->node.off[idx]; j < cells->node.off[idx + 1]; j++, k++) {
            send[i].node[k] = cells->node.idx[j];
        }
    }

    MPI_Datatype type;
    MPI_Type_contiguous(sizeof(Neighbor), MPI_BYTE, &type);
    MPI_Type_commit(&type);

    Neighbor *recv = arena_malloc(tot_recv, sizeof(*recv));
    MPI_Alltoallv(send, num_send, off_send, type, recv, num_recv, off_recv, type, sync.comm);
    MPI_Type_free(&type);

    arena_load(save);

    *num_neighbors = tot_recv;
    return arena_smuggle(recv, sync_lmax(*num_neighbors), sizeof(*recv));
}

/* Rebuild communicator as a dist-graph over adjacent partitions. */
static void decompose_comm(const long *adjncy_part, long num_adjncy)
{
    Arena save = arena_save();

    long *count = arena_calloc(sync.size, sizeof(*count));
    for (long i = 0; i < num_adjncy; i++) {
        if (adjncy_part[i] != sync.rank) {
            count[adjncy_part[i]] += 1;
        }
    }

    long deg = 0;
    for (long i = 0; i < sync.size; i++) {
        if (count[i] > 0) {
            deg += 1;
        }
    }

    int *rank = arena_malloc(deg, sizeof(*rank));
    int *weight = arena_malloc(deg, sizeof(*weight));
    long num = 0;
    for (long i = 0; i < sync.size; i++) {
        if (count[i] > 0) {
            rank[num] = i;
            weight[num] = count[i];
            num += 1;
        }
    }
    assert(num == deg);

    MPI_Comm comm;
    MPI_Dist_graph_create_adjacent(sync.comm, deg, rank, weight, deg, rank, weight, MPI_INFO_NULL,
                                   true, &comm);
    sync_reinit(comm);

    arena_load(save);
}

/* After comm change, swap node arrays and neighbor lists with previous rank. */
static void resolve_comm_reorder(MeshNodes *nodes, Neighbor *neighbor, long *num_neighbors, int dst)
{
    Arena save = arena_save();

    int *rank = arena_malloc(sync.size, sizeof(*rank));
    MPI_Allgather(&dst, 1, MPI_INT, rank, 1, MPI_INT, sync.comm);

    int src = -1;
    for (long i = 0; i < sync.size; i++) {
        if (rank[i] == sync.rank) {
            src = i;
            break;
        }
    }
    assert(src != -1);

    long num_nodes;
    int tag = sync_tag();
    MPI_Sendrecv(&nodes->num, 1, MPI_LONG, dst, tag, &num_nodes, 1, MPI_LONG, src, tag, sync.comm,
                 MPI_STATUS_IGNORE);

    MPI_Datatype type;
    MPI_Type_contiguous(sizeof(vector), MPI_BYTE, &type);
    MPI_Type_commit(&type);

    vector *coord = malloc(num_nodes * sizeof(*coord));
    tag = sync_tag();
    MPI_Sendrecv(nodes->coord, nodes->num, type, dst, tag, coord, num_nodes, type, src, tag,
                 sync.comm, MPI_STATUS_IGNORE);

    nodes->num = num_nodes;
    free(nodes->coord);
    nodes->coord = coord;

    MPI_Type_contiguous(sizeof(Neighbor), MPI_BYTE, &type);
    MPI_Type_commit(&type);

    int count = sync_lmax(*num_neighbors);
    tag = sync_tag();
    MPI_Sendrecv_replace(neighbor, count, type, dst, tag, src, tag, sync.comm, MPI_STATUS_IGNORE);
    MPI_Type_free(&type);

    tag = sync_tag();
    MPI_Sendrecv_replace(num_neighbors, 1, MPI_LONG, dst, tag, src, tag, sync.comm,
                         MPI_STATUS_IGNORE);

    arena_load(save);
}

typedef struct {
    long entity;
    long num;
    long node[MAX_CELL_NODES];
} Cell;

static int cmp_cell(const void *lhs_, const void *rhs_)
{
    const Cell *lhs = lhs_;
    const Cell *rhs = rhs_;
    return cmp_asc(lhs->entity, rhs->entity);
}

/* Move cells to their target partitions and rebuild entity offsets. */
static void redistribute_cells(MeshCells *cells, MeshEntities *entities, const idx_t *part)
{
    Arena save = arena_save();

    int *num_send = arena_calloc(sync.size, sizeof(*num_send));
    for (long i = 0; i < cells->num; i++) {
        num_send[part[i]] += 1;
    }

    int *off_send = arena_malloc(sync.size + 1, sizeof(*off_send));
    off_send[0] = 0;
    for (long i = 0; i < sync.size; i++) {
        off_send[i + 1] = off_send[i] + num_send[i];
    }

    Cell *send = arena_calloc(cells->num, sizeof(*send));
    for (long i = 0; i < sync.size; i++) {
        off_send[i + 1] -= num_send[i];
    }
    for (long entity = 0; entity < entities->num; entity++) {
        for (long i = entities->cell_off[entity]; i < entities->cell_off[entity + 1]; i++) {
            long idx = off_send[part[i] + 1]++;
            send[idx].entity = entity;
            send[idx].num = cells->node.off[i + 1] - cells->node.off[i];
            for (long k = 0, j = cells->node.off[i]; j < cells->node.off[i + 1]; j++, k++) {
                send[idx].node[k] = cells->node.idx[j];
            }
        }
    }
    for (long i = 0; i < sync.size; i++) {
        assert(off_send[i + 1] - off_send[i] == num_send[i]);
    }

    int *num_recv = arena_malloc(sync.size, sizeof(*num_recv));
    MPI_Alltoall(num_send, 1, MPI_INT, num_recv, 1, MPI_INT, sync.comm);

    int *off_recv = arena_malloc(sync.size + 1, sizeof(*off_recv));
    off_recv[0] = 0;
    for (long i = 0; i < sync.size; i++) {
        off_recv[i + 1] = off_recv[i] + num_recv[i];
    }

    MPI_Datatype type;
    MPI_Type_contiguous(sizeof(Cell), MPI_BYTE, &type);
    MPI_Type_commit(&type);

    cells->num = off_recv[sync.size];
    assert(cells->num > 0);

    Cell *recv = arena_malloc(cells->num, sizeof(*recv));
    MPI_Alltoallv(send, num_send, off_send, type, recv, num_recv, off_recv, type, sync.comm);
    MPI_Type_free(&type);

    qsort(recv, cells->num, sizeof(*recv), cmp_cell);

    long *off = realloc(cells->node.off, (cells->num + 1) * sizeof(*off));
    assert(off);

    long *idx = realloc(cells->node.idx, (cells->num * MAX_CELL_NODES) * sizeof(*idx));
    assert(idx);

    memset(entities->cell_off, 0, (entities->num + 1) * sizeof(*entities->cell_off));

    for (long i = 0; i < cells->num; i++) {
        off[i + 1] = off[i] + recv[i].num;
        for (long k = 0, j = off[i]; j < off[i + 1]; j++, k++) {
            idx[j] = recv[i].node[k];
        }
        entities->cell_off[recv[i].entity + 1] += 1;
    }
    for (long i = 0; i < entities->num; i++) {
        entities->cell_off[i + 1] += entities->cell_off[i];
    }

    cells->node.off = off;
    cells->node.idx = realloc(idx, off[cells->num] * sizeof(*idx));
    assert(cells->node.idx);

    arena_load(save);
}

/* Append received neighbor cells to local cells. */
static void append_neighbor_cells(MeshCells *cells, const Neighbor *neighbor, long num_neighbors)
{
    long num_cells = cells->num + num_neighbors;
    assert(num_cells > 0);

    long *off = realloc(cells->node.off, (num_cells + 1) * sizeof(*off));
    assert(off);

    long *idx = realloc(cells->node.idx, (num_cells * MAX_CELL_NODES) * sizeof(*idx));
    assert(idx);

    long num = 0;
    for (long i = cells->num; i < num_cells; i++) {
        off[i + 1] = off[i];
        for (long j = 0; j < neighbor[num].num; j++) {
            idx[off[i + 1]++] = neighbor[num].node[j];
        }
        num += 1;
    }
    assert(num == num_neighbors);

    cells->num = num_cells;
    cells->node.off = off;
    cells->node.idx = realloc(idx, off[cells->num] * sizeof(*idx));
    assert(cells->node.idx);
}

/* Partition cells (build dual, partition, comm rebuild), then append neighbors. */
static void partition_cells(MeshNodes *nodes, MeshCells *cells, MeshEntities *entities)
{
    Arena save = arena_save();

    Dual dual = connect_cells(nodes, cells, entities);

    idx_t *part = arena_malloc(cells->num, sizeof(*part));
    compute_partitioning(cells, &dual, part);

    long num_adjncy;
    long *adjncy = collect_adjncys(cells, &dual, part, &num_adjncy);

    long *adjncy_part = arena_malloc(num_adjncy, sizeof(*adjncy_part));
    collect_parts(cells, part, adjncy, adjncy_part, num_adjncy);

    long num_neighbors;
    Neighbor *neighbor =
        collect_neighbor_cells(cells, adjncy, adjncy_part, num_adjncy, &num_neighbors);

    int dst = sync.rank;
    decompose_comm(adjncy_part, num_adjncy);
    free(dual.xadj);
    free(dual.adjncy);

    if (dst != sync.rank) {
        resolve_comm_reorder(nodes, neighbor, &num_neighbors, dst);
    }

    redistribute_cells(cells, entities, part);
    append_neighbor_cells(cells, neighbor, num_neighbors);

    arena_load(save);
}

/* Derive {inner,ghost,periodic} counts from entity offsets. */
static void compute_cell_counts(MeshCells *cells, const MeshEntities *entities)
{
    long num_inner = 0;
    long num_ghost = 0;
    long num_periodic = 0;
    for (long i = 0; i < entities->num; i++) {
        if (i < entities->num_inner) {
            num_inner += entities->cell_off[i + 1] - entities->cell_off[i];
        }
        else if (i < entities->off_ghost) {
            num_ghost += entities->cell_off[i + 1] - entities->cell_off[i];
        }
        else {
            num_periodic += entities->cell_off[i + 1] - entities->cell_off[i];
        }
    }
    cells->num_inner = num_inner;
    cells->off_ghost = num_inner + num_ghost;
    cells->off_periodic = num_inner + num_ghost + num_periodic;
}

typedef struct {
    long rank;
    long old;
    long new;
} Node;

static int cmp_node(const void *lhs_, const void *rhs_)
{
    const Node *lhs = lhs_;
    const Node *rhs = rhs_;
    return cmp_asc(lhs->rank, rhs->rank);
}

/* Redistribute nodes to match current cells and renumber connectivity. */
static void partition_nodes(MeshNodes *nodes, MeshCells *cells)
{
    Arena save = arena_save();

    long *global = arena_malloc(cells->node.off[cells->num], sizeof(*global));
    long num = 0;
    for (long i = 0; i < cells->num; i++) {
        for (long j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
            global[num++] = cells->node.idx[j];
        }
    }
    array_lunique(global, &num);

    vector *coord = arena_malloc(num, sizeof(*coord));
    collect_coords(nodes, global, coord, num);

    long cap = sync_lmax(num);
    Node *node = arena_calloc(cap, sizeof(*node));
    for (long i = 0; i < num; i++) {
        node[i].rank = sync.rank;
        node[i].old = global[i];
        node[i].new = -1;
    }
    for (long i = num; i < cap; i++) {
        node[i].old = node[i].new = -1;
    }

    MPI_Datatype type;
    MPI_Type_contiguous(sizeof(Node), MPI_BYTE, &type);
    MPI_Type_commit(&type);

    int dst = (sync.rank + 1) % sync.size;
    int src = (sync.rank - 1 + sync.size) % sync.size;
    int tag = sync_tag();
    for (long step = 0; step < sync.size; step++) {
        for (long i = 0; i < cap; i++) {
            long key = node[i].old;
            if (key != -1 && bsearch(&key, global, num, sizeof(*global), cmp_long)) {
                node[i].rank = lmin(node[i].rank, sync.rank);
            }
        }
        MPI_Sendrecv_replace(node, cap, type, dst, tag, src, tag, sync.comm, MPI_STATUS_IGNORE);
    }

    long num_inner = 0;
    for (long i = 0; i < num; i++) {
        if (node[i].rank == sync.rank) {
            node[i].rank = -sync.rank;
            num_inner += 1;
        }
    }
    qsort(node, num, sizeof(*node), cmp_node);

    long (*old2new)[2] = arena_malloc(num_inner, sizeof(*old2new));
    long off_nodes = sync_lexsum(num_inner);
    long idx = 0;
    for (long i = 0; i < num; i++) {
        if (node[i].rank == -sync.rank) {
            node[i].new = i + off_nodes;
            old2new[idx][0] = node[i].old;
            old2new[idx][1] = node[i].new;
            idx += 1;
        }
    }
    assert(idx == num_inner);
    qsort(old2new, num_inner, sizeof(*old2new), cmp_long);

    tag = sync_tag();
    for (long step = 0; step < sync.size; step++) {
        for (long i = 0; i < cap; i++) {
            if (node[i].new == -1) {
                long key = node[i].old;
                long *val = bsearch(&key, old2new, num_inner, sizeof(*old2new), cmp_long);
                if (val) {
                    node[i].new = val[1];
                }
            }
        }
        MPI_Sendrecv_replace(node, cap, type, dst, tag, src, tag, sync.comm, MPI_STATUS_IGNORE);
    }
    MPI_Type_free(&type);

    nodes->num = num;
    assert(nodes->num > 0);
    nodes->num_inner = num_inner;

    free(nodes->coord);
    nodes->coord = malloc(nodes->num * sizeof(*nodes->coord));
    assert(nodes->coord);

    for (long i = 0; i < nodes->num; i++) {
        long key = node[i].old;
        long *val = bsearch(&key, global, num, sizeof(*global), cmp_long);
        assert(val);
        nodes->coord[i] = coord[val - global];
    }

    for (long i = 0; i < nodes->num; i++) {
        assert(node[i].new != -1);
        global[i] = node[i].new;
    }

    long (*old2local)[2] = arena_malloc(nodes->num, sizeof(*old2local));
    for (long i = 0; i < nodes->num; i++) {
        old2local[i][0] = node[i].old;
        old2local[i][1] = i;
    }
    qsort(old2local, nodes->num, sizeof(*old2local), cmp_long);

    for (long i = 0; i < cells->num; i++) {
        for (long j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
            long key = cells->node.idx[j];
            long *val = bsearch(&key, old2local, nodes->num, sizeof(*old2local), cmp_long);
            assert(val);
            cells->node.idx[j] = val[1];
        }
    }

    arena_load(save);

    nodes->global = arena_smuggle(global, nodes->num, sizeof(*global));
}

/* Translation between periodic pair = mean(center_rhs) - mean(center_lhs). */
static void compute_translation(const MeshNodes *nodes, const MeshCells *cells,
                                const MeshEntities *entities, vector *translation, long lhs,
                                long rhs)
{
    long num_cells[2] = {0};
    vector mean[2] = {0};
    for (long i = 0; i < entities->num; i++) {
        if (i == lhs || i == rhs) {
            for (long j = entities->cell_off[i]; j < entities->cell_off[i + 1]; j++) {
                vector center = {0};
                long num_nodes = cells->node.off[j + 1] - cells->node.off[j];
                for (long k = cells->node.off[j]; k < cells->node.off[j + 1]; k++) {
                    vector coord = nodes->coord[cells->node.idx[k]];
                    center = vector_add(center, vector_div(coord, num_nodes));
                }
                Side side = (i == lhs) ? LHS : RHS;
                num_cells[side] += 1;
                mean[side] = vector_add(mean[side], center);
            }
        }
    }

    num_cells[LHS] = sync_lsum(num_cells[LHS]);
    num_cells[RHS] = sync_lsum(num_cells[RHS]);
    assert(num_cells[LHS] == num_cells[RHS]);

    mean[LHS] = vector_div(sync_vsum(mean[LHS]), num_cells[LHS]);
    mean[RHS] = vector_div(sync_vsum(mean[RHS]), num_cells[RHS]);

    translation[lhs] = vector_sub(mean[RHS], mean[LHS]);
    translation[rhs] = vector_sub(mean[LHS], mean[RHS]);
}

/* Compute offsets for all periodic entity pairs. */
static void compute_translations(const MeshNodes *nodes, const MeshCells *cells,
                                 MeshEntities *entities)
{
    vector *translation = arena_calloc(entities->num, sizeof(*translation));
    for (long lhs = 0; lhs < entities->num; lhs++) {
        Name name;
        if (rotate_at_char(name, entities->name[lhs], ':') > 0) {
            for (long rhs = lhs + 1; rhs < entities->num; rhs++) {
                if (!strcmp(entities->name[rhs], name)) {
                    compute_translation(nodes, cells, entities, translation, lhs, rhs);
                }
            }
        }
    }
    entities->translation = translation;
}

typedef struct {
    vector center;
    long entity;
    long rank;
} Recv;

/* Identify owner ranks of periodic outers by matching shifted centers. */
static void collect_periodic_ranks(const MeshCells *cells, const MeshEntities *entities, Recv *recv,
                                   long cap)
{
    Arena save = arena_save();

    Kdtree *centers = kdtree_create(0);
    long num = 0;
    for (long i = entities->off_ghost; i < entities->num; i++) {
        for (long j = entities->cell_off[i]; j < entities->cell_off[i + 1]; j++) {
            vector center = vector_add(recv[num++].center, entities->translation[i]);
            kdtree_insert(centers, center, 0);
        }
    }
    assert(num == cells->off_periodic - cells->off_ghost);

    for (long i = 0; i < num; i++) {
        recv[i].rank = -1;
    }

    MPI_Datatype type;
    MPI_Type_contiguous(sizeof(Recv), MPI_BYTE, &type);
    MPI_Type_commit(&type);

    int dst = (sync.rank + 1) % sync.size;
    int src = (sync.rank - 1 + sync.size) % sync.size;
    int tag = sync_tag();
    for (long step = 0; step < sync.size; step++) {
        for (long i = 0; i < cap; i++) {
            if (recv[i].rank == -1 && kdtree_lookup(centers, recv[i].center)) {
                recv[i].rank = sync.rank;
            }
        }
        MPI_Sendrecv_replace(recv, cap, type, dst, tag, src, tag, sync.comm, MPI_STATUS_IGNORE);
    }
    MPI_Type_free(&type);

    for (long i = 0; i < num; i++) {
        assert(recv[i].rank != -1);
    }

    arena_load(save);
}

/* Identify owner ranks of non-periodic outers by matching inner centers. */
static void collect_neighbor_ranks(const MeshNodes *nodes, const MeshCells *cells, Recv *recv,
                                   long cap, long num)
{
    Arena save = arena_save();

    Kdtree *centers = kdtree_create(0);
    for (long i = 0; i < cells->num_inner; i++) {
        vector center = {0};
        long num_nodes = cells->node.off[i + 1] - cells->node.off[i];
        for (long j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
            vector coord = nodes->coord[cells->node.idx[j]];
            center = vector_add(center, vector_div(coord, num_nodes));
        }
        kdtree_insert(centers, center, 0);
    }

    for (long i = cells->off_periodic - cells->off_ghost; i < num; i++) {
        recv[i].rank = -1;
    }

    MPI_Datatype type;
    MPI_Type_contiguous(sizeof(Recv), MPI_BYTE, &type);
    MPI_Type_commit(&type);

    int dst = (sync.rank + 1) % sync.size;
    int src = (sync.rank - 1 + sync.size) % sync.size;
    int tag = sync_tag();
    for (long step = 0; step < sync.size; step++) {
        for (long i = 0; i < cap; i++) {
            if (recv[i].rank == -1 && kdtree_lookup(centers, recv[i].center)) {
                recv[i].rank = sync.rank;
            }
        }
        MPI_Sendrecv_replace(recv, cap, type, dst, tag, src, tag, sync.comm, MPI_STATUS_IGNORE);
    }
    MPI_Type_free(&type);

    for (long i = cells->off_periodic - cells->off_ghost; i < num; i++) {
        assert(recv[i].rank != -1);
    }

    arena_load(save);
}

static void compute_cell_map(const MeshCells *cells, const MeshEntities *entities, const Recv *recv,
                             long tot, long *map)
{
    Arena save = arena_save();

    long (*off)[sync.size + 1] = arena_calloc(entities->num + 1, sizeof(*off));
    long num = 0;
    for (long i = entities->off_ghost; i < entities->num; i++) {
        for (long j = entities->cell_off[i]; j < entities->cell_off[i + 1]; j++) {
            off[i][recv[num++].rank + 1] += 1;
        }
    }
    assert(num == cells->off_periodic - cells->off_ghost);
    for (long i = num; i < tot; i++) {
        off[entities->num][recv[i].rank + 1] += 1;
    }

    long base = 0;
    for (long i = 0; i < entities->num + 1; ++i) {
        for (int j = 0; j < sync.size; ++j) {
            off[i][j + 1] += off[i][j];
        }
        for (int j = 0; j < sync.size + 1; ++j) {
            off[i][j] += base;
        }
        base = off[i][sync.size];
    }

    num = 0;
    for (long i = entities->off_ghost; i < entities->num; i++) {
        for (long j = entities->cell_off[i]; j < entities->cell_off[i + 1]; j++) {
            map[num] = off[i][recv[num].rank]++;
            num += 1;
        }
    }
    assert(num == cells->off_periodic - cells->off_ghost);
    for (long i = num; i < tot; i++) {
        map[i] = off[entities->num][recv[i].rank]++;
    }

    arena_load(save);
}

/* Reorder outer cells and Recv by (entity,rank). */
static void reorder(MeshCells *cells, const MeshEntities *entities, Recv *recv, long beg, long end)
{
    Arena save = arena_save();

    long tot = end - beg;
    long *map = arena_malloc(tot, sizeof(*map));
    compute_cell_map(cells, entities, recv, tot, map);

    mesh_reorder_cells(cells, 0, beg, end, map);

    Recv *buf = arena_malloc(tot, sizeof(*buf));
    for (long i = 0; i < tot; i++) {
        buf[map[i]] = recv[i];
    }
    memcpy(recv, buf, tot * sizeof(*buf));

    arena_load(save);
}

/* Build neighbor groups from outers sorted by (entity,rank). */
static void create_neighbors(const MeshNodes *nodes, MeshCells *cells, const MeshEntities *entities,
                             MeshNeighbors *neighbors)
{
    Arena save = arena_save();

    long beg = cells->off_ghost;
    long end = cells->num;
    if (beg == end) {
        return;
    }

    long tot = end - beg;
    long cap = sync_lmax(tot);
    Recv *recv = arena_calloc(cap, sizeof(*recv));

    long num = 0;
    for (long i = beg; i < end; i++) {
        vector center = {0};
        long num_nodes = cells->node.off[i + 1] - cells->node.off[i];
        for (long j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
            vector coord = nodes->coord[cells->node.idx[j]];
            center = vector_add(center, vector_div(coord, num_nodes));
        }
        recv[num].center = center;
        recv[num].entity = array_ldigitize(&entities->cell_off[1], i, entities->num);
        num += 1;
    }
    assert(num == tot);

    collect_periodic_ranks(cells, entities, recv, cap);
    collect_neighbor_ranks(nodes, cells, recv, cap, num);

    reorder(cells, entities, recv, beg, end);

    neighbors->num = 1;
    for (long i = 1; i < num; i++) {
        if (recv[i].entity != recv[i - 1].entity || recv[i].rank != recv[i - 1].rank) {
            neighbors->num += 1;
        }
    }

    long *rank = arena_malloc(neighbors->num, sizeof(*rank));
    long *recv_off = arena_malloc(neighbors->num + 1, sizeof(*recv_off));

    rank[0] = recv[0].rank;
    recv_off[0] = beg;

    long idx = 1;
    for (long i = 1; i < num; i++) {
        if (recv[i].entity != recv[i - 1].entity || recv[i].rank != recv[i - 1].rank) {
            rank[idx] = recv[i].rank;
            recv_off[idx] = beg + i;
            idx += 1;
        }
    }
    recv_off[idx] = end;
    assert(idx == neighbors->num);

    arena_load(save);

    neighbors->rank = arena_smuggle(rank, neighbors->num, sizeof(*rank));
    neighbors->recv_off = arena_smuggle(recv_off, neighbors->num + 1, sizeof(*recv_off));
}

/* Convert heap-backed arrays to arena-owned arrays. */
static void convert_allocations(MeshNodes *nodes, MeshCells *cells, MeshEntities *entities)
{
    vector *coord = arena_memdup(nodes->coord, nodes->num, sizeof(*coord));
    free(nodes->coord);
    nodes->coord = coord;

    long *off = arena_memdup(cells->node.off, cells->num + 1, sizeof(*off));
    free(cells->node.off);
    cells->node.off = off;

    long *idx = arena_memdup(cells->node.idx, off[cells->num], sizeof(*idx));
    free(cells->node.idx);
    cells->node.idx = idx;

    Name *name = arena_memdup(entities->name, entities->num, sizeof(*name));
    free(entities->name);
    entities->name = name;

    long *cell_off = arena_memdup(entities->cell_off, entities->num + 1, sizeof(*cell_off));
    free(entities->cell_off);
    entities->cell_off = cell_off;
}

Mesh *mesh_read(const char *fname)
{
    assert(fname);

    Mesh *mesh = arena_calloc(1, sizeof(*mesh));

    read_file(mesh, fname);

    mesh->cells.num_inner = count_inner_cells(&mesh->entities);
    assert(sync.size <= sync_lsum(mesh->cells.num_inner));

    partition_cells(&mesh->nodes, &mesh->cells, &mesh->entities);
    compute_cell_counts(&mesh->cells, &mesh->entities);

    partition_nodes(&mesh->nodes, &mesh->cells);

    compute_translations(&mesh->nodes, &mesh->cells, &mesh->entities);
    create_neighbors(&mesh->nodes, &mesh->cells, &mesh->entities, &mesh->neighbors);

    convert_allocations(&mesh->nodes, &mesh->cells, &mesh->entities);

    return mesh;
}
