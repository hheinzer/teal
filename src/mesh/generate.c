#include <assert.h>
#include <math.h>
#include <metis.h>
#include <stdlib.h>

#include "mesh.h"
#include "reorder.h"
#include "teal/arena.h"
#include "teal/array.h"
#include "teal/kdtree.h"
#include "teal/sync.h"
#include "teal/utils.h"
#include "teal/vector.h"

// Build dual (cell->cell) graph from cell->node with METIS; drop outer-outer edges.
static void connect_cells(const MeshNodes *nodes, MeshCells *cells)
{
    Arena save = arena_save();

    idx_t num_elems = cells->num;
    idx_t num_nodes = nodes->num;

    idx_t *eptr = arena_malloc(num_elems + 1, sizeof(*eptr));
    idx_t *eind = arena_malloc(num_elems * (MAX_CELL_NODES + 1), sizeof(*eind));

    eptr[0] = 0;

    for (long i = 0; i < cells->num; i++) {
        eptr[i + 1] = eptr[i];
        for (long j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
            eind[eptr[i + 1]++] = cells->node.idx[j];
        }
        if (i >= cells->num_inner) {
            eind[eptr[i + 1]++] = num_nodes++;  // add dummy node to enforce ncommon=3
        }
    }

    idx_t ncommon = 3;
    idx_t numflag = 0;

    idx_t *xadj;
    idx_t *adjncy;
    long ret =
        METIS_MeshToDual(&num_elems, &num_nodes, eptr, eind, &ncommon, &numflag, &xadj, &adjncy);
    assert(ret == METIS_OK);

    long *off = arena_malloc(cells->num + 1, sizeof(*off));
    long *idx = arena_malloc(cells->num * MAX_CELL_FACES, sizeof(*idx));

    off[0] = 0;

    for (long i = 0; i < cells->num; i++) {
        off[i + 1] = off[i];
        for (long j = xadj[i]; j < xadj[i + 1]; j++) {
            if (i < cells->num_inner || adjncy[j] < cells->num_inner) {
                idx[off[i + 1]++] = adjncy[j];
            }
        }
        assert(off[i + 1] - off[i] > 0);
    }
    free(xadj);
    free(adjncy);

    arena_load(save);

    cells->cell.off = arena_smuggle(off, cells->num + 1, sizeof(*off));
    cells->cell.idx = arena_smuggle(idx, cells->cell.off[cells->num], sizeof(*idx));
}

// Pick the next unmapped inner cell with the smallest geometric center.
static long find_seed_cell(const MeshNodes *nodes, const MeshCells *cells, const long *map)
{
    long seed = -1;
    vector min_center = {SCALAR_MAX, SCALAR_MAX, SCALAR_MAX};
    for (long i = 0; i < cells->num_inner; i++) {
        if (map[i] == -1) {
            vector center = {0};
            long num_nodes = cells->node.off[i + 1] - cells->node.off[i];
            for (long j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
                vector coord = nodes->coord[cells->node.idx[j]];
                vector_inc(&center, vector_div(coord, num_nodes));
            }
            if (cmp_vector(&center, &min_center) < 0) {
                seed = i;
                min_center = center;
            }
        }
    }
    assert(seed != -1);
    return seed;
}

// BFS-reorder inner cells [0,num_inner); remap cell->cell accordingly.
static void improve_cell_ordering(const MeshNodes *nodes, MeshCells *cells)
{
    Arena save = arena_save();

    long *map = arena_malloc(cells->num_inner, sizeof(*map));
    for (long i = 0; i < cells->num_inner; i++) {
        map[i] = -1;
    }

    long *queue = arena_malloc(cells->num_inner, sizeof(*queue));

    long num = 0;
    while (num < cells->num_inner) {
        long seed = find_seed_cell(nodes, cells, map);
        long beg = 0;
        long end = 0;
        queue[end++] = seed;
        map[seed] = num++;
        while (beg < end) {
            long cur = queue[beg++];
            for (long i = cells->cell.off[cur]; i < cells->cell.off[cur + 1]; i++) {
                long idx = cells->cell.idx[i];
                if (idx < cells->num_inner && map[idx] == -1) {
                    queue[end++] = idx;
                    map[idx] = num++;
                }
            }
        }
    }
    assert(num == cells->num_inner);

    mesh_reorder_cells(cells, 0, 0, cells->num_inner, map);

    arena_load(save);
}

// Collect new global node indices for outer nodes from owning ranks.
static void collect_global(long *global, const long *map, long num_inner, long num)
{
    Arena save = arena_save();

    long *num_nodes = arena_malloc(sync.size + 1, sizeof(*num_nodes));
    num_nodes[0] = 0;
    MPI_Allgather(&num_inner, 1, MPI_LONG, &num_nodes[1], 1, MPI_LONG, sync.comm);
    for (long i = 0; i < sync.size; i++) {
        num_nodes[i + 1] += num_nodes[i];
    }

    int *num_recv = arena_calloc(sync.size, sizeof(*num_recv));
    for (long i = num_inner; i < num; i++) {
        long rank = array_digitize(&num_nodes[1], global[i], sync.size);
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
    for (long i = num_inner; i < num; i++) {
        long rank = array_digitize(&num_nodes[1], global[i], sync.size);
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

    long *global_send = arena_malloc(tot_send, sizeof(*global_send));
    for (long i = 0; i < tot_send; i++) {
        global_send[i] = global[map[idx_send[i]]];
    }

    long tot_recv = off_recv[sync.size];
    long *global_recv = arena_malloc(tot_recv, sizeof(*global_recv));
    MPI_Alltoallv(global_send, num_send, off_send, MPI_LONG, global_recv, num_recv, off_recv,
                  MPI_LONG, sync.comm);

    for (long i = 0; i < sync.size; i++) {
        off_recv[i + 1] -= num_recv[i];
    }
    for (long i = num_inner; i < num; i++) {
        long rank = array_digitize(&num_nodes[1], global[i], sync.size);
        global[i] = global_recv[off_recv[rank + 1]++];
    }
    for (long i = 0; i < sync.size; i++) {
        assert(off_recv[i + 1] - off_recv[i] == num_recv[i]);
    }

    arena_load(save);
}

// Reorder nodes so inner nodes used by cells come first, then outer nodes.
static void improve_node_ordering(MeshNodes *nodes, MeshCells *cells)
{
    Arena save = arena_save();

    long *map = arena_malloc(nodes->num, sizeof(*map));
    for (long i = 0; i < nodes->num; i++) {
        map[i] = -1;
    }

    long num = 0;
    for (long i = 0; i < cells->num; i++) {
        for (long j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
            if (cells->node.idx[j] < nodes->num_inner) {
                if (map[cells->node.idx[j]] == -1) {
                    map[cells->node.idx[j]] = num++;
                }
            }
        }
    }
    assert(num == nodes->num_inner);

    for (long i = 0; i < cells->num; i++) {
        for (long j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
            if (cells->node.idx[j] >= nodes->num_inner) {
                if (map[cells->node.idx[j]] == -1) {
                    map[cells->node.idx[j]] = num++;
                }
            }
        }
    }
    assert(num == nodes->num);

    mesh_reorder_nodes(nodes, cells, map);

    long off_inner = sync_exsum(nodes->num_inner);
    for (long i = 0; i < nodes->num_inner; i++) {
        nodes->global[i] = off_inner + i;
    }

    collect_global(nodes->global, map, nodes->num_inner, nodes->num);

    arena_load(save);
}

// Fill `node` with node indices common to cells left and right and return the count.
static long compute_intersection(long *node, const MeshCells *cells, long left, long right)
{
    long num = 0;
    for (long i = cells->node.off[left]; i < cells->node.off[left + 1]; i++) {
        for (long j = cells->node.off[right]; j < cells->node.off[right + 1]; j++) {
            if (cells->node.idx[i] == cells->node.idx[j]) {
                assert(num < MAX_FACE_NODES);
                node[num++] = cells->node.idx[i];
            }
        }
    }
    return num;
}

// Build faces: one per unique adjacent cell pair; fill face->node and face->cell.
static void create_faces(const MeshCells *cells, MeshFaces *faces)
{
    Arena save = arena_save();

    long cap = cells->num_inner * MAX_CELL_FACES;
    struct {
        long left;
        long right;
        long num;
        long node[MAX_FACE_NODES];
    } *face = arena_malloc(cap, sizeof(*face));

    long num = 0;
    for (long i = 0; i < cells->num_inner; i++) {
        for (long j = cells->cell.off[i]; j < cells->cell.off[i + 1]; j++) {
            long left = i;
            long right = cells->cell.idx[j];
            if (left < right) {
                assert(num < cap);
                face[num].left = left;
                face[num].right = right;
                face[num].num = compute_intersection(face[num].node, cells, left, right);
                num += 1;
            }
        }
    }

    long *off = arena_malloc(num + 1, sizeof(*off));
    long *idx = arena_malloc(num * MAX_FACE_NODES, sizeof(*idx));
    Adjacent *cell = arena_malloc(num, sizeof(*cell));

    off[0] = 0;

    long num_inner = 0;
    long num_ghost = 0;
    for (long i = 0; i < num; i++) {
        assert(face[i].left < cells->num_inner);
        if (face[i].right < cells->num_inner) {
            num_inner += 1;
        }
        else if (face[i].right < cells->off_ghost) {
            num_ghost += 1;
        }
        off[i + 1] = off[i] + face[i].num;
        for (long k = 0, j = off[i]; j < off[i + 1]; j++, k++) {
            idx[j] = face[i].node[k];
        }
        cell[i].left = face[i].left;
        cell[i].right = face[i].right;
    }
    assert(num_ghost == cells->off_ghost - cells->num_inner);

    arena_load(save);

    faces->num = num;
    faces->num_inner = num_inner;
    faces->off_ghost = num_inner + num_ghost;
    faces->node.off = arena_smuggle(off, num + 1, sizeof(*off));
    faces->node.idx = arena_smuggle(idx, faces->node.off[faces->num], sizeof(*idx));
    faces->cell = arena_smuggle(cell, num, sizeof(*cell));
}

// Reorder faces: inner first (by left), then outer (by right).
static void reorder(const MeshCells *cells, MeshFaces *faces)
{
    Arena save = arena_save();

    long *key = arena_malloc(faces->num, sizeof(*key));
    for (long i = 0; i < faces->num; i++) {
        if (faces->cell[i].right < cells->num_inner) {
            key[i] = faces->cell[i].left;
        }
        else {
            key[i] = faces->cell[i].right;
        }
    }
    mesh_reorder_faces(faces, key);

    arena_load(save);
}

// Compute per-entity face offsets: inner faces first, then one face per outer cell by entity.
static void compute_face_entities(const MeshFaces *faces, MeshEntities *entities)
{
    long *face_off = arena_malloc(entities->num + 1, sizeof(*face_off));
    face_off[0] = faces->num_inner;
    for (long i = 0; i < entities->num; i++) {
        face_off[i + 1] = face_off[i];
        if (i >= entities->num_inner) {
            face_off[i + 1] += entities->cell_off[i + 1] - entities->cell_off[i];
        }
    }
    entities->face_off = face_off;
}

// Build send lists: match neighbors' requested receive centers to local cells.
static void compute_send_graph(const MeshNodes *nodes, const MeshCells *cells,
                               const MeshEntities *entities, MeshNeighbors *neighbors)
{
    Arena save = arena_save();

    Kdtree *center2local = kdtree_create(sizeof(long));
    for (long i = cells->off_ghost; i < cells->num; i++) {
        vector center = {0};
        long num_nodes = cells->node.off[i + 1] - cells->node.off[i];
        for (long j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
            vector coord = nodes->coord[cells->node.idx[j]];
            vector_inc(&center, vector_div(coord, num_nodes));
        }
        long local = cells->cell.idx[cells->cell.off[i]];
        kdtree_insert(center2local, center, &local);
    }
    for (long i = cells->off_periodic; i < cells->num; i++) {
        for (long j = cells->cell.off[i]; j < cells->cell.off[i + 1]; j++) {
            long local = cells->cell.idx[j];
            vector center = {0};
            long num_nodes = cells->node.off[local + 1] - cells->node.off[local];
            for (long k = cells->node.off[local]; k < cells->node.off[local + 1]; k++) {
                vector coord = nodes->coord[cells->node.idx[k]];
                vector_inc(&center, vector_div(coord, num_nodes));
            }
            kdtree_insert(center2local, center, &local);
        }
    }

    long tot_recv = cells->num - cells->off_ghost;
    vector *recv = arena_malloc(tot_recv, sizeof(*recv));
    long num = 0;
    for (long i = cells->off_ghost; i < cells->num; i++) {
        vector center = {0};
        long num_nodes = cells->node.off[i + 1] - cells->node.off[i];
        for (long j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
            vector coord = nodes->coord[cells->node.idx[j]];
            vector_inc(&center, vector_div(coord, num_nodes));
        }
        long entity = array_digitize(&entities->cell_off[1], i, entities->num);
        if (entity < entities->num) {
            recv[num] = vector_add(center, entities->translation[entity]);
        }
        else {
            recv[num] = center;
        }
        num += 1;
    }
    assert(num == tot_recv);

    long *off = arena_malloc(neighbors->num + 1, sizeof(*off));
    off[0] = 0;
    long tag = sync_tag();
    MPI_Request *req_recv = arena_malloc(neighbors->num, sizeof(*req_recv));
    MPI_Request *req_send = arena_malloc(neighbors->num, sizeof(*req_send));
    for (long i = 0; i < neighbors->num; i++) {
        long num_recv = neighbors->recv_off[i + 1] - neighbors->recv_off[i];
        MPI_Irecv(&off[i + 1], 1, MPI_LONG, neighbors->rank[i], tag, sync.comm, &req_recv[i]);
        MPI_Isend(&num_recv, 1, MPI_LONG, neighbors->rank[i], tag, sync.comm, &req_send[i]);
    }
    MPI_Waitall(neighbors->num, req_recv, MPI_STATUSES_IGNORE);
    MPI_Waitall(neighbors->num, req_send, MPI_STATUSES_IGNORE);
    for (long i = 0; i < neighbors->num; i++) {
        off[i + 1] += off[i];
    }

    MPI_Datatype type;
    MPI_Type_contiguous(sizeof(vector), MPI_BYTE, &type);
    MPI_Type_commit(&type);

    long tot_send = off[neighbors->num];
    vector *send = arena_malloc(tot_send, sizeof(*send));
    tag = sync_tag();
    long off_recv = 0;
    for (long i = 0; i < neighbors->num; i++) {
        long num_recv = neighbors->recv_off[i + 1] - neighbors->recv_off[i];
        long num_send = off[i + 1] - off[i];
        MPI_Irecv(&send[off[i]], num_send, type, neighbors->rank[i], tag, sync.comm, &req_recv[i]);
        MPI_Isend(&recv[off_recv], num_recv, type, neighbors->rank[i], tag, sync.comm,
                  &req_send[i]);
        off_recv += num_recv;
    }
    assert(off_recv == tot_recv);
    MPI_Waitall(neighbors->num, req_recv, MPI_STATUSES_IGNORE);
    MPI_Waitall(neighbors->num, req_send, MPI_STATUSES_IGNORE);
    MPI_Type_free(&type);

    long *idx = arena_malloc(tot_send, sizeof(*idx));
    for (long i = 0; i < tot_send; i++) {
        long *local = kdtree_lookup(center2local, send[i]);
        assert(local);
        idx[i] = *local;
    }

    arena_load(save);

    neighbors->send.off = arena_smuggle(off, neighbors->num + 1, sizeof(*off));
    neighbors->send.idx = arena_smuggle(idx, tot_send, sizeof(*idx));
}

// Ensure quad vertex order yields a valid planar quad; triangles are always valid.
static void correct_coord_order(vector *coord, long num_nodes)
{
    switch (num_nodes) {
        case 3: return;
        case 4:
            for (long i = 0; i < 3; i++) {
                vector a2b = vector_sub(coord[1], coord[0]);
                vector a2c = vector_sub(coord[2], coord[0]);
                vector a2d = vector_sub(coord[3], coord[0]);
                if (vector_dot(vector_cross(a2b, a2c), vector_cross(a2c, a2d)) > 0) {
                    return;
                }
                if (i + 2 < num_nodes) {
                    swap_vector(&coord[i + 1], &coord[i + 2]);
                }
            }
            error("quad vertex order could not be made valid");
        default: error("invalid number of nodes (%ld)", num_nodes);
    }
}

// Area of a triangular or quadrilateral face.
static scalar compute_face_area(const vector *coord, long num_nodes)
{
    switch (num_nodes) {
        case 3: {
            vector a2b = vector_sub(coord[1], coord[0]);
            vector a2c = vector_sub(coord[2], coord[0]);
            return vector_norm(vector_cross(a2b, a2c)) / 2;
        }
        case 4: {
            vector lhs[3] = {coord[0], coord[1], coord[2]};
            vector rhs[3] = {coord[0], coord[2], coord[3]};
            return compute_face_area(lhs, 3) + compute_face_area(rhs, 3);
        }
        default: error("invalid number of nodes (%ld)", num_nodes);
    }
}

// Weighted average of vectors with scalar weights.
static vector weighted_average(const vector *arr, const scalar *wgt, long num)
{
    vector wsum = {0};
    for (long i = 0; i < num; i++) {
        vector_inc(&wsum, vector_mul(wgt[i], arr[i]));
    }
    return vector_div(wsum, array_fsum(wgt, num));
}

// Face centroid for triangles/quads (quad via area-weighted triangles).
static vector compute_face_center(const vector *coord, long num_nodes)
{
    switch (num_nodes) {
        case 3: return vector_div(vector_sum(coord, 3), 3);
        case 4: {
            vector lhs[3] = {coord[0], coord[1], coord[2]};
            vector rhs[3] = {coord[0], coord[2], coord[3]};
            vector cen[2] = {compute_face_center(lhs, 3), compute_face_center(rhs, 3)};
            scalar area[2] = {compute_face_area(lhs, 3), compute_face_area(rhs, 3)};
            return weighted_average(cen, area, 2);
        }
        default: error("invalid number of nodes (%ld)", num_nodes);
    }
}

// Unit normal for triangles/quads.
static vector compute_face_normal(const vector *coord, long num_nodes)
{
    switch (num_nodes) {
        case 3: {
            vector a2b = vector_sub(coord[1], coord[0]);
            vector a2c = vector_sub(coord[2], coord[0]);
            vector normal = vector_cross(a2b, a2c);
            return vector_div(normal, vector_norm(normal));
        }
        case 4: {
            vector normal = {0};
            for (long i = 0; i < 4; i++) {
                vector_inc(&normal, vector_cross(coord[i], coord[(i + 1) % 4]));
            }
            return vector_div(normal, vector_norm(normal));
        }
        default: error("invalid number of nodes (%ld)", num_nodes);
    }
}

// Orthonormal basis on a face aligned with its normal.
static Basis compute_face_basis(const vector *coord, long num_nodes)
{
    Basis basis;
    vector normal = basis.normal = compute_face_normal(coord, num_nodes);
    scalar nqz = hypot(normal.x, normal.y);
    scalar nqy = hypot(normal.x, normal.z);
    if (nqz > nqy) {
        basis.tangent1 = (vector){-normal.y / nqz, normal.x / nqz, 0};
        basis.tangent2 = (vector){-normal.x * normal.z / nqz, -normal.y * normal.z / nqz, nqz};
    }
    else {
        basis.tangent1 = (vector){-normal.x * normal.y / nqy, nqy, -normal.y * normal.z / nqy};
        basis.tangent2 = (vector){-normal.z / nqy, 0, normal.x / nqy};
    }
    return basis;
}

// Flip the face normal so it points from the left cell to the right cell if necessary.
static void correct_face_basis(const MeshNodes *nodes, const MeshCells *cells, long left,
                               vector center, Basis *basis)
{
    vector mean = {0};
    long num_nodes = cells->node.off[left + 1] - cells->node.off[left];
    for (long i = cells->node.off[left]; i < cells->node.off[left + 1]; i++) {
        vector coord = nodes->coord[cells->node.idx[i]];
        vector_inc(&mean, vector_div(coord, num_nodes));
    }
    if (vector_dot(vector_sub(mean, center), basis->normal) > 0) {
        vector_scale(&basis->normal, -1);
        vector_scale(&basis->tangent1, -1);
    }
}

static void compute_face_geometry(const MeshNodes *nodes, const MeshCells *cells, MeshFaces *faces)
{
    scalar *area = arena_malloc(faces->num, sizeof(*area));
    vector *center = arena_malloc(faces->num, sizeof(*center));
    Basis *basis = arena_malloc(faces->num, sizeof(*basis));

    for (long i = 0; i < faces->num; i++) {
        vector coord[MAX_FACE_NODES] = {0};
        long num_nodes = faces->node.off[i + 1] - faces->node.off[i];
        for (long k = 0, j = faces->node.off[i]; j < faces->node.off[i + 1]; j++, k++) {
            coord[k] = nodes->coord[faces->node.idx[j]];
        }
        correct_coord_order(coord, num_nodes);
        area[i] = compute_face_area(coord, num_nodes);
        center[i] = compute_face_center(coord, num_nodes);
        basis[i] = compute_face_basis(coord, num_nodes);
        correct_face_basis(nodes, cells, faces->cell[i].left, center[i], &basis[i]);
    }

    faces->area = area;
    faces->center = center;
    faces->basis = basis;
}

// Volume of supported cell types by tetrahedral decomposition.
static scalar compute_cell_volume(const vector *coord, long num_nodes)
{
    switch (num_nodes) {
        case 4: {
            vector a2b = vector_sub(coord[1], coord[0]);
            vector a2c = vector_sub(coord[2], coord[0]);
            vector a2d = vector_sub(coord[3], coord[0]);
            return fabs(vector_dot(a2b, vector_cross(a2c, a2d))) / 6;
        }
        case 5: {
            vector lhs[4] = {coord[0], coord[1], coord[2], coord[4]};
            vector rhs[4] = {coord[0], coord[2], coord[3], coord[4]};
            return compute_cell_volume(lhs, 4) + compute_cell_volume(rhs, 4);
        }
        case 6: {
            vector lhs[5] = {coord[0], coord[1], coord[4], coord[3], coord[5]};
            vector rhs[4] = {coord[0], coord[1], coord[2], coord[5]};
            return compute_cell_volume(lhs, 5) + compute_cell_volume(rhs, 4);
        }
        case 8: {
            vector lhs[6] = {coord[0], coord[1], coord[2], coord[4], coord[5], coord[6]};
            vector rhs[6] = {coord[0], coord[2], coord[3], coord[4], coord[6], coord[7]};
            return compute_cell_volume(lhs, 6) + compute_cell_volume(rhs, 6);
        }
        default: error("invalid number of nodes (%ld)", num_nodes);
    }
}

// Cell centroid, volume-weighted for mixed splits.
static vector compute_cell_center(const vector *coord, long num_nodes)
{
    switch (num_nodes) {
        case 4: return vector_div(vector_sum(coord, num_nodes), 4);
        case 5: {
            vector lhs[4] = {coord[0], coord[1], coord[2], coord[4]};
            vector rhs[4] = {coord[0], coord[2], coord[3], coord[4]};
            vector cen[2] = {compute_cell_center(lhs, 4), compute_cell_center(rhs, 4)};
            scalar vol[2] = {compute_cell_volume(lhs, 4), compute_cell_volume(rhs, 4)};
            return weighted_average(cen, vol, 2);
        }
        case 6: {
            vector lhs[5] = {coord[0], coord[1], coord[4], coord[3], coord[5]};
            vector rhs[4] = {coord[0], coord[1], coord[2], coord[5]};
            vector cen[2] = {compute_cell_center(lhs, 5), compute_cell_center(rhs, 4)};
            scalar vol[2] = {compute_cell_volume(lhs, 5), compute_cell_volume(rhs, 4)};
            return weighted_average(cen, vol, 2);
        }
        case 8: {
            vector lhs[6] = {coord[0], coord[1], coord[2], coord[4], coord[5], coord[6]};
            vector rhs[6] = {coord[0], coord[2], coord[3], coord[4], coord[6], coord[7]};
            vector cen[2] = {compute_cell_center(lhs, 6), compute_cell_center(rhs, 6)};
            scalar vol[2] = {compute_cell_volume(lhs, 6), compute_cell_volume(rhs, 6)};
            return weighted_average(cen, vol, 2);
        }
        default: error("invalid number of nodes (%ld)", num_nodes);
    }
}

// Exchange centers for outer cells so neighbor/periodic copies receive an interior center.
static void collect_centers(const MeshNeighbors *neighbors, vector *center)
{
    Arena save = arena_save();

    MPI_Datatype type;
    MPI_Type_contiguous(sizeof(vector), MPI_BYTE, &type);
    MPI_Type_commit(&type);

    vector *send = arena_malloc(neighbors->send.off[neighbors->num], sizeof(*send));
    long tag = sync_tag();
    MPI_Request *req_recv = arena_malloc(neighbors->num, sizeof(*req_recv));
    MPI_Request *req_send = arena_malloc(neighbors->num, sizeof(*req_send));
    for (long i = 0; i < neighbors->num; i++) {
        for (long j = neighbors->send.off[i]; j < neighbors->send.off[i + 1]; j++) {
            send[j] = center[neighbors->send.idx[j]];
        }
        long sendcount = neighbors->send.off[i + 1] - neighbors->send.off[i];
        long recvcount = neighbors->recv_off[i + 1] - neighbors->recv_off[i];
        MPI_Irecv(&center[neighbors->recv_off[i]], recvcount, type, neighbors->rank[i], tag,
                  sync.comm, &req_recv[i]);
        MPI_Isend(&send[neighbors->send.off[i]], sendcount, type, neighbors->rank[i], tag,
                  sync.comm, &req_send[i]);
    }
    MPI_Waitall(neighbors->num, req_recv, MPI_STATUSES_IGNORE);
    MPI_Waitall(neighbors->num, req_send, MPI_STATUSES_IGNORE);
    MPI_Type_free(&type);

    arena_load(save);
}

// Compute cell volumes, centers, and axis-aligned projections.
static void compute_cell_geometry(const MeshNodes *nodes, MeshCells *cells, const MeshFaces *faces,
                                  const MeshEntities *entities, const MeshNeighbors *neighbors)
{
    scalar *volume = arena_calloc(cells->num, sizeof(*volume));
    vector *center = arena_malloc(cells->num, sizeof(*center));
    vector *projection = arena_calloc(cells->num, sizeof(*projection));

    for (long i = 0; i < cells->num_inner; i++) {
        vector coord[MAX_CELL_NODES] = {0};
        long num_nodes = cells->node.off[i + 1] - cells->node.off[i];
        for (long k = 0, j = cells->node.off[i]; j < cells->node.off[i + 1]; j++, k++) {
            coord[k] = nodes->coord[cells->node.idx[j]];
        }
        volume[i] = compute_cell_volume(coord, num_nodes);
        center[i] = compute_cell_center(coord, num_nodes);
    }

    for (long i = 0; i < faces->num; i++) {
        long left = faces->cell[i].left;
        long right = faces->cell[i].right;
        vector normal = faces->basis[i].normal;
        vector inc = vector_mul(faces->area[i] / 2, vector_abs(normal));
        vector_inc(&projection[left], inc);
        if (right < cells->num_inner) {
            vector_inc(&projection[right], inc);
        }
        else if (i < faces->off_ghost) {  // mirror left center across the faces
            vector c2f = vector_sub(faces->center[i], center[left]);
            center[right] =
                vector_add(center[left], vector_mul(2 * vector_dot(c2f, normal), normal));
        }
    }

    collect_centers(neighbors, center);

    for (long i = entities->off_ghost; i < entities->num; i++) {
        for (long j = entities->cell_off[i]; j < entities->cell_off[i + 1]; j++) {
            vector_dec(&center[j], entities->translation[i]);
        }
    }

    cells->volume = volume;
    cells->center = center;
    cells->projection = projection;
    cells->sum_volume = sync_fsum(array_fsum(volume, cells->num_inner));
}

void compute_cell_offsets(const MeshNodes *nodes, MeshCells *cells)
{
    vector *offset = arena_malloc(cells->cell.off[cells->num_inner], sizeof(*offset));
    for (long i = 0; i < cells->num_inner; i++) {
        for (long j = cells->cell.off[i]; j < cells->cell.off[i + 1]; j++) {
            long node[MAX_FACE_NODES];
            long num_nodes = compute_intersection(node, cells, i, cells->cell.idx[j]);

            vector coord[MAX_FACE_NODES];
            for (long k = 0; k < num_nodes; k++) {
                coord[k] = nodes->coord[node[k]];
            }
            correct_coord_order(coord, num_nodes);

            vector center = compute_face_center(coord, num_nodes);
            offset[j] = vector_sub(center, cells->center[i]);
        }
    }
    cells->offset = offset;
}

static vector compute_weight(const scalar *r11, const scalar *r12, const scalar *r22,
                             const scalar *r13, const scalar *r23, const scalar *r33, vector delta,
                             long idx)
{
    vector weight = {0};
    scalar beta = (r12[idx] * r23[idx] - r13[idx] * r22[idx]) / (r11[idx] * r22[idx]);
    vector alpha;
    alpha.x = delta.x / sq(r11[idx]);
    alpha.y = (delta.y - r12[idx] / r11[idx] * delta.x) / sq(r22[idx]);
    alpha.z = (delta.z - r23[idx] / r22[idx] * delta.y + beta * delta.x) / sq(r33[idx]);
    if (!is_close(delta.x, 0)) {
        weight.x += alpha.x;
    }
    if (!is_close(delta.y, 0)) {
        weight.x += -r12[idx] / r11[idx] * alpha.y;
        weight.y += alpha.y;
    }
    if (!is_close(delta.z, 0)) {
        weight.x += beta * alpha.z;
        weight.y += -r23[idx] / r22[idx] * alpha.z;
        weight.z += alpha.z;
    }
    scalar theta2 = 1 / vector_norm2(delta);
    vector_scale(&weight, theta2);
    assert(isfinite(weight.x));
    assert(isfinite(weight.y));
    assert(isfinite(weight.z));
    return weight;
}

// Least-squares weights for gradient reconstruction on faces.
static void compute_face_weights(const MeshCells *cells, MeshFaces *faces)
{
    Arena save = arena_save();

    scalar *r11 = arena_calloc(cells->num_inner, sizeof(*r11));
    scalar *r12 = arena_calloc(cells->num_inner, sizeof(*r12));
    scalar *r22 = arena_calloc(cells->num_inner, sizeof(*r22));
    scalar *r13 = arena_calloc(cells->num_inner, sizeof(*r13));
    scalar *r23 = arena_calloc(cells->num_inner, sizeof(*r23));
    scalar *r33 = arena_calloc(cells->num_inner, sizeof(*r33));
    for (long i = 0; i < faces->num; i++) {
        long left = faces->cell[i].left;
        long right = faces->cell[i].right;
        vector delta = vector_sub(cells->center[right], cells->center[left]);
        scalar theta2 = 1 / vector_norm2(delta);
        r11[left] += theta2 * delta.x * delta.x;
        r12[left] += theta2 * delta.x * delta.y;
        r22[left] += theta2 * delta.y * delta.y;
        r13[left] += theta2 * delta.x * delta.z;
        r23[left] += theta2 * delta.y * delta.z;
        r33[left] += theta2 * delta.z * delta.z;
        if (right < cells->num_inner) {
            r11[right] += theta2 * delta.x * delta.x;
            r12[right] += theta2 * delta.x * delta.y;
            r22[right] += theta2 * delta.y * delta.y;
            r13[right] += theta2 * delta.x * delta.z;
            r23[right] += theta2 * delta.y * delta.z;
            r33[right] += theta2 * delta.z * delta.z;
        }
    }
    for (long i = 0; i < cells->num_inner; i++) {
        r11[i] = sqrt(r11[i]);
        r12[i] = r12[i] / r11[i];
        r22[i] = sqrt(r22[i] - sq(r12[i]));
        r13[i] = r13[i] / r11[i];
        r23[i] = (r23[i] - r12[i] * r13[i]) / r22[i];
        r33[i] = sqrt(r33[i] - (sq(r13[i]) + sq(r23[i])));
    }

    Weight *weight = arena_malloc(faces->num, sizeof(*weight));
    for (long i = 0; i < faces->num; i++) {
        long left = faces->cell[i].left;
        long right = faces->cell[i].right;
        vector delta_l = vector_sub(cells->center[right], cells->center[left]);
        weight[i].left = compute_weight(r11, r12, r22, r13, r23, r33, delta_l, left);
        if (right < cells->num_inner) {
            vector delta_r = vector_sub(cells->center[left], cells->center[right]);
            weight[i].right = compute_weight(r11, r12, r22, r13, r23, r33, delta_r, right);
        }
    }

    arena_load(save);

    faces->weight = arena_smuggle(weight, faces->num, sizeof(*weight));
}

// Compute face-to-cell offset vectors from face center to each adjacent cell center.
static void compute_face_offsets(const MeshCells *cells, MeshFaces *faces)
{
    Offset *offset = arena_malloc(faces->num, sizeof(*offset));
    for (long i = 0; i < faces->num; i++) {
        long left = faces->cell[i].left;
        long right = faces->cell[i].right;
        offset[i].left = vector_sub(faces->center[i], cells->center[left]);
        offset[i].right = vector_sub(faces->center[i], cells->center[right]);
    }
    faces->offset = offset;
}

// Compute face correction vectors (unit direction and length between adjacent cell centers).
static void compute_face_correction(const MeshCells *cells, MeshFaces *faces)
{
    Correction *correction = arena_malloc(faces->num, sizeof(*correction));
    for (long i = 0; i < faces->num; i++) {
        long left = faces->cell[i].left;
        long right = faces->cell[i].right;
        vector delta = vector_sub(cells->center[right], cells->center[left]);
        scalar norm = vector_norm(delta);
        assert(norm > 0);
        correction[i].unit = vector_div(delta, norm);
        correction[i].norm = norm;
    }
    faces->correction = correction;
}

void mesh_generate(Mesh *mesh)
{
    assert(mesh);

    connect_cells(&mesh->nodes, &mesh->cells);

    improve_cell_ordering(&mesh->nodes, &mesh->cells);
    improve_node_ordering(&mesh->nodes, &mesh->cells);

    create_faces(&mesh->cells, &mesh->faces);
    reorder(&mesh->cells, &mesh->faces);

    compute_face_entities(&mesh->faces, &mesh->entities);

    compute_send_graph(&mesh->nodes, &mesh->cells, &mesh->entities, &mesh->neighbors);

    compute_face_geometry(&mesh->nodes, &mesh->cells, &mesh->faces);
    compute_cell_geometry(&mesh->nodes, &mesh->cells, &mesh->faces, &mesh->entities,
                          &mesh->neighbors);

    compute_cell_offsets(&mesh->nodes, &mesh->cells);
    compute_face_weights(&mesh->cells, &mesh->faces);
    compute_face_offsets(&mesh->cells, &mesh->faces);
    compute_face_correction(&mesh->cells, &mesh->faces);
}
