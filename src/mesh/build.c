#include <math.h>
#include <metis.h>
#include <stdlib.h>

#include "mesh.h"
#include "teal/arena.h"
#include "teal/array.h"
#include "teal/dict.h"
#include "teal/kdtree.h"
#include "teal/sync.h"
#include "teal/utils.h"
#include "teal/vector.h"

/* Build cell-to-cell connectivity from cell-to-node connectivity using METIS. */
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
    int ret =
        METIS_MeshToDual(&num_elems, &num_nodes, eptr, eind, &ncommon, &numflag, &xadj, &adjncy);
    assert(ret == METIS_OK);

    arena_load(save);

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
        assert(off[i + 1] - off[i] > 0);  // FIXME: rank has cells that do not belong here
    }
    free(xadj);
    free(adjncy);

    cells->cell.off = off;
    cells->cell.idx = arena_resize(idx, off[cells->num], sizeof(*idx));
    assert(cells->cell.idx);
}

/* Build face-to-node and face-to-cell connectivity. Each unique adjacent cell pair is one face. */
static void create_faces(const MeshCells *cells, MeshFaces *faces)
{
    Arena save = arena_save();

    typedef struct {
        long left;
        long right;
    } Pair;

    typedef struct {
        long num;
        long node[MAX_FACE_NODES];
    } Face;

    Dict pair2face = dict_create(sizeof(Pair), sizeof(Face));
    for (long i = 0; i < cells->num; i++) {
        for (long j = cells->cell.off[i]; j < cells->cell.off[i + 1]; j++) {
            long min = lmin(i, cells->cell.idx[j]);
            long max = lmax(i, cells->cell.idx[j]);
            Face face = {0};
            for (long ii = cells->node.off[min]; ii < cells->node.off[min + 1]; ii++) {
                for (long jj = cells->node.off[max]; jj < cells->node.off[max + 1]; jj++) {
                    if (cells->node.idx[ii] == cells->node.idx[jj]) {
                        face.node[face.num++] = cells->node.idx[ii];  // intersection of cell nodes
                    }
                }
            }
            dict_insert(&pair2face, &(Pair){min, max}, &face);
        }
    }

    long num_faces = pair2face.num;
    long *off = arena_malloc(num_faces + 1, sizeof(*off));
    long *idx = arena_malloc(num_faces * MAX_FACE_NODES, sizeof(*idx));
    MeshFaceCell *cell = arena_malloc(num_faces, sizeof(*cell));

    off[0] = 0;

    long num = 0;
    long num_inner = 0;
    long num_ghost = 0;
    for (DictItem *item = pair2face.beg; item; item = item->next) {
        Pair *pair = item->key;
        Face *face = item->val;
        assert(pair->left < cells->num_inner);
        if (pair->right < cells->num_inner) {
            num_inner += 1;
        }
        else if (pair->right < cells->num_inner + cells->num_ghost) {
            num_ghost += 1;
        }
        off[num + 1] = off[num] + face->num;
        for (long j = 0, i = off[num]; i < off[num + 1]; i++, j++) {
            idx[i] = face->node[j];
        }
        cell[num].left = pair->left;
        cell[num].right = pair->right;
        num += 1;
    }
    assert(num == num_faces);
    assert(num_ghost == cells->num_ghost);

    arena_load(save);

    faces->num = num_faces;
    faces->num_inner = num_inner;
    faces->num_ghost = num_ghost;
    faces->node.off = arena_smuggle(off, num_faces + 1, sizeof(*off));
    faces->node.idx = arena_smuggle(idx, faces->node.off[faces->num], sizeof(*idx));
    faces->cell = arena_smuggle(cell, num_faces, sizeof(*cell));
}

static long *compute_face_map(const MeshCells *cells, const MeshFaces *faces)
{
    long *map = arena_malloc(faces->num, sizeof(*map));
    for (long i = 0; i < faces->num; i++) {
        if (faces->cell[i].right < cells->num_inner) {
            map[i] = faces->cell[i].left;
        }
        else {
            map[i] = faces->cell[i].right;
        }
    }
    return map;
}

/* Reorder face arrays to [inner (sorted by left), outer (sorted by right)]. */
static void reorder_faces(const MeshCells *cells, MeshFaces *faces)
{
    Arena save = arena_save();

    typedef struct {
        long map;
        long num;
        long node[MAX_FACE_NODES];
        MeshFaceCell cell;
    } Face;

    Face *face = arena_malloc(faces->num, sizeof(*face));

    long *map = compute_face_map(cells, faces);
    for (long i = 0; i < faces->num; i++) {
        face[i].map = map[i];
        face[i].num = faces->node.off[i + 1] - faces->node.off[i];
        for (long k = 0, j = faces->node.off[i]; j < faces->node.off[i + 1]; j++, k++) {
            face[i].node[k] = faces->node.idx[j];
        }
        face[i].cell = faces->cell[i];
    }
    qsort(face, faces->num, sizeof(*face), lcmp);

    for (long i = 0; i < faces->num; i++) {
        faces->node.off[i + 1] = faces->node.off[i] + face[i].num;
        for (long k = 0, j = faces->node.off[i]; j < faces->node.off[i + 1]; j++, k++) {
            faces->node.idx[j] = face[i].node[k];
        }
        faces->cell[i] = face[i].cell;
    }

    arena_load(save);
}

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

/* Build neighbor send graph by matching outer-cell centers against neighbor receive windows. */
static void compute_send_graph(const MeshNodes *nodes, const MeshCells *cells,
                               const MeshEntities *entities, MeshNeighbors *neighbors)
{
    Arena save = arena_save();

    Kdtree center2local = kdtree_create(sizeof(long));
    for (long i = 0; i < cells->num; i++) {
        if (cells->num_inner <= i && i < cells->num_inner + cells->num_ghost) {
            continue;  // skip ghost cells
        }
        vector center = {0};
        long num_nodes = cells->node.off[i + 1] - cells->node.off[i];
        for (long j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
            vector coord = nodes->coord[cells->node.idx[j]];
            center = vector_add(center, vector_div(coord, num_nodes));
        }
        long local = (i < cells->num_inner) ? i : cells->cell.idx[cells->cell.off[i]];
        kdtree_insert(&center2local, center, &local);
    }

    long tot_recv = cells->num - cells->num_inner - cells->num_ghost;
    vector *recv = arena_malloc(tot_recv, sizeof(*recv));
    long num = 0;
    for (long i = cells->num_inner + cells->num_ghost; i < cells->num; i++) {
        vector center = {0};
        long num_nodes = cells->node.off[i + 1] - cells->node.off[i];
        for (long j = cells->node.off[i]; j < cells->node.off[i + 1]; j++) {
            vector coord = nodes->coord[cells->node.idx[j]];
            center = vector_add(center, vector_div(coord, num_nodes));
        }
        long entity = array_digitize(&entities->cell_off[1], i, entities->num);
        if (entity < entities->num) {
            recv[num] = vector_add(center, entities->offset[entity]);
        }
        else {
            recv[num] = center;
        }
        num += 1;
    }
    assert(num == tot_recv);

    long *off = arena_malloc(neighbors->num + 1, sizeof(*off));
    off[0] = 0;
    int tag = sync_tag();
    MPI_Request *req = arena_malloc(neighbors->num, sizeof(*req));
    for (long i = 0; i < neighbors->num; i++) {
        long num_recv = neighbors->recv_off[i + 1] - neighbors->recv_off[i];
        MPI_Isendrecv(&num_recv, 1, MPI_LONG, neighbors->rank[i], tag, &off[i + 1], 1, MPI_LONG,
                      neighbors->rank[i], tag, sync.comm, &req[i]);
    }
    MPI_Waitall(neighbors->num, req, MPI_STATUSES_IGNORE);
    for (long i = 0; i < neighbors->num; i++) {
        off[i + 1] += off[i];
    }

    long tot_send = off[neighbors->num];
    vector *send = arena_malloc(tot_send, sizeof(*send));
    tag = sync_tag();
    long off_recv = 0;
    for (long i = 0; i < neighbors->num; i++) {
        long num_recv = neighbors->recv_off[i + 1] - neighbors->recv_off[i];
        long num_send = off[i + 1] - off[i];
        MPI_Isendrecv(&recv[off_recv], num_recv, vector_type, neighbors->rank[i], tag,
                      &send[off[i]], num_send, vector_type, neighbors->rank[i], tag, sync.comm,
                      &req[i]);
        off_recv += num_recv;
    }
    MPI_Waitall(neighbors->num, req, MPI_STATUSES_IGNORE);
    assert(off_recv == tot_recv);

    long *idx = arena_malloc(tot_send, sizeof(*idx));
    for (long i = 0; i < tot_send; i++) {
        long *local = kdtree_lookup(&center2local, send[i]);
        assert(local);
        idx[i] = *local;
    }

    arena_load(save);

    neighbors->send.off = arena_smuggle(off, neighbors->num + 1, sizeof(*off));
    neighbors->send.idx = arena_smuggle(idx, tot_send, sizeof(*idx));
}

/* Ensure quad vertex order yields a valid planar quad; triangles are always valid. */
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
                    vswap(&coord[i + 1], &coord[i + 2]);
                }
            }
            abort();
        default: abort();
    }
}

static double compute_face_area(const vector *coord, long num_nodes)
{
    switch (num_nodes) {
        case 3: {
            vector a2b = vector_sub(coord[1], coord[0]);
            vector a2c = vector_sub(coord[2], coord[0]);
            return vector_len(vector_cross(a2b, a2c)) / 2;
        }
        case 4: {
            vector lhs[3] = {coord[0], coord[1], coord[2]};
            vector rhs[3] = {coord[0], coord[2], coord[3]};
            return compute_face_area(lhs, 3) + compute_face_area(rhs, 3);
        }
        default: abort();
    }
}

static vector weighted_average(const vector *arr, const double *wgt, long num)
{
    vector wsum = {0};
    for (long i = 0; i < num; i++) {
        wsum = vector_add(wsum, vector_mul(arr[i], wgt[i]));
    }
    return vector_div(wsum, array_fsum(wgt, num));
}

static vector compute_face_center(const vector *coord, long num_nodes)
{
    switch (num_nodes) {
        case 3: return vector_div(array_vsum(coord, 3), 3);
        case 4: {
            vector lhs[3] = {coord[0], coord[1], coord[2]};
            vector rhs[3] = {coord[0], coord[2], coord[3]};
            vector cen[2] = {compute_face_center(lhs, 3), compute_face_center(rhs, 3)};
            double area[2] = {compute_face_area(lhs, 3), compute_face_area(rhs, 3)};
            return weighted_average(cen, area, 2);
        }
        default: abort();
    }
}

static vector compute_face_normal(const vector *coord, long num_nodes)
{
    switch (num_nodes) {
        case 3: {
            vector a2b = vector_sub(coord[1], coord[0]);
            vector a2c = vector_sub(coord[2], coord[0]);
            return vector_unit(vector_cross(a2b, a2c));
        }
        case 4: {
            vector sum = {0};
            for (long i = 0; i < 4; i++) {
                sum = vector_add(sum, vector_cross(coord[i], coord[(i + 1) % 4]));
            }
            return vector_unit(sum);
        }
        default: abort();
    }
}

static matrix compute_face_basis(const vector *coord, long num_nodes)
{
    matrix basis;
    vector nrm = basis.x = compute_face_normal(coord, num_nodes);
    double nqz = hypot(nrm.x, nrm.y);
    double nqy = hypot(nrm.x, nrm.z);
    if (nqz > nqy) {
        basis.y = (vector){-nrm.y / nqz, nrm.x / nqz, 0};
        basis.z = (vector){-nrm.x * nrm.z / nqz, -nrm.y * nrm.z / nqz, nqz};
    }
    else {
        basis.y = (vector){-nrm.x * nrm.y / nqy, nqy, -nrm.y * nrm.z / nqy};
        basis.z = (vector){-nrm.z / nqy, 0, nrm.x / nqy};
    }
    return basis;
}

/* Flip the face normal so it points from the left cell to the right cell if necessary. */
static void correct_face_basis(const MeshNodes *nodes, const MeshCells *cells, long left,
                               vector center, matrix *basis)
{
    vector mean = {0};
    long num_nodes = cells->node.off[left + 1] - cells->node.off[left];
    for (long i = cells->node.off[left]; i < cells->node.off[left + 1]; i++) {
        vector coord = nodes->coord[cells->node.idx[i]];
        mean = vector_add(mean, vector_div(coord, num_nodes));
    }
    if (vector_dot(vector_sub(mean, center), basis->x) > 0) {
        basis->x = vector_mul(basis->x, -1);
        basis->y = vector_mul(basis->y, -1);
    }
}

static void compute_face_geometry(const MeshNodes *nodes, const MeshCells *cells, MeshFaces *faces)
{
    double *area = arena_malloc(faces->num, sizeof(*area));
    vector *center = arena_malloc(faces->num, sizeof(*center));
    matrix *basis = arena_malloc(faces->num, sizeof(*basis));

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

static double compute_cell_volume(const vector *coord, long num_nodes)
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
        default: abort();
    }
}

static vector compute_cell_center(const vector *coord, long num_nodes)
{
    switch (num_nodes) {
        case 4: return vector_div(array_vsum(coord, num_nodes), 4);
        case 5: {
            vector lhs[4] = {coord[0], coord[1], coord[2], coord[4]};
            vector rhs[4] = {coord[0], coord[2], coord[3], coord[4]};
            vector cen[2] = {compute_cell_center(lhs, 4), compute_cell_center(rhs, 4)};
            double vol[2] = {compute_cell_volume(lhs, 4), compute_cell_volume(rhs, 4)};
            return weighted_average(cen, vol, 2);
        }
        case 6: {
            vector lhs[5] = {coord[0], coord[1], coord[4], coord[3], coord[5]};
            vector rhs[4] = {coord[0], coord[1], coord[2], coord[5]};
            vector cen[2] = {compute_cell_center(lhs, 5), compute_cell_center(rhs, 4)};
            double vol[2] = {compute_cell_volume(lhs, 5), compute_cell_volume(rhs, 4)};
            return weighted_average(cen, vol, 2);
        }
        case 8: {
            vector lhs[6] = {coord[0], coord[1], coord[2], coord[4], coord[5], coord[6]};
            vector rhs[6] = {coord[0], coord[2], coord[3], coord[4], coord[6], coord[7]};
            vector cen[2] = {compute_cell_center(lhs, 6), compute_cell_center(rhs, 6)};
            double vol[2] = {compute_cell_volume(lhs, 6), compute_cell_volume(rhs, 6)};
            return weighted_average(cen, vol, 2);
        }
        default: abort();
    }
}

/* Exchange centers for outer cells so neighbor/periodic copies receive an interior center. */
static void collect_centers(const MeshNeighbors *neighbors, vector *center)
{
    Arena save = arena_save();

    vector *send = arena_malloc(neighbors->send.off[neighbors->num], sizeof(*send));
    int tag = sync_tag();
    MPI_Request *req = arena_malloc(neighbors->num, sizeof(*req));
    for (long i = 0; i < neighbors->num; i++) {
        for (long j = neighbors->send.off[i]; j < neighbors->send.off[i + 1]; j++) {
            send[j] = center[neighbors->send.idx[j]];
        }
        int sendcount = neighbors->send.off[i + 1] - neighbors->send.off[i];
        int recvcount = neighbors->recv_off[i + 1] - neighbors->recv_off[i];
        MPI_Isendrecv(&send[neighbors->send.off[i]], sendcount, vector_type, neighbors->rank[i],
                      tag, &center[neighbors->recv_off[i]], recvcount, vector_type,
                      neighbors->rank[i], tag, sync.comm, &req[i]);
    }
    MPI_Waitall(neighbors->num, req, MPI_STATUSES_IGNORE);

    arena_load(save);
}

static void compute_cell_geometry(const MeshNodes *nodes, MeshCells *cells, const MeshFaces *faces,
                                  const MeshEntities *entities, const MeshNeighbors *neighbors)
{
    double *volume = arena_calloc(cells->num, sizeof(*volume));
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
        vector nrm = faces->basis[i].x;
        vector inc = vector_mul(vector_abs(nrm), faces->area[i] / 2);
        projection[left] = vector_add(projection[left], inc);
        if (right < cells->num_inner) {
            projection[right] = vector_add(projection[right], inc);
        }
        else if (i < faces->num_inner + faces->num_ghost) {  // mirror left center across the faces
            vector c2f = vector_sub(faces->center[i], center[left]);
            center[right] = vector_add(center[left], vector_mul(nrm, 2 * vector_dot(c2f, nrm)));
        }
    }

    collect_centers(neighbors, center);

    for (long i = entities->num_inner + entities->num_ghost; i < entities->num; i++) {
        for (long j = entities->cell_off[i]; j < entities->cell_off[i + 1]; j++) {
            center[j] = vector_sub(center[j], entities->offset[i]);
        }
    }

    cells->volume = volume;
    cells->center = center;
    cells->projection = projection;
}

static void compute_face_weights(const MeshCells *cells, MeshFaces *faces)
{
    Arena save = arena_save();

    double *r11 = arena_calloc(cells->num_inner, sizeof(*r11));
    double *r12 = arena_calloc(cells->num_inner, sizeof(*r12));
    double *r22 = arena_calloc(cells->num_inner, sizeof(*r22));
    double *r13 = arena_calloc(cells->num_inner, sizeof(*r13));
    double *r23 = arena_calloc(cells->num_inner, sizeof(*r23));
    double *r33 = arena_calloc(cells->num_inner, sizeof(*r33));
    for (long i = 0; i < faces->num; i++) {
        long left = faces->cell[i].left;
        long right = faces->cell[i].right;
        vector delta = vector_sub(cells->center[right], cells->center[left]);
        double theta2 = sq(1.0 / vector_len(delta));
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

    vector *weight = arena_calloc(faces->num, sizeof(*weight));
    for (long i = 0; i < faces->num; i++) {
        long left = faces->cell[i].left;
        long right = faces->cell[i].right;
        vector delta = vector_sub(cells->center[right], cells->center[left]);
        double theta2 = sq(1.0 / vector_len(delta));
        double beta = (r12[left] * r23[left] - r13[left] * r22[left]) / (r11[left] * r22[left]);
        vector alpha;
        alpha.x = delta.x / sq(r11[left]);
        alpha.y = (delta.y - r12[left] / r11[left] * delta.x) / sq(r22[left]);
        alpha.z = (delta.z - r23[left] / r22[left] * delta.y + beta * delta.x) / sq(r33[left]);
        if (!isclose(cells->center[left].x, cells->center[right].x)) {
            weight[i].x += alpha.x;
        }
        if (!isclose(cells->center[left].y, cells->center[right].y)) {
            weight[i].x += -r12[left] / r11[left] * alpha.y;
            weight[i].y += alpha.y;
        }
        if (!isclose(cells->center[left].z, cells->center[right].z)) {
            weight[i].x += beta * alpha.z;
            weight[i].y += -r23[left] / r22[left] * alpha.z;
            weight[i].z += alpha.z;
        }
        weight[i] = vector_mul(weight[i], theta2);
        assert(isfinite(weight[i].x));
        assert(isfinite(weight[i].y));
        assert(isfinite(weight[i].z));
    }

    arena_load(save);

    faces->weight = arena_smuggle(weight, faces->num, sizeof(*weight));
}

void mesh_build(Mesh *mesh)
{
    connect_cells(&mesh->nodes, &mesh->cells);

    create_faces(&mesh->cells, &mesh->faces);
    reorder_faces(&mesh->cells, &mesh->faces);

    compute_face_entities(&mesh->faces, &mesh->entities);

    compute_send_graph(&mesh->nodes, &mesh->cells, &mesh->entities, &mesh->neighbors);

    compute_face_geometry(&mesh->nodes, &mesh->cells, &mesh->faces);
    compute_cell_geometry(&mesh->nodes, &mesh->cells, &mesh->faces, &mesh->entities,
                          &mesh->neighbors);

    compute_face_weights(&mesh->cells, &mesh->faces);
}
