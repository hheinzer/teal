#include <assert.h>
#include <limits.h>
#include <math.h>
#include <metis.h>
#include <stdlib.h>

#include "private.h"
#include "sync2.h"
#include "teal2.h"
#include "utils2.h"

static void connect_cells(Mesh2 *mesh)
{
    idx_t ne = mesh->cells.num;  // NOLINT(readability-identifier-length)
    idx_t nn = mesh->nodes.num;  // NOLINT(readability-identifier-length)

    idx_t *eptr = teal2_calloc(mesh->cells.num + 1, sizeof(*eptr));
    idx_t *eind = teal2_calloc(mesh->cells.num * (MAX_CELL_NODES + 1), sizeof(*eind));

    for (int i = 0; i < mesh->cells.num; i++) {
        eptr[i + 1] = eptr[i];
        for (int j = mesh->cells.node.off[i]; j < mesh->cells.node.off[i + 1]; j++) {
            assert(eptr[i + 1] - eptr[i] < MAX_CELL_NODES + 1);
            eind[eptr[i + 1]++] = mesh->cells.node.idx[j];
        }
        if (i >= mesh->cells.num_inner) {
            assert(eptr[i + 1] - eptr[i] < MAX_CELL_NODES + 1);
            eind[eptr[i + 1]++] = nn++;  // add dummy node to enforce ncommon=3
        }
    }

    idx_t ncommon = 3;
    idx_t numflag = 0;
    idx_t *xadj;
    idx_t *adjncy;
    int ret = METIS_MeshToDual(&ne, &nn, eptr, eind, &ncommon, &numflag, &xadj, &adjncy);
    if (ret != METIS_OK) {
        teal2_error("METIS_MeshToDual failed (%d)", ret);
    }

    int *cell_off = teal2_calloc(mesh->cells.num + 1, sizeof(*cell_off));
    int *cell_idx = teal2_calloc(mesh->cells.num * MAX_CELL_FACES, sizeof(*cell_idx));

    for (int i = 0; i < mesh->cells.num; i++) {
        cell_off[i + 1] = cell_off[i];
        for (idx_t j = xadj[i]; j < xadj[i + 1]; j++) {
            if (i < mesh->cells.num_inner || adjncy[j] < mesh->cells.num_inner) {
                assert(cell_off[i + 1] - cell_off[i] < MAX_CELL_FACES);
                assert(adjncy[j] <= INT_MAX);
                cell_idx[cell_off[i + 1]++] = (int)adjncy[j];
            }
        }
    }

    mesh->cells.cell.off = cell_off;
    mesh->cells.cell.idx = teal2_realloc(cell_idx, cell_off[mesh->cells.num], sizeof(*cell_idx));

    teal2_free(eptr);
    teal2_free(eind);
    METIS_Free(xadj);
    METIS_Free(adjncy);
}

static void reorder_cells(Mesh2 *mesh)
{
    int *map = teal2_calloc(mesh->cells.num_inner, sizeof(*map));
    int *degree = teal2_calloc(mesh->cells.num_inner, sizeof(*degree));

    for (int i = 0; i < mesh->cells.num_inner; i++) {
        map[i] = -1;
        degree[i] = mesh->cells.cell.off[i + 1] - mesh->cells.cell.off[i];
    }

    int *queue = teal2_calloc(mesh->cells.num_inner, sizeof(*queue));
    int beg = 0;
    int end = 0;

    assert(mesh->cells.num_inner > 0);
    int seed = rand() % mesh->cells.num_inner;

    int num = 0;
    map[seed] = num++;
    queue[end++] = seed;
    while (beg < end) {
        int cur = queue[beg++];
        int cell[MAX_CELL_FACES];
        int num_cells = 0;
        for (int i = mesh->cells.cell.off[cur]; i < mesh->cells.cell.off[cur + 1]; i++) {
            int idx = mesh->cells.cell.idx[i];
            if (idx < mesh->cells.num_inner && map[idx] == -1) {
                assert(num_cells < MAX_CELL_FACES);
                cell[num_cells++] = idx;
                map[idx] = -2;
            }
        }
        for (int i = 1; i < num_cells; i++) {
            int idx = cell[i];
            int deg = degree[idx];
            int j = i - 1;  // NOLINT(readability-identifier-length)
            while (j >= 0 && degree[cell[j]] < deg) {
                cell[j + 1] = cell[j];
                j -= 1;
            }
            cell[j + 1] = idx;
        }
        for (int i = 0; i < num_cells; i++) {
            map[cell[i]] = num++;
            queue[end++] = cell[i];
        }
    }
    assert(num == mesh->cells.num_inner);

    mesh2_reorder_cells(mesh, map, 0, mesh->cells.num_inner);

    teal2_free(map);
    teal2_free(degree);
    teal2_free(queue);
}

static void collect_globals(Mesh2 *mesh, const int *map)
{
    long prefix = mesh->nodes.num_inner;
    sync2_prefix(&prefix, 1, MPI_LONG);

    for (int i = 0; i < mesh->nodes.num_inner; i++) {
        mesh->nodes.global[i] = prefix + i;
    }

    long *off_nodes = teal2_calloc(sync2.size + 1, sizeof(*off_nodes));
    sync2_offsets(&(long){mesh->nodes.num_inner}, off_nodes, 1, MPI_LONG);

    int num_nodes = mesh->nodes.num - mesh->nodes.num_inner;
    int *rank = teal2_calloc(num_nodes, sizeof(*rank));

    int num = 0;
    for (int i = mesh->nodes.num_inner; i < mesh->nodes.num; i++) {
        rank[num] = digitize(&mesh->nodes.global[i], off_nodes, sync2.size, sizeof(*off_nodes),
                             compare_long);
        assert(0 <= rank[num] && rank[num] < sync2.size);
        num += 1;
    }
    assert(num == num_nodes);

    int *num_recv = teal2_calloc(sync2.size, sizeof(*num_recv));
    for (int i = 0; i < num_nodes; i++) {
        num_recv[rank[i]] += 1;
    }

    int *off_recv = teal2_calloc(sync2.size + 1, sizeof(*off_recv));
    for (int i = 0; i < sync2.size; i++) {
        off_recv[i + 1] = off_recv[i] + num_recv[i];
    }

    int *cur_recv = teal2_calloc(sync2.size, sizeof(*cur_recv));
    copy(cur_recv, off_recv, sync2.size, sizeof(*cur_recv));

    int *idx_recv = teal2_calloc(num_nodes, sizeof(*idx_recv));
    num = 0;
    for (int i = mesh->nodes.num_inner; i < mesh->nodes.num; i++) {
        long idx_local = mesh->nodes.global[i] - off_nodes[rank[num]];
        assert(idx_local <= INT_MAX);
        idx_recv[cur_recv[rank[num]]++] = (int)idx_local;
        num += 1;
    }
    for (int i = 0; i < sync2.size; i++) {
        assert(cur_recv[i] == off_recv[i + 1]);
    }
    assert(num == num_nodes);

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

    long *send = teal2_calloc(tot_send, sizeof(*send));
    for (int i = 0; i < tot_send; i++) {
        send[i] = mesh->nodes.global[map[idx_send[i]]];
    }

    long *recv = teal2_calloc(num_nodes, sizeof(*recv));
    MPI_Alltoallv(send, num_send, off_send, MPI_LONG, recv, num_recv, off_recv, MPI_LONG,
                  sync2.comm);

    copy(cur_recv, off_recv, sync2.size, sizeof(*cur_recv));

    num = 0;
    for (int i = mesh->nodes.num_inner; i < mesh->nodes.num; i++) {
        mesh->nodes.global[i] = recv[cur_recv[rank[num]]++];
        num += 1;
    }
    for (int i = 0; i < sync2.size; i++) {
        assert(cur_recv[i] == off_recv[i + 1]);
    }
    assert(num == num_nodes);

    teal2_free(off_nodes);
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
}

static void reorder_nodes(Mesh2 *mesh)
{
    int *map = teal2_calloc(mesh->nodes.num, sizeof(*map));
    for (int i = 0; i < mesh->nodes.num; i++) {
        map[i] = -1;
    }

    int num = 0;
    for (int i = 0; i < mesh->cells.num; i++) {
        for (int j = mesh->cells.node.off[i]; j < mesh->cells.node.off[i + 1]; j++) {
            int idx = mesh->cells.node.idx[j];
            if (idx < mesh->nodes.num_inner && map[idx] == -1) {
                map[idx] = num++;
            }
        }
    }
    assert(num == mesh->nodes.num_inner);

    for (int i = 0; i < mesh->cells.num; i++) {
        for (int j = mesh->cells.node.off[i]; j < mesh->cells.node.off[i + 1]; j++) {
            int idx = mesh->cells.node.idx[j];
            if (idx >= mesh->nodes.num_inner && map[idx] == -1) {
                map[idx] = num++;
            }
        }
    }
    assert(num == mesh->nodes.num);

    mesh2_reorder_nodes(mesh, map, 0, mesh->nodes.num);

    collect_globals(mesh, map);

    teal2_free(map);
}

typedef struct {
    int inner;
    int left;
    int right;
    int num_nodes;
    int node[MAX_FACE_NODES];
} Face;

static int compare_face(const void *lhs_, const void *rhs_)
{
    const Face *lhs = lhs_;
    const Face *rhs = rhs_;
    if (lhs->inner != rhs->inner) {
        return (lhs->inner > rhs->inner) ? -1 : +1;
    }
    if (lhs->inner) {
        return (lhs->left > rhs->left) - (lhs->left < rhs->left);
    }
    return (lhs->right > rhs->right) - (lhs->right < rhs->right);
}

static void compute_face_nodes(Face *face, const Mesh2 *mesh, int left, int right)
{
    face->num_nodes = 0;
    for (int i = mesh->cells.node.off[left]; i < mesh->cells.node.off[left + 1]; i++) {
        for (int j = mesh->cells.node.off[right]; j < mesh->cells.node.off[right + 1]; j++) {
            if (mesh->cells.node.idx[i] == mesh->cells.node.idx[j]) {
                assert(face->num_nodes < MAX_FACE_NODES);
                face->node[face->num_nodes++] = mesh->cells.node.idx[i];
            }
        }
    }
}

static void correct_face_node_order(Face *face, const Mesh2 *mesh)
{
    switch (face->num_nodes) {
        case 3: return;
        case 4:
            for (int i = 0; i < 3; i++) {
                Vector coord[4];
                for (int j = 0; j < 4; j++) {
                    coord[j] = mesh->nodes.coord[face->node[j]];
                }
                Vector a2b = vector2_sub(coord[1], coord[0]);
                Vector a2c = vector2_sub(coord[2], coord[0]);
                Vector a2d = vector2_sub(coord[3], coord[0]);
                if (vector2_dot(vector2_cross(a2b, a2c), vector2_cross(a2c, a2d)) > 0) {
                    return;
                }
                if (i + 2 < face->num_nodes) {
                    int swap = face->node[i + 1];
                    face->node[i + 1] = face->node[i + 2];
                    face->node[i + 2] = swap;
                }
            }
            teal2_error("quad vertex order could not be made valid");
        default: teal2_error("invalid number of nodes (%d)", face->num_nodes);
    }
}

static void create_faces(Mesh2 *mesh)
{
    int max_faces = mesh->cells.num_inner * MAX_CELL_FACES;
    Face *face = teal2_calloc(max_faces, sizeof(*face));

    int num_faces = 0;
    for (int i = 0; i < mesh->cells.num_inner; i++) {
        for (int j = mesh->cells.cell.off[i]; j < mesh->cells.cell.off[i + 1]; j++) {
            int left = i;
            int right = mesh->cells.cell.idx[j];
            if (left < right) {
                assert(num_faces < max_faces);
                face[num_faces].inner = (right < mesh->cells.num_inner);
                face[num_faces].left = left;
                face[num_faces].right = right;
                compute_face_nodes(&face[num_faces], mesh, left, right);
                correct_face_node_order(&face[num_faces], mesh);
                num_faces += 1;
            }
        }
    }
    sort(face, num_faces, sizeof(*face), compare_face);

    int *node_off = teal2_calloc(num_faces + 1, sizeof(*node_off));
    int *node_idx = teal2_calloc(num_faces * MAX_FACE_NODES, sizeof(*node_idx));
    IntPair *cell_idx = teal2_calloc(num_faces, sizeof(*cell_idx));

    int num_inner = 0;
    int num_boundary = 0;
    for (int i = 0; i < num_faces; i++) {
        assert(face[i].left < mesh->cells.num_inner);
        if (face[i].right < mesh->cells.num_inner) {
            num_inner += 1;
        }
        else if (face[i].right < mesh->cells.off_boundary) {
            num_boundary += 1;
        }
        node_off[i + 1] = node_off[i];
        for (int j = 0; j < face[i].num_nodes; j++) {
            node_idx[node_off[i + 1]++] = face[i].node[j];
        }
        cell_idx[i].left = face[i].left;
        cell_idx[i].right = face[i].right;
    }
    assert(num_boundary == mesh->cells.off_boundary - mesh->cells.num_inner);

    mesh->faces.num = num_faces;
    mesh->faces.num_inner = num_inner;
    mesh->faces.off_boundary = num_inner + num_boundary;
    mesh->faces.node.off = node_off;
    mesh->faces.node.idx = node_idx;
    mesh->faces.cell_idx = cell_idx;

    teal2_free(face);
}

static void create_face_entities(Mesh2 *mesh)
{
    int *face_off = teal2_calloc(mesh->entities.num + 1, sizeof(*face_off));
    face_off[0] = mesh->faces.num_inner;
    for (int i = 0; i < mesh->entities.num; i++) {
        face_off[i + 1] = face_off[i];
        if (i >= mesh->entities.num_inner) {
            face_off[i + 1] += mesh->entities.cell_off[i + 1] - mesh->entities.cell_off[i];
        }
    }
    mesh->entities.face_off = face_off;
}

static double compute_face_area(const Vector *coord, int num_nodes)
{
    switch (num_nodes) {
        case 3: {
            Vector a2b = vector2_sub(coord[1], coord[0]);
            Vector a2c = vector2_sub(coord[2], coord[0]);
            return vector2_norm(vector2_cross(a2b, a2c)) / 2;
        }
        case 4: {
            Vector lhs[3] = {coord[0], coord[1], coord[2]};
            Vector rhs[3] = {coord[0], coord[2], coord[3]};
            return compute_face_area(lhs, 3) + compute_face_area(rhs, 3);
        }
        default: teal2_error("invalid number of nodes (%d)", num_nodes);
    }
}

static void compute_face_areas(Mesh2 *mesh)
{
    double *area = teal2_calloc(mesh->faces.num, sizeof(*area));
    for (int i = 0; i < mesh->faces.num; i++) {
        Vector coord[MAX_FACE_NODES];
        int num_nodes = 0;
        for (int j = mesh->faces.node.off[i]; j < mesh->faces.node.off[i + 1]; j++) {
            coord[num_nodes++] = mesh->nodes.coord[mesh->faces.node.idx[j]];
        }
        area[i] = compute_face_area(coord, num_nodes);
    }
    mesh->faces.area = area;
}

static Vector compute_face_center(const Vector *coord, int num_nodes)
{
    switch (num_nodes) {
        case 3: {
            Vector sum = vector2_add(coord[0], vector2_add(coord[1], coord[2]));
            return vector2_div(sum, 3);
        }
        case 4: {
            Vector lhs[3] = {coord[0], coord[1], coord[2]};
            Vector rhs[3] = {coord[0], coord[2], coord[3]};
            double area[2] = {compute_face_area(lhs, 3), compute_face_area(rhs, 3)};
            Vector center[2] = {compute_face_center(lhs, 3), compute_face_center(rhs, 3)};
            Vector mul[2] = {vector2_mul(area[0], center[0]), vector2_mul(area[1], center[1])};
            return vector2_div(vector2_add(mul[0], mul[1]), area[0] + area[1]);
        }
        default: teal2_error("invalid number of nodes (%d)", num_nodes);
    }
}

static void compute_face_centers(Mesh2 *mesh)
{
    Vector *center = teal2_calloc(mesh->faces.num, sizeof(*center));
    for (int i = 0; i < mesh->faces.num; i++) {
        Vector coord[MAX_FACE_NODES];
        int num_nodes = 0;
        for (int j = mesh->faces.node.off[i]; j < mesh->faces.node.off[i + 1]; j++) {
            coord[num_nodes++] = mesh->nodes.coord[mesh->faces.node.idx[j]];
        }
        center[i] = compute_face_center(coord, num_nodes);
    }
    mesh->faces.center = center;
}

static Vector compute_face_normal(const Vector *coord, int num_nodes)
{
    switch (num_nodes) {
        case 3: {
            Vector a2b = vector2_sub(coord[1], coord[0]);
            Vector a2c = vector2_sub(coord[2], coord[0]);
            Vector normal = vector2_cross(a2b, a2c);
            return vector2_div(normal, vector2_norm(normal));
        }
        case 4: {
            Vector normal = {0};
            for (int i = 0; i < 4; i++) {
                vector2_iadd(&normal, vector2_cross(coord[i], coord[(i + 1) % 4]));
            }
            return vector2_div(normal, vector2_norm(normal));
        }
        default: teal2_error("invalid number of nodes (%d)", num_nodes);
    }
}

static Matrix compute_face_basis(const Vector *coord, int num_nodes)
{
    Vector normal = compute_face_normal(coord, num_nodes);
    Matrix basis = {.x = normal};
    double nqz = hypot(normal.x, normal.y);
    double nqy = hypot(normal.x, normal.z);
    if (nqz > nqy) {
        basis.y = (Vector){-normal.y / nqz, normal.x / nqz, 0};
        basis.z = (Vector){-normal.x * normal.z / nqz, -normal.y * normal.z / nqz, nqz};
    }
    else {
        basis.y = (Vector){-normal.x * normal.y / nqy, nqy, -normal.y * normal.z / nqy};
        basis.z = (Vector){-normal.z / nqy, 0, normal.x / nqy};
    }
    return basis;
}

static void correct_face_basis(Matrix *basis, const Mesh2 *mesh, int idx)
{
    int left = mesh->faces.cell_idx[idx].left;

    Vector mean = {0};
    for (int i = mesh->cells.node.off[left]; i < mesh->cells.node.off[left + 1]; i++) {
        vector2_iadd(&mean, mesh->nodes.coord[mesh->cells.node.idx[i]]);
    }
    vector2_idiv(&mean, mesh->cells.node.off[left + 1] - mesh->cells.node.off[left]);

    if (vector2_dot(vector2_sub(mean, mesh->faces.center[idx]), basis->x) > 0) {
        vector2_imul(&basis->x, -1);
        vector2_imul(&basis->y, -1);
    }
}

static void compute_face_bases(Mesh2 *mesh)
{
    Matrix *basis = teal2_calloc(mesh->faces.num, sizeof(*basis));
    for (int i = 0; i < mesh->faces.num; i++) {
        Vector coord[MAX_FACE_NODES];
        int num_nodes = 0;
        for (int j = mesh->faces.node.off[i]; j < mesh->faces.node.off[i + 1]; j++) {
            coord[num_nodes++] = mesh->nodes.coord[mesh->faces.node.idx[j]];
        }
        basis[i] = compute_face_basis(coord, num_nodes);
        correct_face_basis(&basis[i], mesh, i);
    }
    mesh->faces.basis = basis;
}

static double compute_cell_volume(const Vector *coord, int num_nodes)
{
    switch (num_nodes) {
        case 4: {
            Vector a2b = vector2_sub(coord[1], coord[0]);
            Vector a2c = vector2_sub(coord[2], coord[0]);
            Vector a2d = vector2_sub(coord[3], coord[0]);
            return fabs(vector2_dot(a2b, vector2_cross(a2c, a2d))) / 6;
        }
        case 5: {
            Vector lhs[4] = {coord[0], coord[1], coord[2], coord[4]};
            Vector rhs[4] = {coord[0], coord[2], coord[3], coord[4]};
            return compute_cell_volume(lhs, 4) + compute_cell_volume(rhs, 4);
        }
        case 6: {
            Vector lhs[5] = {coord[0], coord[1], coord[4], coord[3], coord[5]};
            Vector rhs[4] = {coord[0], coord[1], coord[2], coord[5]};
            return compute_cell_volume(lhs, 5) + compute_cell_volume(rhs, 4);
        }
        case 8: {
            Vector lhs[6] = {coord[0], coord[1], coord[2], coord[4], coord[5], coord[6]};
            Vector rhs[6] = {coord[0], coord[2], coord[3], coord[4], coord[6], coord[7]};
            return compute_cell_volume(lhs, 6) + compute_cell_volume(rhs, 6);
        }
        default: teal2_error("invalid number of nodes (%d)", num_nodes);
    }
}

static void compute_cell_volumes(Mesh2 *mesh)
{
    double *volume = teal2_calloc(mesh->cells.num, sizeof(*volume));
    for (int i = 0; i < mesh->cells.num_inner; i++) {
        Vector coord[MAX_CELL_NODES];
        int num_nodes = 0;
        for (int j = mesh->cells.node.off[i]; j < mesh->cells.node.off[i + 1]; j++) {
            coord[num_nodes++] = mesh->nodes.coord[mesh->cells.node.idx[j]];
        }
        volume[i] = compute_cell_volume(coord, num_nodes);
    }
    mesh->cells.volume = volume;
}

static Vector compute_cell_center(const Vector *coord, int num_nodes)
{
    switch (num_nodes) {
        case 4: {
            Vector sum[2] = {vector2_add(coord[0], coord[1]), vector2_add(coord[2], coord[3])};
            return vector2_div(vector2_add(sum[0], sum[1]), 4);
        }
        case 5: {
            Vector lhs[4] = {coord[0], coord[1], coord[2], coord[4]};
            Vector rhs[4] = {coord[0], coord[2], coord[3], coord[4]};
            double volume[2] = {compute_cell_volume(lhs, 4), compute_cell_volume(rhs, 4)};
            Vector center[2] = {compute_cell_center(lhs, 4), compute_cell_center(rhs, 4)};
            Vector mul[2] = {vector2_mul(volume[0], center[0]), vector2_mul(volume[1], center[1])};
            return vector2_div(vector2_add(mul[0], mul[1]), volume[0] + volume[1]);
        }
        case 6: {
            Vector lhs[5] = {coord[0], coord[1], coord[4], coord[3], coord[5]};
            Vector rhs[4] = {coord[0], coord[1], coord[2], coord[5]};
            double volume[2] = {compute_cell_volume(lhs, 5), compute_cell_volume(rhs, 4)};
            Vector center[2] = {compute_cell_center(lhs, 5), compute_cell_center(rhs, 4)};
            Vector mul[2] = {vector2_mul(volume[0], center[0]), vector2_mul(volume[1], center[1])};
            return vector2_div(vector2_add(mul[0], mul[1]), volume[0] + volume[1]);
        }
        case 8: {
            Vector lhs[6] = {coord[0], coord[1], coord[2], coord[4], coord[5], coord[6]};
            Vector rhs[6] = {coord[0], coord[2], coord[3], coord[4], coord[6], coord[7]};
            double volume[2] = {compute_cell_volume(lhs, 6), compute_cell_volume(rhs, 6)};
            Vector center[2] = {compute_cell_center(lhs, 6), compute_cell_center(rhs, 6)};
            Vector mul[2] = {vector2_mul(volume[0], center[0]), vector2_mul(volume[1], center[1])};
            return vector2_div(vector2_add(mul[0], mul[1]), volume[0] + volume[1]);
        }
        default: teal2_error("invalid number of nodes (%d)", num_nodes);
    }
}

static void collect_centers(Vector *center, const Mesh2 *mesh)
{
    int tot_send = mesh->neighbors.send.off[mesh->neighbors.num];
    Vector *send = teal2_calloc(tot_send, sizeof(*send));

    MPI_Datatype type;
    MPI_Type_contiguous(3, MPI_DOUBLE, &type);
    MPI_Type_commit(&type);

    MPI_Request *req_recv = teal2_calloc(mesh->neighbors.num, sizeof(*req_recv));
    MPI_Request *req_send = teal2_calloc(mesh->neighbors.num, sizeof(*req_send));
    for (int i = 0; i < mesh->neighbors.num; i++) {
        for (int j = mesh->neighbors.send.off[i]; j < mesh->neighbors.send.off[i + 1]; j++) {
            send[j] = center[mesh->neighbors.send.idx[j]];
        }
        int num_recv = mesh->neighbors.recv_off[i + 1] - mesh->neighbors.recv_off[i];
        int num_send = mesh->neighbors.send.off[i + 1] - mesh->neighbors.send.off[i];
        MPI_Irecv(&center[mesh->neighbors.recv_off[i]], num_recv, type, mesh->neighbors.rank[i],
                  mesh->neighbors.tag[i][0], sync2.comm, &req_recv[i]);
        MPI_Isend(&send[mesh->neighbors.send.off[i]], num_send, type, mesh->neighbors.rank[i],
                  mesh->neighbors.tag[i][1], sync2.comm, &req_send[i]);
    }
    MPI_Waitall(mesh->neighbors.num, req_recv, MPI_STATUSES_IGNORE);
    MPI_Waitall(mesh->neighbors.num, req_send, MPI_STATUSES_IGNORE);

    MPI_Type_free(&type);

    teal2_free(send);
    teal2_free(req_recv);
    teal2_free(req_send);
}

static void compute_cell_centers(Mesh2 *mesh)
{
    Vector *center = teal2_calloc(mesh->cells.num, sizeof(*center));

    for (int i = 0; i < mesh->cells.num_inner; i++) {
        Vector coord[MAX_CELL_NODES];
        int num_nodes = 0;
        for (int j = mesh->cells.node.off[i]; j < mesh->cells.node.off[i + 1]; j++) {
            coord[num_nodes++] = mesh->nodes.coord[mesh->cells.node.idx[j]];
        }
        center[i] = compute_cell_center(coord, num_nodes);
    }

    for (int i = mesh->faces.num_inner; i < mesh->faces.off_boundary; i++) {
        int left = mesh->faces.cell_idx[i].left;
        int right = mesh->faces.cell_idx[i].right;
        Vector normal = mesh->faces.basis[i].x;
        Vector offset = vector2_sub(mesh->faces.center[i], center[left]);
        double factor = 2 * vector2_dot(normal, offset);
        center[right] = vector2_add(center[left], vector2_mul(factor, normal));
    }

    collect_centers(center, mesh);

    for (int i = mesh->entities.off_boundary; i < mesh->entities.num; i++) {
        for (int j = mesh->entities.cell_off[i]; j < mesh->entities.cell_off[i + 1]; j++) {
            center[j] = matrix_vector(mesh->entities.rotation[i], center[j]);
            vector2_iadd(&center[j], mesh->entities.translation[i]);
        }
    }

    mesh->cells.center = center;
}

static void compute_cell_projections(Mesh2 *mesh)
{
    Vector *projection = teal2_calloc(mesh->cells.num, sizeof(*projection));
    for (int i = 0; i < mesh->faces.num; i++) {
        int left = mesh->faces.cell_idx[i].left;
        int right = mesh->faces.cell_idx[i].right;
        Vector normal = mesh->faces.basis[i].x;
        Vector increment = vector2_mul(mesh->faces.area[i] / 2, vector2_abs(normal));
        vector2_iadd(&projection[left], increment);
        if (right < mesh->cells.num_inner) {
            vector2_iadd(&projection[right], increment);
        }
    }
    mesh->cells.projection = projection;
}

typedef struct {
    int left, right, idx;
} Map;

static int compare_map(const void *lhs_, const void *rhs_)
{
    const Map *lhs = lhs_;
    const Map *rhs = rhs_;
    if (lhs->left != rhs->left) {
        return (lhs->left < rhs->left) ? -1 : +1;
    }
    return (lhs->right > rhs->right) - (lhs->right < rhs->right);
}

static void compute_cell_offsets(Mesh2 *mesh)
{
    Map *map = teal2_calloc(mesh->faces.num, sizeof(*map));
    for (int i = 0; i < mesh->faces.num; i++) {
        map[i].left = mesh->faces.cell_idx[i].left;
        map[i].right = mesh->faces.cell_idx[i].right;
        map[i].idx = i;
    }
    sort(map, mesh->faces.num, sizeof(*map), compare_map);

    int num_offsets = mesh->cells.cell.off[mesh->cells.num_inner];
    Vector *offset = teal2_calloc(num_offsets, sizeof(*offset));
    for (int i = 0; i < mesh->cells.num_inner; i++) {
        for (int j = mesh->cells.cell.off[i]; j < mesh->cells.cell.off[i + 1]; j++) {
            int idx = mesh->cells.cell.idx[j];
            Map key = {.left = (i < idx) ? i : idx, .right = (i > idx) ? i : idx};
            Map *val = search(&key, map, mesh->faces.num, sizeof(*map), compare_map);
            assert(val);
            offset[j] = vector2_sub(mesh->faces.center[val->idx], mesh->cells.center[i]);
        }
    }
    mesh->cells.offset = offset;

    teal2_free(map);
}

static Vector compute_face_weight(Vector delta, const double *r11, const double *r12,
                                  const double *r22, const double *r13, const double *r23,
                                  const double *r33, int idx)
{
    Vector weight = {0};
    double beta = (r12[idx] * r23[idx] - r13[idx] * r22[idx]) / (r11[idx] * r22[idx]);
    Vector alpha;
    alpha.x = delta.x / sq(r11[idx]);
    alpha.y = (delta.y - r12[idx] / r11[idx] * delta.x) / sq(r22[idx]);
    alpha.z = (delta.z - r23[idx] / r22[idx] * delta.y + beta * delta.x) / sq(r33[idx]);
    if (!isclose(delta.x, 0)) {
        weight.x += alpha.x;
    }
    if (!isclose(delta.y, 0)) {
        weight.x += -r12[idx] / r11[idx] * alpha.y;
        weight.y += alpha.y;
    }
    if (!isclose(delta.z, 0)) {
        weight.x += beta * alpha.z;
        weight.y += -r23[idx] / r22[idx] * alpha.z;
        weight.z += alpha.z;
    }
    double theta2 = 1 / vector2_norm2(delta);
    vector2_imul(&weight, theta2);
    assert(isfinite(weight.x) && isfinite(weight.y) && isfinite(weight.z));
    return weight;
}

static void compute_face_weights(Mesh2 *mesh)
{
    double *r11 = teal2_calloc(mesh->cells.num_inner, sizeof(*r11));
    double *r12 = teal2_calloc(mesh->cells.num_inner, sizeof(*r12));
    double *r22 = teal2_calloc(mesh->cells.num_inner, sizeof(*r22));
    double *r13 = teal2_calloc(mesh->cells.num_inner, sizeof(*r13));
    double *r23 = teal2_calloc(mesh->cells.num_inner, sizeof(*r23));
    double *r33 = teal2_calloc(mesh->cells.num_inner, sizeof(*r33));
    for (int i = 0; i < mesh->faces.num; i++) {
        int left = mesh->faces.cell_idx[i].left;
        int right = mesh->faces.cell_idx[i].right;
        Vector delta = vector2_sub(mesh->cells.center[right], mesh->cells.center[left]);
        double theta2 = 1 / vector2_norm2(delta);
        r11[left] += theta2 * delta.x * delta.x;
        r12[left] += theta2 * delta.x * delta.y;
        r22[left] += theta2 * delta.y * delta.y;
        r13[left] += theta2 * delta.x * delta.z;
        r23[left] += theta2 * delta.y * delta.z;
        r33[left] += theta2 * delta.z * delta.z;
        if (right < mesh->cells.num_inner) {
            r11[right] += theta2 * delta.x * delta.x;
            r12[right] += theta2 * delta.x * delta.y;
            r22[right] += theta2 * delta.y * delta.y;
            r13[right] += theta2 * delta.x * delta.z;
            r23[right] += theta2 * delta.y * delta.z;
            r33[right] += theta2 * delta.z * delta.z;
        }
    }
    for (int i = 0; i < mesh->cells.num_inner; i++) {
        r11[i] = sqrt(r11[i]);
        r12[i] = r12[i] / r11[i];
        r22[i] = sqrt(r22[i] - sq(r12[i]));
        r13[i] = r13[i] / r11[i];
        r23[i] = (r23[i] - r12[i] * r13[i]) / r22[i];
        r33[i] = sqrt(r33[i] - (sq(r13[i]) + sq(r23[i])));
    }

    VectorPair *weight = teal2_calloc(mesh->faces.num, sizeof(*weight));
    for (int i = 0; i < mesh->faces.num; i++) {
        int left = mesh->faces.cell_idx[i].left;
        int right = mesh->faces.cell_idx[i].right;
        Vector delta = vector2_sub(mesh->cells.center[right], mesh->cells.center[left]);
        weight[i].left = compute_face_weight(delta, r11, r12, r22, r13, r23, r33, left);
        if (right < mesh->cells.num_inner) {
            vector2_imul(&delta, -1);
            weight[i].right = compute_face_weight(delta, r11, r12, r22, r13, r23, r33, right);
        }
    }
    mesh->faces.weight = weight;

    teal2_free(r11);
    teal2_free(r12);
    teal2_free(r22);
    teal2_free(r13);
    teal2_free(r23);
    teal2_free(r33);
}

static void compute_face_offsets(Mesh2 *mesh)
{
    VectorPair *offset = teal2_calloc(mesh->faces.num, sizeof(*offset));
    for (int i = 0; i < mesh->faces.num; i++) {
        int left = mesh->faces.cell_idx[i].left;
        int right = mesh->faces.cell_idx[i].right;
        offset[i].left = vector2_sub(mesh->faces.center[i], mesh->cells.center[left]);
        offset[i].right = vector2_sub(mesh->faces.center[i], mesh->cells.center[right]);
    }
    mesh->faces.offset = offset;
}

static void compute_face_corrections(Mesh2 *mesh)
{
    UnitNorm *correction = teal2_calloc(mesh->faces.num, sizeof(*correction));
    for (int i = 0; i < mesh->faces.num; i++) {
        int left = mesh->faces.cell_idx[i].left;
        int right = mesh->faces.cell_idx[i].right;
        Vector delta = vector2_sub(mesh->cells.center[right], mesh->cells.center[left]);
        correction[i].norm = vector2_norm(delta);
        correction[i].unit = vector2_div(delta, correction[i].norm);
    }
    mesh->faces.correction = correction;
}

void mesh2_generate(Mesh2 *mesh)
{
    assert(mesh && !mesh->generated);

    connect_cells(mesh);
    reorder_cells(mesh);
    reorder_nodes(mesh);

    create_faces(mesh);
    create_face_entities(mesh);

    compute_face_areas(mesh);
    compute_face_centers(mesh);
    compute_face_bases(mesh);

    compute_cell_volumes(mesh);
    compute_cell_centers(mesh);
    compute_cell_projections(mesh);
    compute_cell_offsets(mesh);

    compute_face_weights(mesh);
    compute_face_offsets(mesh);
    compute_face_corrections(mesh);

    long tot_nodes = mesh->nodes.num;
    sync2_sum(&tot_nodes, 1, MPI_LONG);

    long tot_cells = mesh->cells.num;
    sync2_sum(&tot_cells, 1, MPI_LONG);

    long tot_indices = mesh->cells.node.off[mesh->cells.num];
    sync2_sum(&tot_indices, 1, MPI_LONG);

    if (tot_nodes > INT_MAX || tot_cells > INT_MAX || tot_indices > INT_MAX) {
        teal2.partitioned = 1;
    }

    mesh->generated = 1;
}
