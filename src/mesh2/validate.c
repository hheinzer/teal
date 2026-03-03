#include <assert.h>
#include <string.h>

#include "kdtree2.h"
#include "mesh2.h"
#include "sync2.h"
#include "teal2.h"
#include "utils2.h"

#define ensure(expr) \
    ((expr) ? (void)0 : teal2_error("`%s` failed (%s:%d: %s)", #expr, __FILE__, __LINE__, __func__))

static void validate_unique(const Vector *point, int num)
{
    Kdtree *tree = kdtree2_init(point, num);
    for (int i = 0; i < num; i++) {
        int idx[2];
        int len = kdtree2_nearest(tree, point[i], idx, 2);
        ensure(len > 0);
        for (int j = 0; j < len; j++) {
            ensure(0 <= idx[j] && idx[j] < num);
            if (idx[j] == i) {
                continue;
            }
            Vector lhs = point[i];
            Vector rhs = point[idx[j]];
            ensure(!(isclose(lhs.x, rhs.x) && isclose(lhs.y, rhs.y) && isclose(lhs.z, rhs.z)));
        }
    }
    kdtree2_deinit(tree);
}

static void validate_nodes(const Mesh *mesh)
{
    ensure(mesh->nodes.num > 0);
    ensure(mesh->nodes.num_inner >= 0);
    ensure(mesh->nodes.num_inner <= mesh->nodes.num);

    long beg_global = mesh->nodes.num_inner;
    sync2_prefix(&beg_global, 1, MPI_LONG);
    long end_global = beg_global + mesh->nodes.num_inner;

    long *global = teal2_calloc(mesh->nodes.num, sizeof(*global));
    for (int i = 0; i < mesh->nodes.num; i++) {
        global[i] = mesh->nodes.global[i];
        ensure(global[i] >= 0);
        if (i < mesh->nodes.num_inner) {
            ensure(beg_global <= global[i] && global[i] < end_global);
        }
        else {
            ensure(!(beg_global <= global[i] && global[i] < end_global));
        }
    }
    ensure(unique(global, mesh->nodes.num, sizeof(*global), compare_long) == mesh->nodes.num);
    teal2_free(global);

    Vector *coord = teal2_calloc(mesh->nodes.num, sizeof(*coord));
    sync2_collect(mesh->nodes.coord, coord, mesh->nodes.global, mesh->nodes.num_inner,
                  mesh->nodes.num, MPI_DOUBLE, 3);
    for (int i = 0; i < mesh->nodes.num; i++) {
        ensure(isclose(mesh->nodes.coord[i].x, coord[i].x));
        ensure(isclose(mesh->nodes.coord[i].y, coord[i].y));
        ensure(isclose(mesh->nodes.coord[i].z, coord[i].z));
    }
    teal2_free(coord);

    validate_unique(mesh->nodes.coord, mesh->nodes.num);
}

static void validate_graph(Graph graph, int num, int min_deg, int max_deg, int min_idx, int max_idx)
{
    ensure(graph.off[0] == 0);
    for (int i = 0; i < num; i++) {
        ensure(graph.off[i] <= graph.off[i + 1]);
        if (min_deg >= 0) {
            ensure(graph.off[i + 1] - graph.off[i] >= min_deg);
        }
        if (max_deg >= 0) {
            ensure(graph.off[i + 1] - graph.off[i] <= max_deg);
        }
        for (int j = graph.off[i]; j < graph.off[i + 1]; j++) {
            ensure(min_idx <= graph.idx[j] && graph.idx[j] <= max_idx);
            for (int k = j + 1; k < graph.off[i + 1]; k++) {
                ensure(graph.idx[j] != graph.idx[k]);
            }
        }
    }
}

static int adjacent(Graph graph, int lhs, int rhs)
{
    for (int i = graph.off[lhs]; i < graph.off[lhs + 1]; i++) {
        if (graph.idx[i] == rhs) {
            return 1;
        }
    }
    return 0;
}

static void validate_cells(const Mesh *mesh)
{
    ensure(mesh->cells.num > 0);
    ensure(mesh->cells.num_inner > 0);
    ensure(mesh->cells.num_inner <= mesh->cells.off_boundary);
    ensure(mesh->cells.off_boundary <= mesh->cells.off_periodic);
    ensure(mesh->cells.off_periodic <= mesh->cells.num);

    validate_graph(mesh->cells.node, mesh->cells.num, 3, MAX_CELL_NODES, 0, mesh->nodes.num - 1);
    validate_graph(mesh->cells.cell, mesh->cells.num, 1, MAX_CELL_FACES, 0, mesh->cells.num - 1);

    for (int i = 0; i < mesh->cells.num; i++) {
        for (int j = mesh->cells.cell.off[i]; j < mesh->cells.cell.off[i + 1]; j++) {
            ensure(mesh->cells.cell.idx[j] != i);
            ensure(adjacent(mesh->cells.cell, mesh->cells.cell.idx[j], i));
        }
    }

    for (int i = 0; i < mesh->cells.num; i++) {
        if (i < mesh->cells.num_inner) {
            ensure(!isclose(mesh->cells.volume[i], 0) && mesh->cells.volume[i] > 0);
        }
        else {
            ensure(isclose(mesh->cells.volume[i], 0));
        }
    }

    validate_unique(mesh->cells.center, mesh->cells.num);
}

static void validate_faces(const Mesh *mesh)
{
    ensure(mesh->faces.num > 0);
    ensure(mesh->faces.num_inner >= 0);
    ensure(mesh->faces.num_inner <= mesh->faces.off_boundary);
    ensure(mesh->faces.off_boundary <= mesh->faces.num);

    validate_graph(mesh->faces.node, mesh->faces.num, 3, MAX_FACE_NODES, 0, mesh->nodes.num - 1);

    for (int i = 0; i < mesh->faces.num; i++) {
        int left = mesh->faces.cell_idx[i].left;
        int right = mesh->faces.cell_idx[i].right;
        ensure(left < right);
        ensure(0 <= left && left < mesh->cells.num_inner);
        if (i < mesh->faces.num_inner) {
            ensure(0 <= right && right < mesh->cells.num_inner);
        }
        else if (i < mesh->faces.off_boundary) {
            ensure(mesh->cells.num_inner <= right && right < mesh->cells.off_boundary);
        }
        else {
            ensure(mesh->cells.off_boundary <= right && right < mesh->cells.num);
        }
        ensure(adjacent(mesh->cells.cell, left, right));
        ensure(adjacent(mesh->cells.cell, right, left));
    }

    for (int i = 0; i < mesh->faces.num; i++) {
        ensure(!isclose(mesh->faces.area[i], 0) && mesh->faces.area[i] > 0);
    }

    validate_unique(mesh->faces.center, mesh->faces.num);

    for (int i = 0; i < mesh->faces.num; i++) {
        ensure(isclose(vector2_norm(mesh->faces.basis[i].x), 1));
        ensure(isclose(vector2_norm(mesh->faces.basis[i].y), 1));
        ensure(isclose(vector2_norm(mesh->faces.basis[i].z), 1));
        ensure(isclose(vector2_dot(mesh->faces.basis[i].x, mesh->faces.basis[i].y), 0));
        ensure(isclose(vector2_dot(mesh->faces.basis[i].y, mesh->faces.basis[i].z), 0));
        ensure(isclose(vector2_dot(mesh->faces.basis[i].z, mesh->faces.basis[i].x), 0));
        int left = mesh->faces.cell_idx[i].left;
        int right = mesh->faces.cell_idx[i].right;
        Vector delta = vector2_sub(mesh->cells.center[right], mesh->cells.center[left]);
        double normal_projection = vector2_dot(mesh->faces.basis[i].x, delta);
        ensure(!isclose(normal_projection, 0) && normal_projection > 0);
    }
}

static void validate_entities(const Mesh *mesh)
{
    ensure(mesh->entities.num > 0);
    ensure(mesh->entities.num_inner > 0);
    ensure(mesh->entities.num_inner <= mesh->entities.off_boundary);
    ensure(mesh->entities.off_boundary <= mesh->entities.num);

    for (int i = 0; i < mesh->entities.num; i++) {
        ensure(strlen(mesh->entities.name[i]) > 0);
    }

    ensure(mesh->entities.cell_off[0] == 0);
    for (int i = 0; i < mesh->entities.num; i++) {
        ensure(mesh->entities.cell_off[i] <= mesh->entities.cell_off[i + 1]);
    }
    ensure(mesh->entities.cell_off[mesh->entities.num] == mesh->cells.off_periodic);
    ensure(mesh->entities.cell_off[mesh->entities.num_inner] == mesh->cells.num_inner);
    ensure(mesh->entities.cell_off[mesh->entities.off_boundary] == mesh->cells.off_boundary);

    for (int i = 0; i < mesh->entities.num; i++) {
        if (i < mesh->entities.num_inner) {
            ensure(mesh->entities.face_off[i] == mesh->faces.num_inner);
        }
        ensure(mesh->entities.face_off[i] <= mesh->entities.face_off[i + 1]);
    }
    ensure(mesh->entities.face_off[mesh->entities.num] <= mesh->faces.num);
    ensure(mesh->entities.face_off[mesh->entities.num_inner] == mesh->faces.num_inner);
    ensure(mesh->entities.face_off[mesh->entities.off_boundary] == mesh->faces.off_boundary);

    for (int i = mesh->entities.num_inner; i < mesh->entities.num; i++) {
        int num_cells = mesh->entities.cell_off[i + 1] - mesh->entities.cell_off[i];
        int num_faces = mesh->entities.face_off[i + 1] - mesh->entities.face_off[i];
        ensure(num_cells == num_faces);
    }

    for (int i = 0; i < mesh->entities.num; i++) {
        if (i < mesh->entities.off_boundary) {
            ensure(isclose(matrix_determinant(mesh->entities.rotation[i]), 0));
            ensure(isclose(vector2_norm(mesh->entities.translation[i]), 0));
        }
        else {
            ensure(isclose(matrix_determinant(mesh->entities.rotation[i]), 1));
        }
    }
}

static void validate_neighbors(const Mesh *mesh)
{
    ensure(mesh->neighbors.num >= 0);

    int *rank = teal2_calloc(mesh->neighbors.num, sizeof(*rank));
    MPI_Request *req_recv = teal2_calloc(mesh->neighbors.num, sizeof(*req_recv));
    MPI_Request *req_send = teal2_calloc(mesh->neighbors.num, sizeof(*req_send));
    for (int i = 0; i < mesh->neighbors.num; i++) {
        ensure(mesh->neighbors.rank[i] >= 0);
        ensure(mesh->neighbors.rank[i] < sync2.size);
        MPI_Irecv(&rank[i], 1, MPI_INT, mesh->neighbors.rank[i], mesh->neighbors.tag[i][0],
                  sync2.comm, &req_recv[i]);
        MPI_Isend(&sync2.rank, 1, MPI_INT, mesh->neighbors.rank[i], mesh->neighbors.tag[i][1],
                  sync2.comm, &req_send[i]);
    }
    MPI_Waitall(mesh->neighbors.num, req_recv, MPI_STATUSES_IGNORE);
    MPI_Waitall(mesh->neighbors.num, req_send, MPI_STATUSES_IGNORE);
    for (int i = 0; i < mesh->neighbors.num; i++) {
        ensure(mesh->neighbors.rank[i] == rank[i]);
    }
    teal2_free(req_recv);
    teal2_free(req_send);
    teal2_free(rank);

    ensure(mesh->neighbors.recv_off[0] == mesh->cells.off_boundary);
    for (int i = 0; i < mesh->neighbors.num; i++) {
        ensure(mesh->neighbors.recv_off[i] <= mesh->neighbors.recv_off[i + 1]);
    }
    ensure(mesh->neighbors.recv_off[mesh->neighbors.num] == mesh->cells.num);

    validate_graph(mesh->neighbors.send, mesh->neighbors.num, -1, -1, 0, mesh->cells.num_inner - 1);
}

void mesh2_validate(const Mesh *mesh)
{
    assert(mesh);
    validate_nodes(mesh);
    validate_cells(mesh);
    validate_faces(mesh);
    validate_entities(mesh);
    validate_neighbors(mesh);
}
