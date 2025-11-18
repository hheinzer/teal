#include "teal/check.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "mesh.h"
#include "teal/arena.h"
#include "teal/array.h"
#include "teal/assert.h"
#include "teal/sync.h"
#include "teal/utils.h"
#include "teal/vector.h"

static void test_nodes(const MeshNodes *nodes)
{
    Arena save = arena_save();

    check(nodes->num > 0);
    check(nodes->num_inner >= 0);
    check(nodes->num_inner <= nodes->num);

    if (nodes->global) {
        int *global = arena_malloc(nodes->num, sizeof(*global));
        for (int i = 0; i < nodes->num; i++) {
            check(nodes->global[i] >= 0);
            global[i] = nodes->global[i];
        }

        qsort(global, nodes->num, sizeof(*global), cmp_int);
        for (int i = 1; i < nodes->num; i++) {
            check(global[i - 1] != global[i]);
        }
    }

    if (nodes->coord) {
        for (int i = 0; i < nodes->num; i++) {
            check(isfinite(nodes->coord[i].x));
            check(isfinite(nodes->coord[i].y));
            check(isfinite(nodes->coord[i].z));
        }
    }

    arena_load(save);
}

static void test_graph(Graph graph, int num, int min_deg, int max_deg, int min_idx, int max_idx)
{
    check(graph.off[0] == 0);
    for (int i = 0; i < num; i++) {
        check(graph.off[i] <= graph.off[i + 1]);
        if (min_deg >= 0) {
            check(graph.off[i + 1] - graph.off[i] >= min_deg);
        }
        if (max_deg >= 0) {
            check(graph.off[i + 1] - graph.off[i] <= max_deg);
        }
        for (int j = graph.off[i]; j < graph.off[i + 1]; j++) {
            check(graph.idx[j] >= min_idx);
            check(graph.idx[j] <= max_idx);
        }
    }
}

static bool has_neighbor(Graph graph, int lhs, int rhs)
{
    for (int i = graph.off[lhs]; i < graph.off[lhs + 1]; i++) {
        if (graph.idx[i] == rhs) {
            return true;
        }
    }
    return false;
}

static void test_cells(const MeshNodes *nodes, const MeshCells *cells)
{
    Arena save = arena_save();

    check(cells->num_inner > 0);
    check(cells->num_inner <= cells->off_ghost);
    check(cells->off_ghost <= cells->off_periodic);
    check(cells->off_periodic <= cells->num);

    if (cells->node.off && cells->node.idx) {
        test_graph(cells->node, cells->num, 3, MAX_CELL_NODES, 0, nodes->num - 1);
    }

    if (cells->cell.off && cells->cell.idx) {
        test_graph(cells->cell, cells->num, 1, MAX_CELL_FACES, 0, cells->num - 1);

        int *local = arena_calloc(cells->num, sizeof(*local));
        for (int i = 0; i < cells->num; i++) {
            for (int j = cells->cell.off[i]; j < cells->cell.off[i + 1]; j++) {
                check(cells->cell.idx[j] != i);
                check(has_neighbor(cells->cell, cells->cell.idx[j], i));
                local[cells->cell.idx[j]] = 1;
            }
        }
        check(array_lsum(local, cells->num) == cells->num);
    }

    if (cells->volume) {
        for (int i = 0; i < cells->num; i++) {
            check(isfinite(cells->volume[i]));
            if (i < cells->num_inner) {
                check(!isclose(cells->volume[i], 0) && cells->volume[i] > 0);
            }
            else {
                check(isclose(cells->volume[i], 0));
            }
        }
    }

    if (cells->center) {
        for (int i = 0; i < cells->num; i++) {
            check(isfinite(cells->center[i].x));
            check(isfinite(cells->center[i].y));
            check(isfinite(cells->center[i].z));
        }
    }

    if (cells->projection) {
        for (int i = 0; i < cells->num; i++) {
            check(isfinite(cells->projection[i].x));
            check(isfinite(cells->projection[i].y));
            check(isfinite(cells->projection[i].z));
            if (i < cells->num_inner) {
                check(!isclose(cells->projection[i].x, 0) && cells->projection[i].x > 0);
                check(!isclose(cells->projection[i].y, 0) && cells->projection[i].y > 0);
                check(!isclose(cells->projection[i].z, 0) && cells->projection[i].z > 0);
            }
            else {
                check(isclose(cells->projection[i].x, 0));
                check(isclose(cells->projection[i].y, 0));
                check(isclose(cells->projection[i].z, 0));
            }
        }
    }

    arena_load(save);
}

static void test_faces(const MeshNodes *nodes, const MeshCells *cells, const MeshFaces *faces)
{
    Arena save = arena_save();

    check(faces->num_inner > 0);
    check(faces->num_inner <= faces->off_ghost);
    check(faces->off_ghost <= faces->num);

    if (faces->node.off && faces->node.idx) {
        test_graph(faces->node, faces->num, 3, MAX_FACE_NODES, 0, nodes->num - 1);
    }

    if (faces->cell) {
        int *left = arena_calloc(cells->num_inner, sizeof(*left));
        int *right = arena_calloc(cells->num, sizeof(*right));
        for (int i = 0; i < faces->num; i++) {
            check(faces->cell[i].left >= 0);
            check(faces->cell[i].left < cells->num_inner);
            if (i < faces->num_inner) {
                check(faces->cell[i].right >= 0);
                check(faces->cell[i].right < cells->num_inner);
            }
            else if (i < faces->off_ghost) {
                check(faces->cell[i].right >= cells->num_inner);
                check(faces->cell[i].right < cells->off_ghost);
            }
            else {
                check(faces->cell[i].right >= cells->off_ghost);
                check(faces->cell[i].right < cells->num);
            }
            check(faces->cell[i].left != faces->cell[i].right);
            left[faces->cell[i].left] = 1;
            right[faces->cell[i].right] = 1;
        }
        check(array_lsum(left, cells->num_inner) <= cells->num_inner);
        check(array_lsum(right, cells->num) >= cells->num - cells->num_inner);
    }

    if (faces->area) {
        for (int i = 0; i < faces->num; i++) {
            check(isfinite(faces->area[i]));
            check(!isclose(faces->area[i], 0) && faces->area[i] > 0);
        }
    }

    if (faces->center) {
        for (int i = 0; i < faces->num; i++) {
            check(isfinite(faces->center[i].x));
            check(isfinite(faces->center[i].y));
            check(isfinite(faces->center[i].z));
        }
    }

    if (faces->basis) {
        for (int i = 0; i < faces->num; i++) {
            check(isfinite(faces->basis[i].x.x));
            check(isfinite(faces->basis[i].x.y));
            check(isfinite(faces->basis[i].x.z));
            check(isfinite(faces->basis[i].y.x));
            check(isfinite(faces->basis[i].y.y));
            check(isfinite(faces->basis[i].y.z));
            check(isfinite(faces->basis[i].z.x));
            check(isfinite(faces->basis[i].z.y));
            check(isfinite(faces->basis[i].z.z));
            check(isclose(vector_norm(faces->basis[i].x), 1));
            check(isclose(vector_norm(faces->basis[i].y), 1));
            check(isclose(vector_norm(faces->basis[i].z), 1));
            check(isclose(vector_dot(faces->basis[i].x, faces->basis[i].y), 0));
            check(isclose(vector_dot(faces->basis[i].y, faces->basis[i].z), 0));
            check(isclose(vector_dot(faces->basis[i].z, faces->basis[i].x), 0));
        }
    }

    if (faces->cell && faces->basis && cells->center) {
        for (int i = 0; i < faces->num; i++) {
            int left = faces->cell[i].left;
            int right = faces->cell[i].right;
            vector normal = faces->basis[i].x;
            vector l2r = vector_sub(cells->center[right], cells->center[left]);
            check(!isclose(vector_norm(l2r), 0) && vector_norm(l2r) > 0);
            check(!isclose(vector_dot(normal, l2r), 0) && vector_dot(normal, l2r) > 0);
        }
    }

    if (faces->weight) {
        for (int i = 0; i < faces->num; i++) {
            check(isfinite(faces->weight[i].x));
            check(isfinite(faces->weight[i].y));
            check(isfinite(faces->weight[i].z));
        }
    }

    arena_load(save);
}

static void test_entities(const MeshCells *cells, const MeshFaces *faces,
                          const MeshEntities *entities)
{
    check(entities->num > 0);
    check(entities->num_inner <= entities->off_ghost);
    check(entities->off_ghost <= entities->num);

    if (entities->name) {
        for (int i = 0; i < entities->num; i++) {
            check(strlen(entities->name[i]) > 0);
        }
    }

    if (entities->cell_off) {
        check(entities->cell_off[0] == 0);
        for (int i = 0; i < entities->num; i++) {
            check(entities->cell_off[i] <= entities->cell_off[i + 1]);
        }
        check(entities->cell_off[entities->num_inner] == cells->num_inner);
        check(entities->cell_off[entities->off_ghost] == cells->off_ghost);
        check(entities->cell_off[entities->num] <= cells->num);
    }

    if (entities->face_off) {
        check(entities->face_off[0] == faces->num_inner);
        for (int i = 0; i < entities->num; i++) {
            check(entities->face_off[i] <= entities->face_off[i + 1]);
        }
        check(entities->face_off[entities->num_inner] == faces->num_inner);
        check(entities->face_off[entities->off_ghost] == faces->off_ghost);
        check(entities->face_off[entities->num] <= faces->num);

        if (entities->cell_off) {
            for (int i = entities->num_inner; i < entities->num; i++) {
                int num_cells = entities->cell_off[i + 1] - entities->cell_off[i];
                int num_faces = entities->face_off[i + 1] - entities->face_off[i];
                check(num_cells == num_faces);
            }
        }
    }

    if (entities->translation) {
        for (int i = 0; i < entities->num; i++) {
            check(isfinite(entities->translation[i].x));
            check(isfinite(entities->translation[i].y));
            check(isfinite(entities->translation[i].z));
            if (i < entities->off_ghost) {
                check(isclose(entities->translation[i].x, 0));
                check(isclose(entities->translation[i].y, 0));
                check(isclose(entities->translation[i].z, 0));
            }
            else {
                check(!isclose(vector_norm(entities->translation[i]), 0));
                if (entities->name) {
                    check(strchr(entities->name[i], ':'));
                }
            }
        }
    }
}

static void test_neighbors(const MeshCells *cells, const MeshNeighbors *neighbors)
{
    Arena save = arena_save();

    check(neighbors->num >= 0);

    if (neighbors->rank) {
        int *rank = arena_malloc(neighbors->num, sizeof(*rank));
        int tag = sync_tag();
        MPI_Request *req = arena_malloc(neighbors->num, sizeof(*req));
        for (int i = 0; i < neighbors->num; i++) {
            check(neighbors->rank[i] >= 0);
            check(neighbors->rank[i] < sync.size);
            MPI_Isendrecv(&sync.rank, 1, MPI_INT, neighbors->rank[i], tag, &rank[i], 1, MPI_INT,
                          neighbors->rank[i], tag, sync.comm, &req[i]);
        }
        MPI_Waitall(neighbors->num, req, MPI_STATUSES_IGNORE);
        for (int i = 0; i < neighbors->num; i++) {
            check(neighbors->rank[i] == rank[i]);
        }
    }

    if (neighbors->recv_off) {
        check(neighbors->recv_off[0] == cells->off_ghost);
        for (int i = 0; i < neighbors->num; i++) {
            check(neighbors->recv_off[i] <= neighbors->recv_off[i + 1]);
        }
        check(neighbors->recv_off[neighbors->num] == cells->num);
    }

    if (neighbors->send.off && neighbors->send.idx) {
        test_graph(neighbors->send, neighbors->num, -1, -1, 0, cells->num_inner - 1);
    }

    arena_load(save);
}

void mesh_check(const Mesh *mesh)
{
    assert(mesh);
    test_nodes(&mesh->nodes);
    test_cells(&mesh->nodes, &mesh->cells);
    test_faces(&mesh->nodes, &mesh->cells, &mesh->faces);
    test_entities(&mesh->cells, &mesh->faces, &mesh->entities);
    test_neighbors(&mesh->cells, &mesh->neighbors);
}
