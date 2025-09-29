#include <math.h>
#include <stdio.h>
#include <string.h>

#include "mesh.h"
#include "teal/arena.h"
#include "teal/dict.h"
#include "teal/sync.h"
#include "teal/utils.h"
#include "teal/vector.h"

#define test(expr) ((expr) ? (void)0 : test_fail(__FILE__, __LINE__, __func__, #expr))

static void test_fail(const char *file, long line, const char *func, const char *expr)
{
    fprintf(stderr, "[%d] %s:%ld: %s: Test `%s` failed.\n", sync.rank, file, line, func, expr);
}

static void test_nodes(const MeshNodes *nodes)
{
    Arena save = arena_save();

    test(nodes->num > 0);
    test(nodes->num_inner >= 0);
    test(nodes->num_inner <= nodes->num);

    if (nodes->global) {
        Dict *global = dict_create(sizeof(long), 0);
        for (long i = 0; i < nodes->num; i++) {
            test(nodes->global[i] >= 0);
            dict_insert(global, &nodes->global[i], 0);
        }
        test(global->num == nodes->num);
    }

    if (nodes->coord) {
        for (long i = 0; i < nodes->num; i++) {
            test(isfinite(nodes->coord[i].x));
            test(isfinite(nodes->coord[i].y));
            test(isfinite(nodes->coord[i].z));
        }
    }

    arena_load(save);
}

static void test_graph(const MeshGraph *graph, long num, long min_deg, long max_deg, long min_idx,
                       long max_idx)
{
    test(graph->off[0] == 0);
    for (long i = 0; i < num; i++) {
        test(graph->off[i] <= graph->off[i + 1]);
        if (min_deg >= 0) {
            test(graph->off[i + 1] - graph->off[i] >= min_deg);
        }
        if (max_deg >= 0) {
            test(graph->off[i + 1] - graph->off[i] <= max_deg);
        }
        for (long j = graph->off[i]; j < graph->off[i + 1]; j++) {
            test(graph->idx[j] >= min_idx);
            test(graph->idx[j] <= max_idx);
        }
    }
}

static bool has_neighbor(const MeshGraph *graph, long lhs, long rhs)
{
    for (long i = graph->off[lhs]; i < graph->off[lhs + 1]; i++) {
        if (graph->idx[i] == rhs) {
            return true;
        }
    }
    return false;
}

static void test_cells(const MeshNodes *nodes, const MeshCells *cells)
{
    Arena save = arena_save();

    test(cells->num > 0);
    test(cells->num_inner > 0);
    test(cells->num_inner < cells->num);

    long num_outer = cells->num - cells->num_inner;
    test(cells->num_ghost >= 0);
    test(cells->num_ghost <= num_outer);
    test(cells->num_periodic >= 0);
    test(cells->num_periodic <= num_outer);
    test(cells->num_ghost + cells->num_periodic <= num_outer);

    if (cells->node.off && cells->node.idx) {
        test_graph(&cells->node, cells->num, 3, MAX_CELL_NODES, 0, nodes->num - 1);
    }

    if (cells->cell.off && cells->cell.idx) {
        test_graph(&cells->cell, cells->num, 1, MAX_CELL_FACES, 0, cells->num - 1);

        Dict *local = dict_create(sizeof(long), 0);
        for (long i = 0; i < cells->num; i++) {
            for (long j = cells->cell.off[i]; j < cells->cell.off[i + 1]; j++) {
                test(cells->cell.idx[j] != i);
                test(has_neighbor(&cells->cell, cells->cell.idx[j], i));
                dict_insert(local, &cells->cell.idx[j], 0);
            }
        }
        test(local->num == cells->num);
    }

    if (cells->volume) {
        for (long i = 0; i < cells->num; i++) {
            test(isfinite(cells->volume[i]));
            if (i < cells->num_inner) {
                test(!isclose(cells->volume[i], 0) && cells->volume[i] > 0);
            }
            else {
                test(isclose(cells->volume[i], 0));
            }
        }
    }

    if (cells->center) {
        for (long i = 0; i < cells->num; i++) {
            test(isfinite(cells->center[i].x));
            test(isfinite(cells->center[i].y));
            test(isfinite(cells->center[i].z));
        }
    }

    if (cells->projection) {
        for (long i = 0; i < cells->num; i++) {
            test(isfinite(cells->projection[i].x));
            test(isfinite(cells->projection[i].y));
            test(isfinite(cells->projection[i].z));
            if (i < cells->num_inner) {
                test(!isclose(cells->projection[i].x, 0) && cells->projection[i].x > 0);
                test(!isclose(cells->projection[i].y, 0) && cells->projection[i].y > 0);
                test(!isclose(cells->projection[i].z, 0) && cells->projection[i].z > 0);
            }
            else {
                test(isclose(cells->projection[i].x, 0));
                test(isclose(cells->projection[i].y, 0));
                test(isclose(cells->projection[i].z, 0));
            }
        }
    }

    arena_load(save);
}

static void test_faces(const MeshNodes *nodes, const MeshCells *cells, const MeshFaces *faces)
{
    Arena save = arena_save();

    test(faces->num > 0);
    test(faces->num_inner > 0);
    test(faces->num_inner < faces->num);
    test(faces->num_ghost >= 0);
    test(faces->num_ghost <= faces->num - faces->num_inner);

    if (faces->node.off && faces->node.idx) {
        test_graph(&faces->node, faces->num, 3, MAX_FACE_NODES, 0, nodes->num - 1);
    }

    if (faces->cell) {
        Dict *left = dict_create(sizeof(long), 0);
        Dict *right = dict_create(sizeof(long), 0);
        for (long i = 0; i < faces->num; i++) {
            test(faces->cell[i].left >= 0);
            test(faces->cell[i].left < cells->num_inner);
            if (i < faces->num_inner) {
                test(faces->cell[i].right >= 0);
                test(faces->cell[i].right < cells->num_inner);
            }
            else if (i < faces->num_inner + faces->num_ghost) {
                test(faces->cell[i].right >= cells->num_inner);
                test(faces->cell[i].right < cells->num_inner + cells->num_ghost);
            }
            else {
                test(faces->cell[i].right >= cells->num_inner + cells->num_ghost);
                test(faces->cell[i].right < cells->num);
            }
            test(faces->cell[i].left != faces->cell[i].right);
            dict_insert(left, &faces->cell[i].left, 0);
            dict_insert(right, &faces->cell[i].right, 0);
        }
        test(left->num <= cells->num_inner);
        test(right->num >= cells->num - cells->num_inner);
    }

    if (faces->area) {
        for (long i = 0; i < faces->num; i++) {
            test(isfinite(faces->area[i]));
            test(!isclose(faces->area[i], 0) && faces->area[i] > 0);
        }
    }

    if (faces->center) {
        for (long i = 0; i < faces->num; i++) {
            test(isfinite(faces->center[i].x));
            test(isfinite(faces->center[i].y));
            test(isfinite(faces->center[i].z));
        }
    }

    if (faces->basis) {
        for (long i = 0; i < faces->num; i++) {
            test(isfinite(faces->basis[i].x.x));
            test(isfinite(faces->basis[i].x.y));
            test(isfinite(faces->basis[i].x.z));
            test(isfinite(faces->basis[i].y.x));
            test(isfinite(faces->basis[i].y.y));
            test(isfinite(faces->basis[i].y.z));
            test(isfinite(faces->basis[i].z.x));
            test(isfinite(faces->basis[i].z.y));
            test(isfinite(faces->basis[i].z.z));
            test(isclose(vector_len(faces->basis[i].x), 1));
            test(isclose(vector_len(faces->basis[i].y), 1));
            test(isclose(vector_len(faces->basis[i].z), 1));
            test(isclose(vector_dot(faces->basis[i].x, faces->basis[i].y), 0));
            test(isclose(vector_dot(faces->basis[i].y, faces->basis[i].z), 0));
            test(isclose(vector_dot(faces->basis[i].z, faces->basis[i].x), 0));
        }
    }

    if (faces->cell && faces->basis && cells->center) {
        for (long i = 0; i < faces->num; i++) {
            long left = faces->cell[i].left;
            long right = faces->cell[i].right;
            vector normal = faces->basis[i].x;
            vector l2r = vector_sub(cells->center[right], cells->center[left]);
            test(!isclose(vector_len(l2r), 0) && vector_len(l2r) > 0);
            test(!isclose(vector_dot(normal, l2r), 0) && vector_dot(normal, l2r) > 0);
        }
    }

    if (faces->weight) {
        for (long i = 0; i < faces->num; i++) {
            test(isfinite(faces->weight[i].x));
            test(isfinite(faces->weight[i].y));
            test(isfinite(faces->weight[i].z));
        }
    }

    arena_load(save);
}

static void test_entities(const MeshCells *cells, const MeshFaces *faces,
                          const MeshEntities *entities)
{
    test(entities->num > 0);
    test(entities->num_inner > 0);
    test(entities->num_inner <= entities->num);
    test(entities->num_ghost >= 0);
    test(entities->num_ghost <= entities->num - entities->num_inner);

    if (entities->name) {
        for (long i = 0; i < entities->num; i++) {
            test(entities->name[i] && strlen(entities->name[i]) > 0);
        }
    }

    if (entities->cell_off) {
        test(entities->cell_off[0] == 0);
        for (long i = 0; i < entities->num; i++) {
            test(entities->cell_off[i] <= entities->cell_off[i + 1]);
        }
        test(entities->cell_off[entities->num_inner] == cells->num_inner);
        test(entities->cell_off[entities->num_inner + entities->num_ghost] ==
             cells->num_inner + cells->num_ghost);
        test(entities->cell_off[entities->num] <= cells->num);
    }

    if (entities->face_off) {
        test(entities->face_off[0] == faces->num_inner);
        for (long i = 0; i < entities->num; i++) {
            test(entities->face_off[i] <= entities->face_off[i + 1]);
        }
        test(entities->face_off[entities->num_inner] == faces->num_inner);
        test(entities->face_off[entities->num_inner + entities->num_ghost] ==
             faces->num_inner + faces->num_ghost);
        test(entities->face_off[entities->num] <= faces->num);

        if (entities->cell_off) {
            for (long i = entities->num_inner; i < entities->num; i++) {
                long num_cells = entities->cell_off[i + 1] - entities->cell_off[i];
                long num_faces = entities->face_off[i + 1] - entities->face_off[i];
                test(num_cells == num_faces);
            }
        }
    }

    if (entities->offset) {
        for (long i = 0; i < entities->num; i++) {
            test(isfinite(entities->offset[i].x));
            test(isfinite(entities->offset[i].y));
            test(isfinite(entities->offset[i].z));
            if (i < entities->num_inner + entities->num_ghost) {
                test(isclose(entities->offset[i].x, 0));
                test(isclose(entities->offset[i].y, 0));
                test(isclose(entities->offset[i].z, 0));
            }
            else {
                test(!isclose(vector_len(entities->offset[i]), 0));
                if (entities->name) {
                    test(strchr(entities->name[i], ':'));
                }
            }
        }
    }
}

static void test_neighbors(const MeshCells *cells, const MeshNeighbors *neighbors)
{
    Arena save = arena_save();

    test(neighbors->num >= 0);

    if (neighbors->rank) {
        int *rank = arena_malloc(neighbors->num, sizeof(*rank));
        int tag = sync_tag();
        MPI_Request *req = arena_malloc(neighbors->num, sizeof(*req));
        for (long i = 0; i < neighbors->num; i++) {
            test(neighbors->rank[i] >= 0);
            test(neighbors->rank[i] < sync.size);
            MPI_Isendrecv(&sync.rank, 1, MPI_INT, neighbors->rank[i], tag, &rank[i], 1, MPI_INT,
                          neighbors->rank[i], tag, sync.comm, &req[i]);
        }
        MPI_Waitall(neighbors->num, req, MPI_STATUSES_IGNORE);
        for (long i = 0; i < neighbors->num; i++) {
            test(neighbors->rank[i] == rank[i]);
        }
    }

    if (neighbors->recv_off) {
        test(neighbors->recv_off[0] == cells->num_inner + cells->num_ghost);
        for (long i = 0; i < neighbors->num; i++) {
            test(neighbors->recv_off[i] <= neighbors->recv_off[i + 1]);
        }
        test(neighbors->recv_off[neighbors->num] == cells->num);
    }

    if (neighbors->send.off && neighbors->send.idx) {
        test_graph(&neighbors->send, neighbors->num, -1, -1, 0, cells->num_inner - 1);
    }

    arena_load(save);
}

void mesh_test(const Mesh *mesh)
{
    assert(mesh);
    test_nodes(&mesh->nodes);
    test_cells(&mesh->nodes, &mesh->cells);
    test_faces(&mesh->nodes, &mesh->cells, &mesh->faces);
    test_entities(&mesh->cells, &mesh->faces, &mesh->entities);
    test_neighbors(&mesh->cells, &mesh->neighbors);
}
