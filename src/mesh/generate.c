#include <math.h>

#include "connectivity.h"
#include "core/array.h"
#include "core/dict.h"
#include "core/memory.h"
#include "core/sync.h"
#include "core/utils.h"
#include "mesh.h"
#include "partition.h"

static void compute_face_geometry(Mesh *mesh);
static void correct_face_node_order(const Mesh *mesh, long *node, long n_nodes);
static void compute_face_area(const double (*x)[N_DIMS], double *area, long n_nodes);
static void compute_face_center(const double (*x)[N_DIMS], double *center, long n_nodes);
static void compute_face_normal(const double (*x)[N_DIMS], double *normal, long n_nodes);
static void correct_face_normal_direction(const Mesh *mesh, const long *cell, const double *center,
                                          double *normal);
static void compute_face_tangents(double *normal);

static void compute_cell_geometry(Mesh *mesh);
static void compute_cell_volume(const double (*x)[N_DIMS], double *volume, long n_nodes);
static void compute_cell_center(const double (*x)[N_DIMS], double *center, long n_nodes);
static void compute_cell_projections(Mesh *mesh);

static void compute_ghost_cell_centers(const Mesh *mesh);

static void find_periodic_connections(const Mesh *mesh, Dict *periodic);
static void compute_face_reconstruction(Mesh *mesh, const Dict *periodic);
static void compute_cell_reconstruction(Mesh *mesh, const Dict *periodic);

void mesh_generate(Mesh *mesh)
{
    partition(mesh);
    connectivity_faces(mesh);

    compute_face_geometry(mesh);
    compute_cell_geometry(mesh);
    compute_cell_projections(mesh);

    sync_all(mesh, N_DIMS, mesh->cell.center);
    compute_ghost_cell_centers(mesh);

    fcleanup(dict_free) Dict periodic = dict_create(mesh->n_ghost_cells);
    find_periodic_connections(mesh, &periodic);
    compute_face_reconstruction(mesh, &periodic);
    compute_cell_reconstruction(mesh, &periodic);

    mesh->volume = sync_sum(array_sum(mesh->cell.volume, mesh->n_inner_cells));
}

static void compute_face_geometry(Mesh *mesh)
{
    const ALIAS(i_node, mesh->face.i_node);
    const ALIAS(coord, mesh->node.coord);
    const ALIAS(cell, mesh->face.cell);
    ALIAS(node, mesh->face.node);
    double *area = memory_calloc(mesh->n_faces, sizeof(*area));
    double(*center)[N_DIMS] = memory_calloc(mesh->n_faces, sizeof(*center));
    double(*normal)[N_DIMS * N_DIMS] = memory_calloc(mesh->n_faces, sizeof(*normal));
    for (long j = 0; j < mesh->n_faces; ++j) {
        const long n_nodes = i_node[j + 1] - i_node[j];
        correct_face_node_order(mesh, &node[i_node[j]], n_nodes);

        double x[MAX_FACE_NODES][N_DIMS] = {0};
        for (long n = 0, i = i_node[j]; i < i_node[j + 1]; ++i, ++n)
            for (long d = 0; d < N_DIMS; ++d) x[n][d] = coord[node[i]][d];

        compute_face_area(x, &area[j], n_nodes);
        compute_face_center(x, center[j], n_nodes);
        compute_face_normal(x, normal[j], n_nodes);
        correct_face_normal_direction(mesh, cell[j], center[j], normal[j]);
        compute_face_tangents(normal[j]);
    }
    mesh->face.area = area;
    mesh->face.center = center;
    mesh->face.normal = normal;
}

static void correct_face_node_order(const Mesh *mesh, long *node, long n_nodes)
{
    const ALIAS(x, mesh->node.coord);
    switch (n_nodes) {
        case 3: break;  // node order of triangle does not need to be corrected
        case 4: {
            double u[N_DIMS], v[N_DIMS], w[N_DIMS], n0[N_DIMS], n1[N_DIMS];
            for (long d = 0; d < N_DIMS; ++d) {
                u[d] = x[node[1]][d] - x[node[0]][d];
                v[d] = x[node[2]][d] - x[node[0]][d];
                w[d] = x[node[3]][d] - x[node[0]][d];
            }
            array_cross(u, v, n0);
            array_cross(v, w, n1);
            if (array_dot(n0, n1, N_DIMS) > 0) return;

            const long swap1 = node[1];
            node[1] = node[2];
            node[2] = swap1;
            for (long d = 0; d < N_DIMS; ++d) {
                u[d] = x[node[1]][d] - x[node[0]][d];
                v[d] = x[node[2]][d] - x[node[0]][d];
                w[d] = x[node[3]][d] - x[node[0]][d];
            }
            array_cross(u, v, n0);
            array_cross(v, w, n1);
            if (array_dot(n0, n1, N_DIMS) > 0) return;

            const long swap2 = node[2];
            node[2] = node[3];
            node[3] = swap2;
            for (long d = 0; d < N_DIMS; ++d) {
                u[d] = x[node[1]][d] - x[node[0]][d];
                v[d] = x[node[2]][d] - x[node[0]][d];
                w[d] = x[node[3]][d] - x[node[0]][d];
            }
            array_cross(u, v, n0);
            array_cross(v, w, n1);
            ensure(array_dot(n0, n1, N_DIMS) > 0);
            break;
        }
        default: error("unsupported number of face nodes '%ld'", n_nodes);
    }
}

static void compute_face_area(const double (*x)[N_DIMS], double *area, long n_nodes)
{
    switch (n_nodes) {
        case 3: {
            double AB[N_DIMS], AC[N_DIMS], cross[N_DIMS];
            for (long d = 0; d < N_DIMS; ++d) {
                AB[d] = x[1][d] - x[0][d];
                AC[d] = x[2][d] - x[0][d];
            }
            *area = array_norm(array_cross(AB, AC, cross), N_DIMS) / 2;
            break;
        }
        case 4: {
            double x0[3][N_DIMS], x1[3][N_DIMS], a0, a1;
            for (long d = 0; d < N_DIMS; ++d) {
                x0[0][d] = x[0][d];
                x0[1][d] = x[1][d];
                x0[2][d] = x[2][d];
                x1[0][d] = x[0][d];
                x1[1][d] = x[2][d];
                x1[2][d] = x[3][d];
            }
            compute_face_area(x0, &a0, 3);
            compute_face_area(x1, &a1, 3);
            *area = a0 + a1;
            break;
        }
        default: error("unsupported number of face nodes '%ld'", n_nodes);
    }
    ensure(*area > 0);
}

static void compute_face_center(const double (*x)[N_DIMS], double *center, long n_nodes)
{
    switch (n_nodes) {
        case 3: {
            for (long d = 0; d < N_DIMS; ++d) center[d] = 0;
            for (long i = 0; i < 3; ++i)
                for (long d = 0; d < N_DIMS; ++d) center[d] += x[i][d] / 3;
            break;
        }
        case 4: {
            double x0[3][N_DIMS], x1[3][N_DIMS], a0, a1, c0[N_DIMS], c1[N_DIMS];
            for (long d = 0; d < N_DIMS; ++d) {
                x0[0][d] = x[0][d];
                x0[1][d] = x[1][d];
                x0[2][d] = x[2][d];
                x1[0][d] = x[0][d];
                x1[1][d] = x[2][d];
                x1[2][d] = x[3][d];
            }
            compute_face_area(x0, &a0, 3);
            compute_face_area(x1, &a1, 3);
            compute_face_center(x0, c0, 3);
            compute_face_center(x1, c1, 3);
            for (long d = 0; d < N_DIMS; ++d) center[d] = (a0 * c0[d] + a1 * c1[d]) / (a0 + a1);
            break;
        }
        default: error("unsupported number of face nodes '%ld'", n_nodes);
    }
}

static void compute_face_normal(const double (*x)[N_DIMS], double *normal, long n_nodes)
{
    switch (n_nodes) {
        case 3: {
            double AB[N_DIMS], AC[N_DIMS];
            for (long d = 0; d < N_DIMS; ++d) {
                AB[d] = x[1][d] - x[0][d];
                AC[d] = x[2][d] - x[0][d];
            }
            array_normalize(array_cross(AB, AC, normal), N_DIMS);
            break;
        }
        case 4: {
            double c[N_DIMS];
            compute_face_center(x, c, n_nodes);
            for (long d = 0; d < N_DIMS; ++d) normal[d] = 0;
            for (long i = 0; i < 4; ++i) {
                double xi[3][N_DIMS], ni[N_DIMS];
                for (long d = 0; d < N_DIMS; ++d) {
                    xi[0][d] = c[d];
                    xi[1][d] = x[i][d];
                    xi[2][d] = x[(i + 1) % 4][d];
                }
                compute_face_normal(xi, ni, 3);
                for (long d = 0; d < N_DIMS; ++d) normal[d] += ni[d] / 4;
            }
            array_normalize(normal, N_DIMS);
            break;
        }
        default: error("unsupported number of face nodes '%ld'", n_nodes);
    }
    ensure(isclose(array_norm(normal, N_DIMS), 1));
}

static void correct_face_normal_direction(const Mesh *mesh, const long *cell, const double *center,
                                          double *normal)
{
    const ALIAS(i_node, mesh->cell.i_node);
    const ALIAS(node, mesh->cell.node);
    const ALIAS(coord, mesh->node.coord);
    const long n_nodes = i_node[cell[L] + 1] - i_node[cell[L]];
    double C[N_DIMS] = {0};
    for (long i = i_node[cell[L]]; i < i_node[cell[L] + 1]; ++i)
        for (long d = 0; d < N_DIMS; ++d) C[d] += coord[node[i]][d] / n_nodes;

    double FC[N_DIMS];
    for (long d = 0; d < N_DIMS; ++d) FC[d] = C[d] - center[d];

    const double dot = array_dot(FC, normal, N_DIMS);
    ensure(!isclose(dot, 0));
    if (dot > 0)  // make sure that normal points outwards wrt. inner cell
        for (long d = 0; d < N_DIMS; ++d) normal[d] *= -1;
}

static void compute_face_tangents(double *normal)
{
    // Billett and Toro 1998
    double(*rotate)[N_DIMS] = (void *)normal;
    const double nx = rotate[X][X], ny = rotate[X][Y], nz = rotate[X][Z];
    const double nqz = sqrt(1 - sq(nz));
    const double nqy = sqrt(1 - sq(ny));
    if (nqz > nqy) {
        rotate[Y][X] = -ny / nqz;
        rotate[Y][Y] = nx / nqz;
        rotate[Y][Z] = 0;
        rotate[Z][X] = -nx * nz / nqz;
        rotate[Z][Y] = -ny * nz / nqz;
        rotate[Z][Z] = nqz;
    }
    else {
        rotate[Y][X] = -nx * ny / nqy;
        rotate[Y][Y] = nqy;
        rotate[Y][Z] = -ny * nz / nqy;
        rotate[Z][X] = -nz / nqy;
        rotate[Z][Y] = 0;
        rotate[Z][Z] = nx / nqy;
    }
    ensure(isclose(array_norm(rotate[X], N_DIMS), 1));
    ensure(isclose(array_norm(rotate[Y], N_DIMS), 1));
    ensure(isclose(array_norm(rotate[Z], N_DIMS), 1));
}

static void compute_cell_geometry(Mesh *mesh)
{
    const ALIAS(i_node, mesh->cell.i_node);
    const ALIAS(node, mesh->cell.node);
    const ALIAS(coord, mesh->node.coord);
    double *volume = memory_calloc(mesh->n_inner_cells, sizeof(*volume));
    double(*center)[N_DIMS] = memory_calloc(mesh->n_cells, sizeof(*center));
    for (long j = 0; j < mesh->n_inner_cells; ++j) {
        const long n_nodes = i_node[j + 1] - i_node[j];

        double x[MAX_CELL_NODES][N_DIMS] = {0};
        for (long n = 0, i = i_node[j]; i < i_node[j + 1]; ++i, ++n)
            for (long d = 0; d < N_DIMS; ++d) x[n][d] = coord[node[i]][d];

        compute_cell_volume(x, &volume[j], n_nodes);
        compute_cell_center(x, center[j], n_nodes);
    }
    mesh->cell.volume = volume;
    mesh->cell.center = center;
}

static void compute_cell_volume(const double (*x)[N_DIMS], double *volume, long n_nodes)
{
    switch (n_nodes) {
        case 4: {
            double AB[N_DIMS], AC[N_DIMS], AD[N_DIMS], cross[N_DIMS];
            for (long d = 0; d < N_DIMS; ++d) {
                AB[d] = x[1][d] - x[0][d];
                AC[d] = x[2][d] - x[0][d];
                AD[d] = x[3][d] - x[0][d];
            }
            *volume = fabs(array_dot(AB, array_cross(AC, AD, cross), N_DIMS)) / 6;
            break;
        }
        case 5: {
            double x0[4][N_DIMS], x1[4][N_DIMS], v0, v1;
            for (long d = 0; d < N_DIMS; ++d) {
                x0[0][d] = x[0][d];
                x0[1][d] = x[1][d];
                x0[2][d] = x[2][d];
                x0[3][d] = x[4][d];
                x1[0][d] = x[0][d];
                x1[1][d] = x[2][d];
                x1[2][d] = x[3][d];
                x1[3][d] = x[4][d];
            }
            compute_cell_volume(x0, &v0, 4);
            compute_cell_volume(x1, &v1, 4);
            *volume = v0 + v1;
            break;
        }
        case 6: {
            double x0[5][N_DIMS], x1[4][N_DIMS], v0, v1;
            for (long d = 0; d < N_DIMS; ++d) {
                x0[0][d] = x[0][d];
                x0[1][d] = x[1][d];
                x0[2][d] = x[4][d];
                x0[3][d] = x[3][d];
                x0[4][d] = x[5][d];
                x1[0][d] = x[0][d];
                x1[1][d] = x[1][d];
                x1[2][d] = x[2][d];
                x1[3][d] = x[5][d];
            }
            compute_cell_volume(x0, &v0, 5);
            compute_cell_volume(x1, &v1, 4);
            *volume = v0 + v1;
            break;
        }
        case 8: {
            double x0[6][N_DIMS], x1[6][N_DIMS], v0, v1;
            for (long d = 0; d < N_DIMS; ++d) {
                x0[0][d] = x[0][d];
                x0[1][d] = x[1][d];
                x0[2][d] = x[2][d];
                x0[3][d] = x[4][d];
                x0[4][d] = x[5][d];
                x0[5][d] = x[6][d];
                x1[0][d] = x[0][d];
                x1[1][d] = x[2][d];
                x1[2][d] = x[3][d];
                x1[3][d] = x[4][d];
                x1[4][d] = x[6][d];
                x1[5][d] = x[7][d];
            }
            compute_cell_volume(x0, &v0, 6);
            compute_cell_volume(x1, &v1, 6);
            *volume = v0 + v1;
            break;
        }
        default: error("unsupported number of cell nodes '%ld'", n_nodes);
    }
    ensure(*volume > 0);
}

static void compute_cell_center(const double (*x)[N_DIMS], double *center, long n_nodes)
{
    switch (n_nodes) {
        case 4: {
            for (long d = 0; d < N_DIMS; ++d) center[d] = 0;
            for (long i = 0; i < 4; ++i)
                for (long d = 0; d < N_DIMS; ++d) center[d] += x[i][d] / 4;
            break;
        }
        case 5: {
            double x0[4][N_DIMS], x1[4][N_DIMS], v0, v1, c0[N_DIMS], c1[N_DIMS];
            for (long d = 0; d < N_DIMS; ++d) {
                x0[0][d] = x[0][d];
                x0[1][d] = x[1][d];
                x0[2][d] = x[2][d];
                x0[3][d] = x[4][d];
                x1[0][d] = x[0][d];
                x1[1][d] = x[2][d];
                x1[2][d] = x[3][d];
                x1[3][d] = x[4][d];
            }
            compute_cell_volume(x0, &v0, 4);
            compute_cell_volume(x1, &v1, 4);
            compute_cell_center(x0, c0, 4);
            compute_cell_center(x1, c1, 4);
            for (long d = 0; d < N_DIMS; ++d) center[d] = (v0 * c0[d] + v1 * c1[d]) / (v0 + v1);
            break;
        }
        case 6: {
            double x0[5][N_DIMS], x1[4][N_DIMS], v0, v1, c0[N_DIMS], c1[N_DIMS];
            for (long d = 0; d < N_DIMS; ++d) {
                x0[0][d] = x[0][d];
                x0[1][d] = x[1][d];
                x0[2][d] = x[4][d];
                x0[3][d] = x[3][d];
                x0[4][d] = x[5][d];
                x1[0][d] = x[0][d];
                x1[1][d] = x[1][d];
                x1[2][d] = x[2][d];
                x1[3][d] = x[5][d];
            }
            compute_cell_volume(x0, &v0, 5);
            compute_cell_volume(x1, &v1, 4);
            compute_cell_center(x0, c0, 5);
            compute_cell_center(x1, c1, 4);
            for (long d = 0; d < N_DIMS; ++d) center[d] = (v0 * c0[d] + v1 * c1[d]) / (v0 + v1);
            break;
        }
        case 8: {
            double x0[6][N_DIMS], x1[6][N_DIMS], v0, v1, c0[N_DIMS], c1[N_DIMS];
            for (long d = 0; d < N_DIMS; ++d) {
                x0[0][d] = x[0][d];
                x0[1][d] = x[1][d];
                x0[2][d] = x[2][d];
                x0[3][d] = x[4][d];
                x0[4][d] = x[5][d];
                x0[5][d] = x[6][d];
                x1[0][d] = x[0][d];
                x1[1][d] = x[2][d];
                x1[2][d] = x[3][d];
                x1[3][d] = x[4][d];
                x1[4][d] = x[6][d];
                x1[5][d] = x[7][d];
            }
            compute_cell_volume(x0, &v0, 6);
            compute_cell_volume(x1, &v1, 6);
            compute_cell_center(x0, c0, 6);
            compute_cell_center(x1, c1, 6);
            for (long d = 0; d < N_DIMS; ++d) center[d] = (v0 * c0[d] + v1 * c1[d]) / (v0 + v1);
            break;
        }
        default: error("unsupported number of cell nodes '%ld'", n_nodes);
    }
}

static void compute_cell_projections(Mesh *mesh)
{
    const ALIAS(cell, mesh->face.cell);
    const ALIAS(area, mesh->face.area);
    const ALIAS(normal, mesh->face.normal);
    double(*projection)[N_DIMS] = memory_calloc(mesh->n_inner_cells, sizeof(*projection));
    for (long i = 0; i < mesh->n_faces; ++i)
        for (long s = 0; s < N_SIDES && cell[i][s] < mesh->n_inner_cells; ++s)
            for (long d = 0; d < N_DIMS; ++d)
                projection[cell[i][s]][d] += area[i] * fabs(normal[i][d]) / 2;
    mesh->cell.projection = projection;
}

static void compute_ghost_cell_centers(const Mesh *mesh)
{
    const ALIAS(j_cell, mesh->entity.j_cell);
    const ALIAS(i_cell, mesh->cell.i_cell);
    const ALIAS(cell, mesh->cell.cell);
    const ALIAS(i_node, mesh->cell.i_node);
    const ALIAS(coord, mesh->node.coord);
    const ALIAS(offset, mesh->entity.offset);
    ALIAS(node, mesh->cell.node);
    ALIAS(center, mesh->cell.center);
    for (long e = 0; e < mesh->n_entities; ++e) {
        if (j_cell[e] < mesh->n_inner_cells) continue;
        for (long j = j_cell[e]; j < j_cell[e + 1]; ++j) {
            const long n_cells = i_cell[j + 1] - i_cell[j];
            switch (n_cells) {
                case 1: {  // ghost cell: mirror inner cell center at face
                    const long n_nodes = i_node[j + 1] - i_node[j];
                    correct_face_node_order(mesh, &node[i_node[j]], n_nodes);

                    double x[MAX_FACE_NODES][N_DIMS] = {0};
                    for (long n = 0, i = i_node[j]; i < i_node[j + 1]; ++i, ++n)
                        for (long d = 0; d < N_DIMS; ++d) x[n][d] = coord[node[i]][d];

                    double cf[N_DIMS], n[N_DIMS];
                    compute_face_center(x, cf, n_nodes);
                    compute_face_normal(x, n, n_nodes);

                    const long inner = cell[i_cell[j]];
                    const double *ci = center[inner];

                    double CF[N_DIMS];
                    for (long d = 0; d < N_DIMS; ++d) CF[d] = cf[d] - ci[d];
                    const double dot = array_dot(CF, n, N_DIMS);

                    for (long d = 0; d < N_DIMS; ++d) center[j][d] = ci[d] + 2 * dot * n[d];
                    break;
                }
                case 2: {  // periodic cell: offset outer cell center
                    const long outer = cell[i_cell[j] + 1];
                    const double *co = center[outer];

                    for (long d = 0; d < N_DIMS; ++d) center[j][d] = co[d] + offset[e][d];
                    break;
                }
                default: error("unsupported number of ghost cell connections '%ld'", n_cells);
            }
        }
    }
}

static void find_periodic_connections(const Mesh *mesh, Dict *periodic)
{
    const ALIAS(i_cell, mesh->cell.i_cell);
    for (long j = mesh->n_inner_cells; j < mesh->n_inner_cells + mesh->n_ghost_cells; ++j) {
        if (i_cell[j + 1] - i_cell[j] != N_SIDES) continue;
        const long *cell = &mesh->cell.cell[i_cell[j]];
        dict_insert(periodic, cell, N_SIDES, &j, 1);
    }
}

static void compute_face_reconstruction(Mesh *mesh, const Dict *periodic)
{
    const ALIAS(cell, mesh->face.cell);
    const ALIAS(center, mesh->cell.center);
    cleanup double(*dx)[N_DIMS] = memory_calloc(mesh->n_faces, sizeof(*dx));
    cleanup double *theta = memory_calloc(mesh->n_faces, sizeof(*theta));
    for (long i = 0; i < mesh->n_faces; ++i) {
        const long inner = cell[i][L];
        long outer = cell[i][R];

        long *c;
        if (dict_lookup(periodic, cell[i], N_SIDES, &c)) outer = *c;

        for (long d = 0; d < N_DIMS; ++d) dx[i][d] = center[outer][d] - center[inner][d];
        theta[i] = 1 / array_norm(dx[i], N_DIMS);
    }

    cleanup double *r11 = memory_calloc(mesh->n_inner_cells, sizeof(*r11));
    cleanup double *r12 = memory_calloc(mesh->n_inner_cells, sizeof(*r12));
    cleanup double *r22 = memory_calloc(mesh->n_inner_cells, sizeof(*r22));
    cleanup double *r13 = memory_calloc(mesh->n_inner_cells, sizeof(*r13));
    cleanup double *r23 = memory_calloc(mesh->n_inner_cells, sizeof(*r23));
    cleanup double *r33 = memory_calloc(mesh->n_inner_cells, sizeof(*r33));
    for (long i = 0; i < mesh->n_faces; ++i) {
        for (long s = 0; s < N_SIDES && cell[i][s] < mesh->n_inner_cells; ++s) {
            r11[cell[i][s]] += sq(theta[i] * dx[i][X]);
            r12[cell[i][s]] += sq(theta[i]) * dx[i][X] * dx[i][Y];
            r22[cell[i][s]] += sq(theta[i] * dx[i][Y]);
            r13[cell[i][s]] += sq(theta[i]) * dx[i][X] * dx[i][Z];
            r23[cell[i][s]] += sq(theta[i]) * dx[i][Y] * dx[i][Z];
            r33[cell[i][s]] += sq(theta[i] * dx[i][Z]);
        }
    }
    for (long i = 0; i < mesh->n_inner_cells; ++i) {
        r11[i] = sqrt(r11[i]) + EPS;
        r12[i] = r12[i] / r11[i];
        r22[i] = sqrt(r22[i] - sq(r12[i])) + EPS;
        r13[i] = r13[i] / r11[i];
        r23[i] = (r23[i] - r12[i] * r13[i]) / r22[i];
        r33[i] = sqrt(r33[i] - (sq(r13[i]) + sq(r23[i]))) + EPS;
    }

    double(*weight)[N_DIMS] = memory_calloc(mesh->n_faces, sizeof(*weight));
    double(*recon)[N_SIDES][N_DIMS] = memory_calloc(mesh->n_faces, sizeof(*recon));
    double(*correction)[N_DIMS + 1] = memory_calloc(mesh->n_faces, sizeof(*correction));
    for (long i = 0; i < mesh->n_faces; ++i) {
        const long cL = cell[i][L];
        const double a1 = dx[i][X] / sq(r11[cL]);
        const double a2 = (dx[i][Y] - r12[cL] / r11[cL] * dx[i][X]) / sq(r22[cL]);
        const double b = (r12[cL] * r23[cL] - r13[cL] * r22[cL]) / (r11[cL] * r22[cL]);
        const double a3 = (dx[i][Z] - r23[cL] / r22[cL] * dx[i][Y] + b * dx[i][X]) / sq(r33[cL]);
        weight[i][X] = sq(theta[i]) * (a1 - r12[cL] / r11[cL] * a2 + b * a3);
        weight[i][Y] = sq(theta[i]) * (a2 - r23[cL] / r22[cL] * a3);
        weight[i][Z] = sq(theta[i]) * (a3);
        ensure(isfinite(weight[i][X]) && isfinite(weight[i][Y]) && isfinite(weight[i][Z]));

        const double *cf = mesh->face.center[i];
        const double *ci = center[cL];
        for (long d = 0; d < N_DIMS; ++d) {
            recon[i][L][d] = cf[d] - ci[d];
            recon[i][R][d] = cf[d] - (ci[d] + dx[i][d]);
        }

        const double l = array_norm(dx[i], N_DIMS);
        for (long d = 0; d < N_DIMS; ++d) correction[i][d] = dx[i][d] / l;
        correction[i][N_DIMS] = l;
    }
    mesh->face.gradient_weight = weight;
    mesh->face.reconstruction = recon;
    mesh->face.gradient_correction = correction;
}

static void compute_cell_reconstruction(Mesh *mesh, const Dict *periodic)
{
    const ALIAS(center, mesh->cell.center);
    const ALIAS(i_cell, mesh->cell.i_cell);
    const ALIAS(cell, mesh->cell.cell);
    const ALIAS(coord, mesh->node.coord);
    double(*recon)[N_DIMS] = memory_calloc(i_cell[mesh->n_inner_cells], sizeof(*recon));
    for (long j = 0; j < mesh->n_inner_cells; ++j) {
        for (long i = i_cell[j]; i < i_cell[j + 1]; ++i) {
            long outer = cell[i];

            long *c;
            if (dict_lookup(periodic, (long[]){j, outer}, N_SIDES, &c)) outer = *c;

            long node[MAX_FACE_NODES] = {0};
            const long n_nodes = connectivity_nodes(mesh, node, j, outer);
            ensure(n_nodes >= N_DIMS);
            correct_face_node_order(mesh, node, n_nodes);

            double xf[MAX_FACE_NODES][N_DIMS] = {0};
            for (long n = 0; n < n_nodes; ++n)
                for (long d = 0; d < N_DIMS; ++d) xf[n][d] = coord[node[n]][d];

            double cf[N_DIMS];
            compute_face_center(xf, cf, n_nodes);

            for (long d = 0; d < N_DIMS; ++d) recon[i][d] = cf[d] - center[j][d];
        }
    }
    mesh->cell.reconstruction = recon;
}
