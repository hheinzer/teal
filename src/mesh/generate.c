#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "connectivity.h"
#include "mesh.h"
#include "partition.h"
#include "periodic.h"
#include "sync.h"
#include "teal/array.h"
#include "teal/dict.h"
#include "teal/isclose.h"
#include "teal/memory.h"
#include "teal/sync.h"
#include "teal/utils.h"

static void compute_face_geometry(Mesh *mesh);

static void compute_cell_geometry(Mesh *mesh);

static void compute_cell_projections(Mesh *mesh);

static void compute_ghost_cell_centers(Mesh *mesh);

static void compute_face_reconstruction(Mesh *mesh, const Dict *periodic);

static void compute_cell_reconstruction(Mesh *mesh, const Dict *periodic);

static void correct_face_node_order(const Mesh *mesh, long *node, long n_nodes);

static double compute_face_area(const Vector3d *x, long n_nodes);

static void compute_face_center(Vector3d center, const Vector3d *x, long n_nodes);

static void compute_face_normal(Vector3d normal, const Vector3d *x, long n_nodes);

static void correct_face_normal_direction(const Mesh *mesh, const Vector2l cell,
                                          const Vector3d center, Vector3d normal);

static void compute_face_basis(Matrix3d basis);

static double compute_cell_volume(const Vector3d *x, long n_nodes);

static void compute_cell_center(Vector3d center, const Vector3d *x, long n_nodes);

void mesh_generate(Mesh **mesh)
{
    partition(mesh);
    connectivity_faces(*mesh);

    compute_face_geometry(*mesh);
    compute_cell_geometry(*mesh);
    compute_cell_projections(*mesh);

    sync_all(*mesh, *(*mesh)->cell.center, N_DIMS);
    compute_ghost_cell_centers(*mesh);

    defer(dict_free) Dict *periodic = dict_create((*mesh)->n_ghost_cells);
    periodic_cell_to_cell(*mesh, periodic);
    compute_face_reconstruction(*mesh, periodic);
    compute_cell_reconstruction(*mesh, periodic);

    (*mesh)->volume = sync_sum(array_sum((*mesh)->cell.volume, (*mesh)->n_inner_cells));
}

static void compute_face_geometry(Mesh *mesh)
{
    const alias(coord, mesh->node.coord);
    const alias(i_node, mesh->face.i_node);
    const alias(cell, mesh->face.cell);
    alias(node, mesh->face.node);
    double *area = memory_calloc(mesh->n_faces, sizeof(*area));
    Vector3d *center = memory_calloc(mesh->n_faces, sizeof(*center));
    Matrix3d *basis = memory_calloc(mesh->n_faces, sizeof(*basis));
    for (long j = 0; j < mesh->n_faces; ++j) {
        const long n_nodes = i_node[j + 1] - i_node[j];
        correct_face_node_order(mesh, &node[i_node[j]], n_nodes);

        Vector3d x[MAX_FACE_NODES] = {0};
        for (long n = 0, i = i_node[j]; i < i_node[j + 1]; ++i, ++n)
            for (long d = 0; d < N_DIMS; ++d) x[n][d] = coord[node[i]][d];

        area[j] = compute_face_area(x, n_nodes);
        compute_face_center(center[j], x, n_nodes);
        compute_face_normal(basis[j][X], x, n_nodes);
        correct_face_normal_direction(mesh, cell[j], center[j], basis[j][X]);
        compute_face_basis(basis[j]);
    }
    mesh->face.area = area;
    mesh->face.center = center;
    mesh->face.basis = basis;
}

static void compute_cell_geometry(Mesh *mesh)
{
    const alias(coord, mesh->node.coord);
    const alias(i_node, mesh->cell.i_node);
    const alias(node, mesh->cell.node);
    double *volume = memory_calloc(mesh->n_inner_cells, sizeof(*volume));
    Vector3d *center = memory_calloc(mesh->n_cells, sizeof(*center));
    for (long j = 0; j < mesh->n_inner_cells; ++j) {
        const long n_nodes = i_node[j + 1] - i_node[j];

        Vector3d x[MAX_CELL_NODES] = {0};
        for (long n = 0, i = i_node[j]; i < i_node[j + 1]; ++i, ++n)
            for (long d = 0; d < N_DIMS; ++d) x[n][d] = coord[node[i]][d];

        volume[j] = compute_cell_volume(x, n_nodes);
        compute_cell_center(center[j], x, n_nodes);
    }
    mesh->cell.volume = volume;
    mesh->cell.center = center;
}

static void compute_cell_projections(Mesh *mesh)
{
    const alias(cell, mesh->face.cell);
    const alias(area, mesh->face.area);
    const alias(basis, mesh->face.basis);
    Vector3d *projection = memory_calloc(mesh->n_inner_cells, sizeof(*projection));
    for (long i = 0; i < mesh->n_faces; ++i)
        for (long s = 0; s < N_SIDES && cell[i][s] < mesh->n_inner_cells; ++s)
            for (long d = 0; d < N_DIMS; ++d)
                projection[cell[i][s]][d] += area[i] * fabs(basis[i][X][d]) / 2;
    mesh->cell.projection = projection;
}

static void compute_ghost_cell_centers(Mesh *mesh)
{
    const alias(coord, mesh->node.coord);
    const alias(i_node, mesh->cell.i_node);
    const alias(i_cell, mesh->cell.i_cell);
    const alias(cell, mesh->cell.cell);
    const alias(j_cell, mesh->entity.j_cell);
    const alias(offset, mesh->entity.offset);
    alias(node, mesh->cell.node);
    alias(center, mesh->cell.center);
    for (long e = 0; e < mesh->n_entities; ++e) {
        if (j_cell[e] < mesh->n_inner_cells) continue;
        for (long j = j_cell[e]; j < j_cell[e + 1]; ++j) {
            const long n_cells = i_cell[j + 1] - i_cell[j];
            switch (n_cells) {
                case 1: {  // ghost cell: mirror inner cell center at face
                    const long n_nodes = i_node[j + 1] - i_node[j];
                    correct_face_node_order(mesh, &node[i_node[j]], n_nodes);

                    Vector3d x[MAX_FACE_NODES] = {0};
                    for (long n = 0, i = i_node[j]; i < i_node[j + 1]; ++i, ++n)
                        for (long d = 0; d < N_DIMS; ++d) x[n][d] = coord[node[i]][d];

                    Vector3d f, n;
                    compute_face_center(f, x, n_nodes);
                    compute_face_normal(n, x, n_nodes);

                    const long inner = cell[i_cell[j]];

                    Vector3d cf;
                    for (long d = 0; d < N_DIMS; ++d) cf[d] = f[d] - center[inner][d];
                    const double dot = array_dot(cf, n, N_DIMS);

                    for (long d = 0; d < N_DIMS; ++d)
                        center[j][d] = center[inner][d] + 2 * dot * n[d];
                    break;
                }
                case 2: {  // periodic cell: offset outer cell center
                    const long outer = cell[i_cell[j] + 1];

                    for (long d = 0; d < N_DIMS; ++d)
                        center[j][d] = center[outer][d] + offset[e][d];
                    break;
                }
                default: abort();
            }
        }
    }
}

static void compute_face_reconstruction(Mesh *mesh, const Dict *periodic)
{
    const alias(center, mesh->cell.center);
    const alias(cell, mesh->face.cell);
    smart Vector3d *dx = memory_calloc(mesh->n_faces, sizeof(*dx));
    smart double *theta = memory_calloc(mesh->n_faces, sizeof(*theta));
    for (long i = 0; i < mesh->n_faces; ++i) {
        const long inner = cell[i][L];
        long outer = cell[i][R];

        DictItem *item = dict_lookup(periodic, cell[i], N_SIDES);
        if (item) outer = *item->val;

        for (long d = 0; d < N_DIMS; ++d) dx[i][d] = center[outer][d] - center[inner][d];
        theta[i] = 1 / array_norm(dx[i], N_DIMS);
    }

    smart double *r11 = memory_calloc(mesh->n_inner_cells, sizeof(*r11));
    smart double *r12 = memory_calloc(mesh->n_inner_cells, sizeof(*r12));
    smart double *r22 = memory_calloc(mesh->n_inner_cells, sizeof(*r22));
    smart double *r13 = memory_calloc(mesh->n_inner_cells, sizeof(*r13));
    smart double *r23 = memory_calloc(mesh->n_inner_cells, sizeof(*r23));
    smart double *r33 = memory_calloc(mesh->n_inner_cells, sizeof(*r33));
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
        r11[i] = sqrt(r11[i]);
        r12[i] = r12[i] / r11[i];
        r22[i] = sqrt(r22[i] - sq(r12[i]));
        r13[i] = r13[i] / r11[i];
        r23[i] = (r23[i] - r12[i] * r13[i]) / r22[i];
        r33[i] = sqrt(r33[i] - (sq(r13[i]) + sq(r23[i])));
    }

    Vector3d(*to_cell)[N_SIDES] = memory_calloc(mesh->n_faces, sizeof(*to_cell));
    Vector3d *weight = memory_calloc(mesh->n_faces, sizeof(*weight));
    Vector4d *correction = memory_calloc(mesh->n_faces, sizeof(*correction));
    for (long i = 0; i < mesh->n_faces; ++i) {
        const long cL = cell[i][L];
        const double *cf = mesh->face.center[i];
        const double *ci = center[cL];
        for (long d = 0; d < N_DIMS; ++d) {
            to_cell[i][L][d] = cf[d] - ci[d];
            to_cell[i][R][d] = cf[d] - (ci[d] + dx[i][d]);
        }

        const double a1 = dx[i][X] / sq(r11[cL]);
        const double a2 = (dx[i][Y] - r12[cL] / r11[cL] * dx[i][X]) / sq(r22[cL]);
        const double b = (r12[cL] * r23[cL] - r13[cL] * r22[cL]) / (r11[cL] * r22[cL]);
        const double a3 = (dx[i][Z] - r23[cL] / r22[cL] * dx[i][Y] + b * dx[i][X]) / sq(r33[cL]);
        if (!is_close(dx[i][X], 0)) {
            weight[i][X] += a1;
        }
        if (!is_close(dx[i][Y], 0)) {
            weight[i][X] += -r12[cL] / r11[cL] * a2;
            weight[i][Y] += a2;
        }
        if (!is_close(dx[i][Z], 0)) {
            weight[i][X] += b * a3;
            weight[i][Y] += -r23[cL] / r22[cL] * a3;
            weight[i][Z] += a3;
        }
        weight[i][X] *= sq(theta[i]);
        weight[i][Y] *= sq(theta[i]);
        weight[i][Z] *= sq(theta[i]);
        assert(isfinite(weight[i][X]) && isfinite(weight[i][Y]) && isfinite(weight[i][Z]));

        const double l = array_norm(dx[i], N_DIMS);
        for (long d = 0; d < N_DIMS; ++d) correction[i][d] = dx[i][d] / l;
        correction[i][N_DIMS] = l;
    }
    mesh->face.weight = weight;
    mesh->face.to_cell = to_cell;
    mesh->face.correction = correction;
}

static void compute_cell_reconstruction(Mesh *mesh, const Dict *periodic)
{
    const alias(coord, mesh->node.coord);
    const alias(i_cell, mesh->cell.i_cell);
    const alias(cell, mesh->cell.cell);
    const alias(center, mesh->cell.center);
    Vector3d *to_cell = memory_calloc(i_cell[mesh->n_inner_cells], sizeof(*to_cell));
    for (long j = 0; j < mesh->n_inner_cells; ++j) {
        for (long i = i_cell[j]; i < i_cell[j + 1]; ++i) {
            long outer = cell[i];

            DictItem *item = dict_lookup(periodic, (Vector2l){j, outer}, N_SIDES);
            if (item) outer = *item->val;

            long node[MAX_FACE_NODES];
            const long n_nodes = connectivity_nodes(mesh, node, j, outer);
            assert(n_nodes >= N_DIMS);
            correct_face_node_order(mesh, node, n_nodes);

            Vector3d xf[MAX_FACE_NODES];
            for (long n = 0; n < n_nodes; ++n)
                for (long d = 0; d < N_DIMS; ++d) xf[n][d] = coord[node[n]][d];

            Vector3d cf;
            compute_face_center(cf, xf, n_nodes);

            for (long d = 0; d < N_DIMS; ++d) to_cell[i][d] = cf[d] - center[j][d];
        }
    }
    mesh->cell.to_cell = to_cell;
}

static void correct_face_node_order(const Mesh *mesh, long *node, long n_nodes)
{
    const alias(x, mesh->node.coord);
    switch (n_nodes) {
        case 3: break;  // node order of triangle does not need to be corrected
        case 4: {
            Vector3d ab, ac, ad, n0, n1;
            for (long d = 0; d < N_DIMS; ++d) {
                ab[d] = x[node[1]][d] - x[node[0]][d];
                ac[d] = x[node[2]][d] - x[node[0]][d];
                ad[d] = x[node[3]][d] - x[node[0]][d];
            }
            array_cross(ab, ac, n0);
            array_cross(ac, ad, n1);
            if (array_dot(n0, n1, N_DIMS) > 0) return;

            const long swap1 = node[1];
            node[1] = node[2];
            node[2] = swap1;
            for (long d = 0; d < N_DIMS; ++d) {
                ab[d] = x[node[1]][d] - x[node[0]][d];
                ac[d] = x[node[2]][d] - x[node[0]][d];
                ad[d] = x[node[3]][d] - x[node[0]][d];
            }
            array_cross(ab, ac, n0);
            array_cross(ac, ad, n1);
            if (array_dot(n0, n1, N_DIMS) > 0) return;

            const long swap2 = node[2];
            node[2] = node[3];
            node[3] = swap2;
            for (long d = 0; d < N_DIMS; ++d) {
                ab[d] = x[node[1]][d] - x[node[0]][d];
                ac[d] = x[node[2]][d] - x[node[0]][d];
                ad[d] = x[node[3]][d] - x[node[0]][d];
            }
            array_cross(ab, ac, n0);
            array_cross(ac, ad, n1);
            assert(array_dot(n0, n1, N_DIMS) > 0);
            break;
        }
        default: abort();
    }
}

static double compute_face_area(const Vector3d *x, long n_nodes)
{
    switch (n_nodes) {
        case 3: {
            Vector3d ab, ac, cross;
            for (long d = 0; d < N_DIMS; ++d) {
                ab[d] = x[1][d] - x[0][d];
                ac[d] = x[2][d] - x[0][d];
            }
            const double area = array_norm(array_cross(ab, ac, cross), N_DIMS) / 2;
            assert(area > 0);
            return area;
        }
        case 4: {
            Vector3d xa[3], xb[3];
            for (long d = 0; d < N_DIMS; ++d) {
                xa[0][d] = x[0][d];
                xa[1][d] = x[1][d];
                xa[2][d] = x[2][d];
                xb[0][d] = x[0][d];
                xb[1][d] = x[2][d];
                xb[2][d] = x[3][d];
            }
            return compute_face_area(xa, 3) + compute_face_area(xb, 3);
        }
        default: abort();
    }
}

static void compute_face_center(Vector3d center, const Vector3d *x, long n_nodes)
{
    switch (n_nodes) {
        case 3: {
            for (long d = 0; d < N_DIMS; ++d) center[d] = 0;
            for (long i = 0; i < 3; ++i)
                for (long d = 0; d < N_DIMS; ++d) center[d] += x[i][d] / 3;
            break;
        }
        case 4: {
            Vector3d xa[3], xb[3], ca, cb;
            for (long d = 0; d < N_DIMS; ++d) {
                xa[0][d] = x[0][d];
                xa[1][d] = x[1][d];
                xa[2][d] = x[2][d];
                xb[0][d] = x[0][d];
                xb[1][d] = x[2][d];
                xb[2][d] = x[3][d];
            }
            const double aa = compute_face_area(xa, 3);
            const double ab = compute_face_area(xb, 3);
            compute_face_center(ca, xa, 3);
            compute_face_center(cb, xb, 3);
            for (long d = 0; d < N_DIMS; ++d) center[d] = (aa * ca[d] + ab * cb[d]) / (aa + ab);
            break;
        }
        default: abort();
    }
}

static void compute_face_normal(Vector3d normal, const Vector3d *x, long n_nodes)
{
    switch (n_nodes) {
        case 3: {
            Vector3d ab, ac;
            for (long d = 0; d < N_DIMS; ++d) {
                ab[d] = x[1][d] - x[0][d];
                ac[d] = x[2][d] - x[0][d];
            }
            array_normalize(array_cross(ab, ac, normal), N_DIMS);
            assert(is_close(array_norm(normal, N_DIMS), 1));
            break;
        }
        case 4: {
            Vector3d c;
            compute_face_center(c, x, n_nodes);
            for (long d = 0; d < N_DIMS; ++d) normal[d] = 0;
            for (long i = 0; i < 4; ++i) {
                Vector3d xi[3], ni;
                for (long d = 0; d < N_DIMS; ++d) {
                    xi[0][d] = c[d];
                    xi[1][d] = x[i][d];
                    xi[2][d] = x[(i + 1) % 4][d];
                }
                compute_face_normal(ni, xi, 3);
                for (long d = 0; d < N_DIMS; ++d) normal[d] += ni[d] / 4;
            }
            array_normalize(normal, N_DIMS);
            break;
        }
        default: abort();
    }
}

static void correct_face_normal_direction(const Mesh *mesh, const Vector2l cell,
                                          const Vector3d center, Vector3d normal)
{
    const alias(coord, mesh->node.coord);
    const alias(i_node, mesh->cell.i_node);
    const alias(node, mesh->cell.node);
    const long n_nodes = i_node[cell[L] + 1] - i_node[cell[L]];
    Vector3d cmean = {0};
    for (long i = i_node[cell[L]]; i < i_node[cell[L] + 1]; ++i)
        for (long d = 0; d < N_DIMS; ++d) cmean[d] += coord[node[i]][d] / n_nodes;

    Vector3d fc;
    for (long d = 0; d < N_DIMS; ++d) fc[d] = cmean[d] - center[d];

    const double dot = array_dot(fc, normal, N_DIMS);
    assert(!is_close(dot, 0));
    if (dot > 0)  // make sure that normal points outwards wrt. inner cell
        for (long d = 0; d < N_DIMS; ++d) normal[d] *= -1;
}

static void compute_face_basis(Matrix3d basis)
{
    // Billett and Toro 1998
    const double nx = basis[X][X], ny = basis[X][Y], nz = basis[X][Z];
    const double nqz = sqrt(1 - sq(nz));
    const double nqy = sqrt(1 - sq(ny));
    if (nqz > nqy) {
        basis[Y][X] = -ny / nqz;
        basis[Y][Y] = nx / nqz;
        basis[Y][Z] = 0;
        basis[Z][X] = -nx * nz / nqz;
        basis[Z][Y] = -ny * nz / nqz;
        basis[Z][Z] = nqz;
    }
    else {
        basis[Y][X] = -nx * ny / nqy;
        basis[Y][Y] = nqy;
        basis[Y][Z] = -ny * nz / nqy;
        basis[Z][X] = -nz / nqy;
        basis[Z][Y] = 0;
        basis[Z][Z] = nx / nqy;
    }
    for (long d = 0; d < N_DIMS; ++d) assert(is_close(array_norm(basis[d], N_DIMS), 1));
}

static double compute_cell_volume(const Vector3d *x, long n_nodes)
{
    switch (n_nodes) {
        case 4: {
            Vector3d ab, ac, ad, cross;
            for (long d = 0; d < N_DIMS; ++d) {
                ab[d] = x[1][d] - x[0][d];
                ac[d] = x[2][d] - x[0][d];
                ad[d] = x[3][d] - x[0][d];
            }
            const double volume = fabs(array_dot(ab, array_cross(ac, ad, cross), N_DIMS)) / 6;
            assert(volume > 0);
            return volume;
        }
        case 5: {
            Vector3d xa[4], xb[4];
            for (long d = 0; d < N_DIMS; ++d) {
                xa[0][d] = x[0][d];
                xa[1][d] = x[1][d];
                xa[2][d] = x[2][d];
                xa[3][d] = x[4][d];
                xb[0][d] = x[0][d];
                xb[1][d] = x[2][d];
                xb[2][d] = x[3][d];
                xb[3][d] = x[4][d];
            }
            return compute_cell_volume(xa, 4) + compute_cell_volume(xb, 4);
        }
        case 6: {
            Vector3d xa[5], xb[4];
            for (long d = 0; d < N_DIMS; ++d) {
                xa[0][d] = x[0][d];
                xa[1][d] = x[1][d];
                xa[2][d] = x[4][d];
                xa[3][d] = x[3][d];
                xa[4][d] = x[5][d];
                xb[0][d] = x[0][d];
                xb[1][d] = x[1][d];
                xb[2][d] = x[2][d];
                xb[3][d] = x[5][d];
            }
            return compute_cell_volume(xa, 5) + compute_cell_volume(xb, 4);
        }
        case 8: {
            Vector3d xa[6], xb[6];
            for (long d = 0; d < N_DIMS; ++d) {
                xa[0][d] = x[0][d];
                xa[1][d] = x[1][d];
                xa[2][d] = x[2][d];
                xa[3][d] = x[4][d];
                xa[4][d] = x[5][d];
                xa[5][d] = x[6][d];
                xb[0][d] = x[0][d];
                xb[1][d] = x[2][d];
                xb[2][d] = x[3][d];
                xb[3][d] = x[4][d];
                xb[4][d] = x[6][d];
                xb[5][d] = x[7][d];
            }
            return compute_cell_volume(xa, 6) + compute_cell_volume(xb, 6);
        }
        default: abort();
    }
}

static void compute_cell_center(Vector3d center, const Vector3d *x, long n_nodes)
{
    switch (n_nodes) {
        case 4: {
            for (long d = 0; d < N_DIMS; ++d) center[d] = 0;
            for (long i = 0; i < 4; ++i)
                for (long d = 0; d < N_DIMS; ++d) center[d] += x[i][d] / 4;
            break;
        }
        case 5: {
            Vector3d xa[4], xb[4], ca, cb;
            for (long d = 0; d < N_DIMS; ++d) {
                xa[0][d] = x[0][d];
                xa[1][d] = x[1][d];
                xa[2][d] = x[2][d];
                xa[3][d] = x[4][d];
                xb[0][d] = x[0][d];
                xb[1][d] = x[2][d];
                xb[2][d] = x[3][d];
                xb[3][d] = x[4][d];
            }
            const double va = compute_cell_volume(xa, 4);
            const double vb = compute_cell_volume(xb, 4);
            compute_cell_center(ca, xa, 4);
            compute_cell_center(cb, xb, 4);
            for (long d = 0; d < N_DIMS; ++d) center[d] = (va * ca[d] + vb * cb[d]) / (va + vb);
            break;
        }
        case 6: {
            Vector3d xa[5], xb[4], ca, cb;
            for (long d = 0; d < N_DIMS; ++d) {
                xa[0][d] = x[0][d];
                xa[1][d] = x[1][d];
                xa[2][d] = x[4][d];
                xa[3][d] = x[3][d];
                xa[4][d] = x[5][d];
                xb[0][d] = x[0][d];
                xb[1][d] = x[1][d];
                xb[2][d] = x[2][d];
                xb[3][d] = x[5][d];
            }
            const double va = compute_cell_volume(xa, 5);
            const double vb = compute_cell_volume(xb, 4);
            compute_cell_center(ca, xa, 5);
            compute_cell_center(cb, xb, 4);
            for (long d = 0; d < N_DIMS; ++d) center[d] = (va * ca[d] + vb * cb[d]) / (va + vb);
            break;
        }
        case 8: {
            Vector3d xa[6], xb[6], ca, cb;
            for (long d = 0; d < N_DIMS; ++d) {
                xa[0][d] = x[0][d];
                xa[1][d] = x[1][d];
                xa[2][d] = x[2][d];
                xa[3][d] = x[4][d];
                xa[4][d] = x[5][d];
                xa[5][d] = x[6][d];
                xb[0][d] = x[0][d];
                xb[1][d] = x[2][d];
                xb[2][d] = x[3][d];
                xb[3][d] = x[4][d];
                xb[4][d] = x[6][d];
                xb[5][d] = x[7][d];
            }
            const double v0 = compute_cell_volume(xa, 6);
            const double v1 = compute_cell_volume(xb, 6);
            compute_cell_center(ca, xa, 6);
            compute_cell_center(cb, xb, 6);
            for (long d = 0; d < N_DIMS; ++d) center[d] = (v0 * ca[d] + v1 * cb[d]) / (v0 + v1);
            break;
        }
        default: abort();
    }
}
