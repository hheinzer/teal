#include <assert.h>
#include <string.h>

#include "private.h"
#include "sync2.h"
#include "teal2.h"
#include "utils2.h"

static int compute_dims(int *dims, Triple num_cells)
{
    long count[3] = {1, 1, 1};
    count[0] = (count[0] > num_cells.x) ? count[0] : num_cells.x;
    count[1] = (count[1] > num_cells.y) ? count[1] : num_cells.y;
    count[2] = (count[2] > num_cells.z) ? count[2] : num_cells.z;

    if (sync2.size > count[0] * count[1] * count[2]) {
        teal2_error("using more rank (%d) than cells (%ld) is not supported", sync2.size,
                    count[0] * count[1] * count[2]);
    }

    int perm[3] = {0, 1, 2};
    if (count[perm[1]] > count[perm[0]]) {
        int swap = perm[0];
        perm[0] = perm[1];
        perm[1] = swap;
    }
    if (count[perm[2]] > count[perm[1]]) {
        int swap = perm[1];
        perm[1] = perm[2];
        perm[2] = swap;
    }
    if (count[perm[1]] > count[perm[0]]) {
        int swap = perm[0];
        perm[0] = perm[1];
        perm[1] = swap;
    }

    dims[0] = dims[1] = dims[2] = 1;
    int ndims = (num_cells.x > 1) + (num_cells.y > 1) + (num_cells.z > 1);

    int prime = sync2.size;
    while (prime > 1) {
        int factor = 2;
        while (factor <= prime / factor && prime % factor != 0) {
            factor += 1;
        }

        if (factor > prime / factor) {
            factor = prime;
        }

        int best_axis = -1;
        long best_score = -1;
        for (int i = 0; i < ndims; i++) {
            int axis = perm[i];
            if (dims[axis] > count[axis] / factor) {
                continue;
            }
            long score = count[axis] / dims[axis];
            if (score > best_score) {
                best_axis = axis;
                best_score = score;
            }
        }
        assert(best_axis != -1);

        dims[best_axis] *= factor;
        prime /= factor;
    }

    return ndims;
}

static void split_bounds(double *min_coord, double *max_coord, int *num_cells, int dim,
                         const int *dims, const int *coords, int ndims)
{
    if (dim < ndims) {
        assert(!isclose(*min_coord, *max_coord) && *num_cells > 1);
        double width = (*max_coord - *min_coord) / *num_cells;
        int base = *num_cells / dims[dim];
        int extra = *num_cells % dims[dim];
        int offset = (coords[dim] * base) + (coords[dim] < extra ? coords[dim] : extra);
        *num_cells = base + (coords[dim] < extra);
        *min_coord = *min_coord + (offset * width);
        *max_coord = *min_coord + (*num_cells * width);
    }
    else {
        *num_cells = 1;
        if (isclose(*min_coord, *max_coord)) {
            *max_coord = *min_coord + 1;
        }
    }
}

static long *compute_globals(Mesh2 *mesh, Triple num_cells, const int *dims, const int *coords,
                             const int *neighbor)
{
    long *global = teal2_calloc(mesh->nodes.num, sizeof(*global));

    long prefix = mesh->nodes.num_inner;
    sync2_prefix(&prefix, 1, MPI_LONG);

    int num = 0;
    int num_inner = 0;
    for (int k = 0; k < num_cells.z + 1; k++) {
        for (int j = 0; j < num_cells.y + 1; j++) {
            for (int i = 0; i < num_cells.x + 1; i++) {
                if (!(i == 0 && coords[0]) && !(j == 0 && coords[1]) && !(k == 0 && coords[2])) {
                    global[num] = prefix + num_inner++;
                }
                if (i == 0 && coords[0]) {
                    MPI_Recv(&global[num], 1, MPI_LONG, neighbor[0], 0, sync2.comm,
                             MPI_STATUS_IGNORE);
                }
                if (j == 0 && coords[1]) {
                    MPI_Recv(&global[num], 1, MPI_LONG, neighbor[2], 1, sync2.comm,
                             MPI_STATUS_IGNORE);
                }
                if (k == 0 && coords[2]) {
                    MPI_Recv(&global[num], 1, MPI_LONG, neighbor[4], 2, sync2.comm,
                             MPI_STATUS_IGNORE);
                }
                if (i == num_cells.x && coords[0] < dims[0] - 1) {
                    MPI_Send(&global[num], 1, MPI_LONG, neighbor[1], 0, sync2.comm);
                }
                if (j == num_cells.y && coords[1] < dims[1] - 1) {
                    MPI_Send(&global[num], 1, MPI_LONG, neighbor[3], 1, sync2.comm);
                }
                if (k == num_cells.z && coords[2] < dims[2] - 1) {
                    MPI_Send(&global[num], 1, MPI_LONG, neighbor[5], 2, sync2.comm);
                }
                num += 1;
            }
        }
    }
    assert(num == mesh->nodes.num);
    assert(num_inner == mesh->nodes.num_inner);

    return global;
}

static void create_nodes(Mesh2 *mesh, Vector min_coord, Vector max_coord, Triple num_cells,
                         const int *dims, const int *coords, const int *neighbor)
{
    int num_nodes = (num_cells.x + 1) * (num_cells.y + 1) * (num_cells.z + 1);
    Vector *coord = teal2_calloc(num_nodes, sizeof(*coord));

    Vector del_coord = vector2_sub(max_coord, min_coord);
    del_coord.x /= num_cells.x;
    del_coord.y /= num_cells.y;
    del_coord.z /= num_cells.z;

    int num = 0;
    int num_inner = 0;
    for (int k = 0; k < num_cells.z + 1; k++) {
        for (int j = 0; j < num_cells.y + 1; j++) {
            for (int i = 0; i < num_cells.x + 1; i++) {
                if (!(i == 0 && coords[0]) && !(j == 0 && coords[1]) && !(k == 0 && coords[2])) {
                    num_inner += 1;
                }
                coord[num].x = min_coord.x + (i * del_coord.x);
                coord[num].y = min_coord.y + (j * del_coord.y);
                coord[num].z = min_coord.z + (k * del_coord.z);
                num += 1;
            }
        }
    }
    assert(num == num_nodes);

    mesh->nodes.num = num_nodes;
    mesh->nodes.num_inner = num_inner;
    mesh->nodes.global = compute_globals(mesh, num_cells, dims, coords, neighbor);
    mesh->nodes.coord = coord;
}

static int num_cells_side(Triple num_cells, int idx)
{
    switch (idx / 2) {
        case 0: return num_cells.y * num_cells.z;
        case 1: return num_cells.z * num_cells.x;
        case 2: return num_cells.x * num_cells.y;
        default: teal2_error("invalid index (%d)", idx);
    }
}

static int is_outer_side(const int *dims, const int *coords, int idx)
{
    if (idx % 2 == 0) {
        return coords[idx / 2] == 0;
    }
    return coords[idx / 2] == dims[idx / 2] - 1;
}

static void create_cells(Mesh2 *mesh, Triple num_cells, const int *dims, const int *coords,
                         const int *neighbor, int ndims)
{
    int num_inner = num_cells.x * num_cells.y * num_cells.z;
    int num_boundary = 0;
    int num_periodic = 0;
    int num_neighbor = 0;
    for (int i = 0; i < 2 * ndims; i++) {
        if (neighbor[i] == MPI_PROC_NULL) {
            num_boundary += num_cells_side(num_cells, i);
        }
        else if (is_outer_side(dims, coords, i)) {
            num_periodic += num_cells_side(num_cells, i);
        }
        else {
            num_neighbor += num_cells_side(num_cells, i);
        }
    }
    mesh->cells.num = num_inner + num_boundary + num_periodic + num_neighbor;

    int *node_off = teal2_calloc(mesh->cells.num + 1, sizeof(*node_off));
    int *node_idx = teal2_calloc(mesh->cells.num * MAX_CELL_NODES, sizeof(*node_idx));

    int nnx = num_cells.x + 1;
    int nny = num_cells.y + 1;

    int num = 0;
    for (int k = 0; k < num_cells.z; k++) {
        for (int j = 0; j < num_cells.y; j++) {
            for (int i = 0; i < num_cells.x; i++) {
                node_off[num + 1] = node_off[num];
                node_idx[node_off[num + 1]++] = (i + 0) + (nnx * ((j + 0) + nny * (k + 0)));
                node_idx[node_off[num + 1]++] = (i + 1) + (nnx * ((j + 0) + nny * (k + 0)));
                node_idx[node_off[num + 1]++] = (i + 1) + (nnx * ((j + 1) + nny * (k + 0)));
                node_idx[node_off[num + 1]++] = (i + 0) + (nnx * ((j + 1) + nny * (k + 0)));
                node_idx[node_off[num + 1]++] = (i + 0) + (nnx * ((j + 0) + nny * (k + 1)));
                node_idx[node_off[num + 1]++] = (i + 1) + (nnx * ((j + 0) + nny * (k + 1)));
                node_idx[node_off[num + 1]++] = (i + 1) + (nnx * ((j + 1) + nny * (k + 1)));
                node_idx[node_off[num + 1]++] = (i + 0) + (nnx * ((j + 1) + nny * (k + 1)));
                num += 1;
            }
        }
    }
    if (ndims >= 1) {
        for (int i = 0; i <= num_cells.x; i += num_cells.x) {
            for (int j = 0; j < num_cells.y; j++) {
                for (int k = 0; k < num_cells.z; k++) {
                    node_off[num + 1] = node_off[num];
                    node_idx[node_off[num + 1]++] = (i + 0) + (nnx * ((j + 0) + nny * (k + 0)));
                    node_idx[node_off[num + 1]++] = (i + 0) + (nnx * ((j + 1) + nny * (k + 0)));
                    node_idx[node_off[num + 1]++] = (i + 0) + (nnx * ((j + 1) + nny * (k + 1)));
                    node_idx[node_off[num + 1]++] = (i + 0) + (nnx * ((j + 0) + nny * (k + 1)));
                    num += 1;
                }
            }
        }
    }
    if (ndims >= 2) {
        for (int j = 0; j <= num_cells.y; j += num_cells.y) {
            for (int k = 0; k < num_cells.z; k++) {
                for (int i = 0; i < num_cells.x; i++) {
                    node_off[num + 1] = node_off[num];
                    node_idx[node_off[num + 1]++] = (i + 0) + (nnx * ((j + 0) + nny * (k + 0)));
                    node_idx[node_off[num + 1]++] = (i + 1) + (nnx * ((j + 0) + nny * (k + 0)));
                    node_idx[node_off[num + 1]++] = (i + 1) + (nnx * ((j + 0) + nny * (k + 1)));
                    node_idx[node_off[num + 1]++] = (i + 0) + (nnx * ((j + 0) + nny * (k + 1)));
                    num += 1;
                }
            }
        }
    }
    if (ndims >= 3) {
        for (int k = 0; k <= num_cells.z; k += num_cells.z) {
            for (int i = 0; i < num_cells.x; i++) {
                for (int j = 0; j < num_cells.y; j++) {
                    node_off[num + 1] = node_off[num];
                    node_idx[node_off[num + 1]++] = (i + 0) + (nnx * ((j + 0) + nny * (k + 0)));
                    node_idx[node_off[num + 1]++] = (i + 1) + (nnx * ((j + 0) + nny * (k + 0)));
                    node_idx[node_off[num + 1]++] = (i + 1) + (nnx * ((j + 1) + nny * (k + 0)));
                    node_idx[node_off[num + 1]++] = (i + 0) + (nnx * ((j + 1) + nny * (k + 0)));
                    num += 1;
                }
            }
        }
    }
    assert(num == mesh->cells.num);

    mesh->cells.num_inner = num_inner;
    mesh->cells.off_boundary = num_inner + num_boundary;
    mesh->cells.off_periodic = num_inner + num_boundary + num_periodic;
    mesh->cells.node.off = node_off;
    mesh->cells.node.idx = teal2_realloc(node_idx, node_off[mesh->cells.num], sizeof(*node_idx));
}

static Vector compute_translation(Vector del_coord, int idx)
{
    double sign = (idx % 2 == 0) ? -1 : +1;
    switch (idx / 2) {
        case 0: return (Vector){sign * del_coord.x, 0, 0};
        case 1: return (Vector){0, sign * del_coord.y, 0};
        case 2: return (Vector){0, 0, sign * del_coord.z};
        default: teal2_error("invalid index (%d)", idx);
    }
}

static void create_entities(Mesh2 *mesh, Triple num_cells, Vector del_coord, const int *dims,
                            const int *periods, const int *coords, const int *neighbor, int ndims)
{
    int num_inner = 1;
    int num_entities = num_inner + (2 * ndims);

    Name *name = teal2_calloc(num_entities, sizeof(*name));
    int *cell_off = teal2_calloc(num_entities + 1, sizeof(*cell_off));
    Matrix *rotation = teal2_calloc(num_entities, sizeof(*rotation));
    Vector *translation = teal2_calloc(num_entities, sizeof(*translation));

    char *entity[] = {"domain", "left", "right", "bottom", "top", "back", "front"};
    strcpy(name[0], entity[0]);
    cell_off[1] = mesh->cells.num_inner;

    int num = 1;
    int num_boundary = 0;
    for (int i = 0; i < 2 * ndims; i++) {
        if (!periods[i / 2]) {
            strcpy(name[num], entity[i + 1]);
            cell_off[num + 1] = cell_off[num];
            if (neighbor[i] == MPI_PROC_NULL) {
                cell_off[num + 1] += num_cells_side(num_cells, i);
            }
            num_boundary += 1;
            num += 1;
        }
    }
    for (int i = 0; i < 2 * ndims; i++) {
        if (periods[i / 2]) {
            strcpy(name[num], entity[i + 1]);
            cell_off[num + 1] = cell_off[num];
            if (is_outer_side(dims, coords, i)) {
                cell_off[num + 1] += num_cells_side(num_cells, i);
            }
            rotation[num].x.x = rotation[num].y.y = rotation[num].z.z = 1;
            translation[num] = compute_translation(del_coord, i);
            num += 1;
        }
    }
    assert(num == num_entities);

    mesh->entities.num = num_entities;
    mesh->entities.num_inner = num_inner;
    mesh->entities.off_boundary = num_inner + num_boundary;
    mesh->entities.name = name;
    mesh->entities.cell_off = cell_off;
    mesh->entities.rotation = rotation;
    mesh->entities.translation = translation;
}

static void append_side_send_idx(int *send_idx, int *send_off, Triple num_cells, int idx)
{
    switch (idx / 2) {
        case 0: {
            int i = (idx % 2 == 0) ? 0 : num_cells.x - 1;  // NOLINT(readability-identifier-length)
            for (int j = 0; j < num_cells.y; j++) {
                for (int k = 0; k < num_cells.z; k++) {
                    send_idx[(*send_off)++] = i + (num_cells.x * (j + (num_cells.y * k)));
                }
            }
            return;
        }
        case 1: {
            int j = (idx % 2 == 0) ? 0 : num_cells.y - 1;  // NOLINT(readability-identifier-length)
            for (int k = 0; k < num_cells.z; k++) {
                for (int i = 0; i < num_cells.x; i++) {
                    send_idx[(*send_off)++] = i + (num_cells.x * (j + (num_cells.y * k)));
                }
            }
            return;
        }
        case 2: {
            int k = (idx % 2 == 0) ? 0 : num_cells.z - 1;  // NOLINT(readability-identifier-length)
            for (int i = 0; i < num_cells.x; i++) {
                for (int j = 0; j < num_cells.y; j++) {
                    send_idx[(*send_off)++] = i + (num_cells.x * (j + (num_cells.y * k)));
                }
            }
            return;
        }
        default: teal2_error("invalid index (%d)", idx);
    }
}

static void create_neighbors(Mesh2 *mesh, Triple num_cells, const int *dims, const int *coords,
                             const int *neighbor, int ndims)
{
    int num_neighbors = 0;
    int num_indices = 0;
    for (int i = 0; i < 2 * ndims; i++) {
        if (neighbor[i] != MPI_PROC_NULL) {
            num_neighbors += 1;
            num_indices += num_cells_side(num_cells, i);
        }
    }

    int (*tag)[2] = teal2_calloc(num_neighbors, sizeof(*tag));
    int *rank = teal2_calloc(num_neighbors, sizeof(*rank));
    int *recv_off = teal2_calloc(num_neighbors + 1, sizeof(*recv_off));
    int *send_off = teal2_calloc(num_neighbors + 1, sizeof(*send_off));
    int *send_idx = teal2_calloc(num_indices, sizeof(*send_idx));

    recv_off[0] = mesh->cells.off_boundary;

    int num = 0;
    for (int i = 0; i < 2 * ndims; i++) {
        if (neighbor[i] != MPI_PROC_NULL && is_outer_side(dims, coords, i)) {
            tag[num][0] = i;
            tag[num][1] = (i % 2 == 0) ? i + 1 : i - 1;
            rank[num] = neighbor[i];
            recv_off[num + 1] = recv_off[num] + num_cells_side(num_cells, i);
            send_off[num + 1] = send_off[num];
            append_side_send_idx(send_idx, &send_off[num + 1], num_cells, i);
            num += 1;
        }
    }
    for (int i = 0; i < 2 * ndims; i++) {
        if (neighbor[i] != MPI_PROC_NULL && !is_outer_side(dims, coords, i)) {
            tag[num][0] = mesh->entities.num;
            tag[num][1] = mesh->entities.num;
            rank[num] = neighbor[i];
            recv_off[num + 1] = recv_off[num] + num_cells_side(num_cells, i);
            send_off[num + 1] = send_off[num];
            append_side_send_idx(send_idx, &send_off[num + 1], num_cells, i);
            num += 1;
        }
    }
    assert(num == num_neighbors);
    assert(send_off[num_neighbors] == num_indices);

    mesh->neighbors.num = num_neighbors;
    mesh->neighbors.tag = tag;
    mesh->neighbors.rank = rank;
    mesh->neighbors.recv_off = recv_off;
    mesh->neighbors.send.off = send_off;
    mesh->neighbors.send.idx = send_idx;
}

static void reorder_nodes(Mesh2 *mesh, Triple num_cells, const int *coords)
{
    int *map = teal2_calloc(mesh->nodes.num, sizeof(*map));

    int num = 0;
    int num_inner = 0;
    int num_outer = 0;
    for (int k = 0; k < num_cells.z + 1; k++) {
        for (int j = 0; j < num_cells.y + 1; j++) {
            for (int i = 0; i < num_cells.x + 1; i++) {
                if (!(i == 0 && coords[0]) && !(j == 0 && coords[1]) && !(k == 0 && coords[2])) {
                    map[num] = num_inner++;
                }
                else {
                    map[num] = mesh->nodes.num_inner + num_outer++;
                }
                num += 1;
            }
        }
    }
    assert(num == mesh->nodes.num);
    assert(num_inner == mesh->nodes.num_inner);
    assert(num_outer == mesh->nodes.num - mesh->nodes.num_inner);

    mesh2_reorder_nodes(mesh, map, 0, mesh->nodes.num);

    teal2_free(map);
}

static void reorder_cells(Mesh2 *mesh, Triple num_cells, const int *dims, const int *coords,
                          const int *neighbor, int ndims)
{
    int *map = teal2_calloc(mesh->cells.num, sizeof(*map));

    int num = 0;
    int num_inner = 0;
    int num_boundary = 0;
    int num_periodic = 0;
    int num_neighbor = 0;
    for (int i = 0; i < mesh->cells.num_inner; i++) {
        map[num++] = num_inner++;
    }
    for (int i = 0; i < 2 * ndims; i++) {
        for (int j = 0; j < num_cells_side(num_cells, i); j++) {
            if (neighbor[i] == MPI_PROC_NULL) {
                map[num] = mesh->cells.num_inner + num_boundary++;
            }
            else if (is_outer_side(dims, coords, i)) {
                map[num] = mesh->cells.off_boundary + num_periodic++;
            }
            else {
                map[num] = mesh->cells.off_periodic + num_neighbor++;
            }
            num += 1;
        }
    }
    assert(num == mesh->cells.num);
    assert(num_inner == mesh->cells.num_inner);
    assert(num_boundary == mesh->cells.off_boundary - mesh->cells.num_inner);
    assert(num_periodic == mesh->cells.off_periodic - mesh->cells.off_boundary);
    assert(num_neighbor == mesh->cells.num - mesh->cells.off_periodic);

    mesh2_reorder_cells(mesh, map, 0, mesh->cells.num);

    teal2_free(map);
}

Mesh2 *mesh2_create(Vector min_coord, Vector max_coord, Triple num_cells, Triple periodic)
{
    int dims[3];
    int ndims = compute_dims(dims, num_cells);

    int periods[3] = {periodic.x, periodic.y, periodic.z};
    int reorder = 1;  // WARN: I have not found a cluster yet which actually does this
    MPI_Comm comm;
    MPI_Cart_create(sync2.comm, ndims, dims, periods, reorder, &comm);
    sync2_reinit(comm);

    int coords[3] = {0};
    MPI_Cart_coords(sync2.comm, sync2.rank, ndims, coords);

    int neighbor[6];
    memset(neighbor, MPI_PROC_NULL, sizeof(neighbor));
    for (int i = 0; i < ndims; i++) {
        MPI_Cart_shift(sync2.comm, i, 1, &neighbor[(2 * i) + 0], &neighbor[(2 * i) + 1]);
    }

    Vector del_coord = vector2_sub(max_coord, min_coord);
    split_bounds(&min_coord.x, &max_coord.x, &num_cells.x, 0, dims, coords, ndims);
    split_bounds(&min_coord.y, &max_coord.y, &num_cells.y, 1, dims, coords, ndims);
    split_bounds(&min_coord.z, &max_coord.z, &num_cells.z, 2, dims, coords, ndims);

    Mesh2 *mesh = teal2_calloc(1, sizeof(*mesh));
    create_nodes(mesh, min_coord, max_coord, num_cells, dims, coords, neighbor);
    create_cells(mesh, num_cells, dims, coords, neighbor, ndims);
    create_entities(mesh, num_cells, del_coord, dims, periods, coords, neighbor, ndims);
    create_neighbors(mesh, num_cells, dims, coords, neighbor, ndims);

    reorder_nodes(mesh, num_cells, coords);
    reorder_cells(mesh, num_cells, dims, coords, neighbor, ndims);

    return mesh;
}
