#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "mesh.h"
#include "reorder.h"
#include "teal/arena.h"
#include "teal/sync.h"
#include "teal/utils.h"
#include "teal/vector.h"

/* Factor `sync.size` into a Cartesian grid, biased by global cell counts to balance work. */
static void compute_dims(tuple num_cells, int ndims, int *dims)
{
    int cells[3] = {lmax(1, num_cells.x), lmax(1, num_cells.y), lmax(1, num_cells.z)};
    dims[0] = dims[1] = dims[2] = 1;

    // sort axes by descending global cell count (try larger cuts first)
    int perm[3] = {0, 1, 2};
    if (cells[perm[1]] > cells[perm[0]]) {
        swap_int(&perm[0], &perm[1]);
    }
    if (cells[perm[2]] > cells[perm[1]]) {
        swap_int(&perm[1], &perm[2]);
    }
    if (cells[perm[1]] > cells[perm[0]]) {
        swap_int(&perm[0], &perm[1]);
    }

    int prime = sync.size;
    while (prime > 1) {
        int factor = 2;
        while ((factor <= prime / factor) && (prime % factor)) {
            factor += 1;
        }
        if (factor > prime / factor) {
            factor = prime;
        }

        // choose axis with largest cell count
        scalar best_score = -1;
        int best_axis = -1;
        for (int i = 0; i < ndims; i++) {
            int axis = perm[i];
            if (dims[axis] > cells[axis] / factor) {
                continue;  // would over-subdivide axis
            }
            scalar score = (scalar)cells[axis] / dims[axis];
            if (score > best_score) {
                best_score = score;
                best_axis = axis;
            }
        }
        assert(best_axis != -1);
        dims[best_axis] *= factor;
        prime /= factor;
    }
}

/* Split physical bounds and cell counts for this rank's Cartesian slice. If axis is inactive,
 * collapse that axis to 1 cell; otherwise compute the rank's local sub-interval and cells count. */
static void split_bounds(scalar *min_coord, scalar *max_coord, int *num_cells, int dim, int ndims,
                         const int *dims, const int *coords)
{
    if (ndims > dim) {
        assert(!isclose(*min_coord, *max_coord));
        scalar width = (*max_coord - *min_coord) / *num_cells;
        int base = *num_cells / dims[dim];
        int extra = *num_cells % dims[dim];
        int offset = (coords[dim] * base) + lmin(coords[dim], extra);
        *num_cells = base + (coords[dim] < extra);
        *min_coord = offset * width;
        *max_coord = *min_coord + (*num_cells * width);
    }
    else {
        *num_cells = 1;
        if (isclose(*min_coord, *max_coord)) {
            *max_coord = *min_coord + 1;  // non-zero extent
        }
    }
}

/* Allocate node coordinates of the local Cartesian grid and count inner nodes. */
static void create_nodes(MeshNodes *nodes, vector min_coord, vector max_coord, tuple num_cells,
                         tuple num_nodes, const int *coords)
{
    nodes->num = num_nodes.x * num_nodes.y * num_nodes.z;
    vector *coord = arena_malloc(nodes->num, sizeof(*coord));
    int num = 0;
    int num_inner = 0;
    for (int k = 0; k < num_nodes.z; k++) {
        for (int j = 0; j < num_nodes.y; j++) {
            for (int i = 0; i < num_nodes.x; i++) {
                if (!(i == 0 && coords[0]) && !(j == 0 && coords[1]) && !(k == 0 && coords[2])) {
                    num_inner += 1;
                }
                coord[num].x = min_coord.x + (i * (max_coord.x - min_coord.x) / num_cells.x);
                coord[num].y = min_coord.y + (j * (max_coord.y - min_coord.y) / num_cells.y);
                coord[num].z = min_coord.z + (k * (max_coord.z - min_coord.z) / num_cells.z);
                num += 1;
            }
        }
    }
    assert(num == nodes->num);
    nodes->num_inner = num_inner;
    nodes->coord = coord;
}

/* Fill local inner nodes lexicographically, then exchange boundary indices with neighbor ranks. */
static void compute_globals(MeshNodes *nodes, tuple num_nodes, const int *dims, const int *coords,
                            const int *neighbor)
{
    int *global = arena_malloc(nodes->num, sizeof(*global));
    int off_inner = sync_lexsum(nodes->num_inner);
    int tag_x = sync_tag();
    int tag_y = sync_tag();
    int tag_z = sync_tag();
    int num = 0;
    int num_inner = 0;
    for (int k = 0; k < num_nodes.z; k++) {
        for (int j = 0; j < num_nodes.y; j++) {
            for (int i = 0; i < num_nodes.x; i++) {
                if (!(i == 0 && coords[0]) && !(j == 0 && coords[1]) && !(k == 0 && coords[2])) {
                    global[num] = off_inner + num_inner++;
                }
                if (i == 0 && coords[0]) {
                    MPI_Recv(&global[num], 1, MPI_INT, neighbor[0], tag_x, sync.comm,
                             MPI_STATUS_IGNORE);
                }
                if (j == 0 && coords[1]) {
                    MPI_Recv(&global[num], 1, MPI_INT, neighbor[2], tag_y, sync.comm,
                             MPI_STATUS_IGNORE);
                }
                if (k == 0 && coords[2]) {
                    MPI_Recv(&global[num], 1, MPI_INT, neighbor[4], tag_z, sync.comm,
                             MPI_STATUS_IGNORE);
                }
                if (i == num_nodes.x - 1 && coords[0] < dims[0] - 1) {
                    MPI_Send(&global[num], 1, MPI_INT, neighbor[1], tag_x, sync.comm);
                }
                if (j == num_nodes.y - 1 && coords[1] < dims[1] - 1) {
                    MPI_Send(&global[num], 1, MPI_INT, neighbor[3], tag_y, sync.comm);
                }
                if (k == num_nodes.z - 1 && coords[2] < dims[2] - 1) {
                    MPI_Send(&global[num], 1, MPI_INT, neighbor[5], tag_z, sync.comm);
                }
                num += 1;
            }
        }
    }
    assert(num == nodes->num);
    assert(num_inner == nodes->num_inner);
    nodes->global = global;
}

static int num_cells_side(tuple num_cells, int idx)
{
    switch (idx / 2) {
        case 0: return num_cells.y * num_cells.z;
        case 1: return num_cells.z * num_cells.x;
        case 2: return num_cells.x * num_cells.y;
        default: error("invalid index -- '%d'", idx);
    }
}

static bool is_edge_side(const int *dims, const int *coords, int idx)
{
    return (idx % 2 == 0) ? (coords[idx / 2] == 0) : (coords[idx / 2] == dims[idx / 2] - 1);
}

/* Build cell-to-node connectivity for inner, ghost, periodic, and neighbor cells. */
static void create_cells(MeshCells *cells, tuple num_cells, tuple num_nodes, int ndims,
                         const int *dims, const int *coords, const int *neighbor)
{
    cells->num_inner = num_cells.x * num_cells.y * num_cells.z;
    cells->num = cells->num_inner;
    int num_ghost = 0;
    int num_periodic = 0;
    for (int i = 0; i < 2 * ndims; i++) {
        if (neighbor[i] == MPI_PROC_NULL) {
            num_ghost += num_cells_side(num_cells, i);
        }
        else if (is_edge_side(dims, coords, i)) {
            num_periodic += num_cells_side(num_cells, i);
        }
        cells->num += num_cells_side(num_cells, i);
    }
    cells->off_ghost = cells->num_inner + num_ghost;
    cells->off_periodic = cells->num_inner + num_ghost + num_periodic;

    int *off = arena_malloc(cells->num + 1, sizeof(*off));
    int *idx = arena_malloc(cells->num * MAX_CELL_NODES, sizeof(*idx));

    off[0] = 0;

    int num = 0;
    for (int k = 0; k < num_cells.z; k++) {
        for (int j = 0; j < num_cells.y; j++) {
            for (int i = 0; i < num_cells.x; i++) {
                off[num + 1] = off[num];
                idx[off[num + 1]++] = (i + 0) + (num_nodes.x * ((j + 0) + num_nodes.y * (k + 0)));
                idx[off[num + 1]++] = (i + 1) + (num_nodes.x * ((j + 0) + num_nodes.y * (k + 0)));
                idx[off[num + 1]++] = (i + 1) + (num_nodes.x * ((j + 1) + num_nodes.y * (k + 0)));
                idx[off[num + 1]++] = (i + 0) + (num_nodes.x * ((j + 1) + num_nodes.y * (k + 0)));
                idx[off[num + 1]++] = (i + 0) + (num_nodes.x * ((j + 0) + num_nodes.y * (k + 1)));
                idx[off[num + 1]++] = (i + 1) + (num_nodes.x * ((j + 0) + num_nodes.y * (k + 1)));
                idx[off[num + 1]++] = (i + 1) + (num_nodes.x * ((j + 1) + num_nodes.y * (k + 1)));
                idx[off[num + 1]++] = (i + 0) + (num_nodes.x * ((j + 1) + num_nodes.y * (k + 1)));
                num += 1;
            }
        }
    }
    if (ndims >= 1) {
        for (int i = 0; i <= num_cells.x; i += num_cells.x) {
            for (int j = 0; j < num_cells.y; j++) {
                for (int k = 0; k < num_cells.z; k++) {
                    off[num + 1] = off[num];
                    idx[off[num + 1]++] =
                        (i + 0) + (num_nodes.x * ((j + 0) + num_nodes.y * (k + 0)));
                    idx[off[num + 1]++] =
                        (i + 0) + (num_nodes.x * ((j + 1) + num_nodes.y * (k + 0)));
                    idx[off[num + 1]++] =
                        (i + 0) + (num_nodes.x * ((j + 1) + num_nodes.y * (k + 1)));
                    idx[off[num + 1]++] =
                        (i + 0) + (num_nodes.x * ((j + 0) + num_nodes.y * (k + 1)));
                    num += 1;
                }
            }
        }
    }
    if (ndims >= 2) {
        for (int j = 0; j <= num_cells.y; j += num_cells.y) {
            for (int k = 0; k < num_cells.z; k++) {
                for (int i = 0; i < num_cells.x; i++) {
                    off[num + 1] = off[num];
                    idx[off[num + 1]++] =
                        (i + 0) + (num_nodes.x * ((j + 0) + num_nodes.y * (k + 0)));
                    idx[off[num + 1]++] =
                        (i + 1) + (num_nodes.x * ((j + 0) + num_nodes.y * (k + 0)));
                    idx[off[num + 1]++] =
                        (i + 1) + (num_nodes.x * ((j + 0) + num_nodes.y * (k + 1)));
                    idx[off[num + 1]++] =
                        (i + 0) + (num_nodes.x * ((j + 0) + num_nodes.y * (k + 1)));
                    num += 1;
                }
            }
        }
    }
    if (ndims >= 3) {
        for (int k = 0; k <= num_cells.z; k += num_cells.z) {
            for (int i = 0; i < num_cells.x; i++) {
                for (int j = 0; j < num_cells.y; j++) {
                    off[num + 1] = off[num];
                    idx[off[num + 1]++] =
                        (i + 0) + (num_nodes.x * ((j + 0) + num_nodes.y * (k + 0)));
                    idx[off[num + 1]++] =
                        (i + 1) + (num_nodes.x * ((j + 0) + num_nodes.y * (k + 0)));
                    idx[off[num + 1]++] =
                        (i + 1) + (num_nodes.x * ((j + 1) + num_nodes.y * (k + 0)));
                    idx[off[num + 1]++] =
                        (i + 0) + (num_nodes.x * ((j + 1) + num_nodes.y * (k + 0)));
                    num += 1;
                }
            }
        }
    }
    assert(num == cells->num);

    cells->node.off = off;
    cells->node.idx = arena_resize(idx, off[cells->num], sizeof(*idx));
    assert(cells->node.idx);
}

/* Build node reorder map to [inner, neighbor] for the current rank. */
static void compute_node_map(const MeshNodes *nodes, tuple num_nodes, const int *coords, int *map)
{
    int num = 0;
    int num_inner = 0;
    int num_neighbors = 0;
    for (int k = 0; k < num_nodes.z; k++) {
        for (int j = 0; j < num_nodes.y; j++) {
            for (int i = 0; i < num_nodes.x; i++) {
                if (!(i == 0 && coords[0]) && !(j == 0 && coords[1]) && !(k == 0 && coords[2])) {
                    map[num] = num_inner++;
                }
                else {
                    map[num] = nodes->num_inner + num_neighbors++;
                }
                num += 1;
            }
        }
    }
    assert(num == nodes->num);
    assert(num_inner == nodes->num_inner);
    assert(num_neighbors == num - num_inner);
}

/* Build cell reorder map to [inner, ghost, periodic, neighbor] for the current rank. */
static void compute_cell_map(const MeshCells *cells, tuple num_cells, int ndims, const int *dims,
                             const int *coords, const int *neighbor, int *map)
{
    int num = 0;
    int num_inner = 0;
    int num_ghost = 0;
    int num_periodic = 0;
    int num_neighbors = 0;
    for (int i = 0; i < cells->num_inner; i++) {
        map[num++] = num_inner++;
    }
    for (int i = 0; i < 2 * ndims; i++) {
        for (int j = 0; j < num_cells_side(num_cells, i); j++) {
            if (neighbor[i] == MPI_PROC_NULL) {
                map[num] = cells->num_inner + num_ghost++;
            }
            else if (is_edge_side(dims, coords, i)) {
                map[num] = cells->off_ghost + num_periodic++;
            }
            else {
                map[num] = cells->off_periodic + num_neighbors++;
            }
            num += 1;
        }
    }
    assert(num == cells->num);
    assert(num_inner == cells->num_inner);
    assert(num_ghost == cells->off_ghost - cells->num_inner);
    assert(num_periodic == cells->off_periodic - cells->off_ghost);
    assert(num_neighbors == num - num_inner - num_ghost - num_periodic);
}

/* Reorder nodes then cells; fix connectivity using in-place maps. */
static void reorder(MeshNodes *nodes, MeshCells *cells, tuple num_cells, tuple num_nodes, int ndims,
                    const int *dims, const int *coords, const int *neighbor)
{
    Arena save = arena_save();

    int num = lmax(nodes->num, cells->num);
    int *map = arena_malloc(num, sizeof(*map));

    compute_node_map(nodes, num_nodes, coords, map);
    mesh_reorder_nodes(nodes, cells, map);

    compute_cell_map(cells, num_cells, ndims, dims, coords, neighbor, map);
    mesh_reorder_cells(cells, 0, 0, cells->num, map);

    arena_load(save);
}

static vector compute_translation(vector del_coord, int idx)
{
    scalar sign = (idx % 2) ? -1 : +1;
    switch (idx / 2) {
        case 0: return (vector){sign * del_coord.x, 0, 0};
        case 1: return (vector){0, sign * del_coord.y, 0};
        case 2: return (vector){0, 0, sign * del_coord.z};
        default: error("invalid index -- '%d'", idx);
    }
}

/* Create entities in the order [inner, ghost, periodic]. */
static void create_entities(MeshEntities *entities, tuple num_cells, vector del_coord, int ndims,
                            const int *dims, const int *periods, const int *coords,
                            const int *neighbor)
{
    entities->num_inner = 1;
    entities->num = entities->num_inner + (2 * ndims);

    Name *name = arena_calloc(entities->num, sizeof(*name));
    int *cell_off = arena_malloc(entities->num + 1, sizeof(*cell_off));
    vector *translation = arena_calloc(entities->num, sizeof(*translation));

    strcpy(name[0], "domain");
    cell_off[0] = 0;
    cell_off[1] = num_cells.x * num_cells.y * num_cells.z;

    Name side[6] = {"left", "right", "bottom", "top", "back", "front"};
    int num = 1;
    int num_ghost = 0;
    for (int i = 0; i < 2 * ndims; i++) {
        if (!periods[i / 2]) {
            strcpy(name[num], side[i]);
            cell_off[num + 1] = cell_off[num];
            if (neighbor[i] == MPI_PROC_NULL) {
                cell_off[num + 1] += num_cells_side(num_cells, i);
            }
            num_ghost += 1;
            num += 1;
        }
    }
    for (int i = 0; i < 2 * ndims; i++) {
        if (periods[i / 2]) {
            sprintf(name[num], "%s:%s", side[i], side[(i % 2 == 0) ? (i + 1) : (i - 1)]);
            cell_off[num + 1] = cell_off[num];
            if (is_edge_side(dims, coords, i)) {
                cell_off[num + 1] += num_cells_side(num_cells, i);
            }
            translation[num] = compute_translation(del_coord, i);
            num += 1;
        }
    }
    assert(num == entities->num);

    entities->off_ghost = entities->num_inner + num_ghost;
    entities->name = name;
    entities->cell_off = cell_off;
    entities->translation = translation;
}

/* Create neighbor metadata in the order [periodic, neighbor]. */
static void create_neighbors(const MeshCells *cells, MeshNeighbors *neighbors, tuple num_cells,
                             int ndims, const int *dims, const int *coords, const int *neighbor)
{
    neighbors->num = 0;
    for (int i = 0; i < 2 * ndims; i++) {
        if (neighbor[i] != MPI_PROC_NULL) {
            neighbors->num += 1;
        }
    }

    int *rank = arena_malloc(neighbors->num, sizeof(*rank));
    int *recv_off = arena_malloc(neighbors->num + 1, sizeof(*recv_off));

    recv_off[0] = cells->off_ghost;

    int num = 0;
    for (int i = 0; i < 2 * ndims; i++) {
        if (neighbor[i] != MPI_PROC_NULL && is_edge_side(dims, coords, i)) {
            rank[num] = neighbor[i];
            recv_off[num + 1] = recv_off[num] + num_cells_side(num_cells, i);
            num += 1;
        }
    }
    for (int i = 0; i < 2 * ndims; i++) {
        if (neighbor[i] != MPI_PROC_NULL && !is_edge_side(dims, coords, i)) {
            rank[num] = neighbor[i];
            recv_off[num + 1] = recv_off[num] + num_cells_side(num_cells, i);
            num += 1;
        }
    }
    assert(num == neighbors->num);

    neighbors->rank = rank;
    neighbors->recv_off = recv_off;
}

Mesh *mesh_create(vector min_coord, vector max_coord, tuple num_cells, flags periodic)
{
    assert(sync.size <= lmax(1, num_cells.x) * lmax(1, num_cells.y) * lmax(1, num_cells.z));

    Mesh *mesh = arena_calloc(1, sizeof(*mesh));

    int dims[3] = {0};
    int ndims = (num_cells.x > 1) + (num_cells.y > 1) + (num_cells.z > 1);
    compute_dims(num_cells, ndims, dims);

    int periods[3] = {periodic.x, periodic.y, periodic.z};
    MPI_Comm comm;
    MPI_Cart_create(sync.comm, ndims, dims, periods, true, &comm);
    sync_reinit(comm);

    int coords[3] = {0};
    MPI_Cart_coords(sync.comm, sync.rank, ndims, coords);

    int neighbor[6];
    MPI_Cart_shift(sync.comm, 0, 1, &neighbor[0], &neighbor[1]);
    MPI_Cart_shift(sync.comm, 1, 1, &neighbor[2], &neighbor[3]);
    MPI_Cart_shift(sync.comm, 2, 1, &neighbor[4], &neighbor[5]);

    vector del_coord = vector_sub(max_coord, min_coord);
    split_bounds(&min_coord.x, &max_coord.x, &num_cells.x, 0, ndims, dims, coords);
    split_bounds(&min_coord.y, &max_coord.y, &num_cells.y, 1, ndims, dims, coords);
    split_bounds(&min_coord.z, &max_coord.z, &num_cells.z, 2, ndims, dims, coords);

    tuple num_nodes = {num_cells.x + 1, num_cells.y + 1, num_cells.z + 1};
    create_nodes(&mesh->nodes, min_coord, max_coord, num_cells, num_nodes, coords);
    compute_globals(&mesh->nodes, num_nodes, dims, coords, neighbor);

    create_cells(&mesh->cells, num_cells, num_nodes, ndims, dims, coords, neighbor);
    reorder(&mesh->nodes, &mesh->cells, num_cells, num_nodes, ndims, dims, coords, neighbor);

    create_entities(&mesh->entities, num_cells, del_coord, ndims, dims, periods, coords, neighbor);
    create_neighbors(&mesh->cells, &mesh->neighbors, num_cells, ndims, dims, coords, neighbor);

    return mesh;
}
