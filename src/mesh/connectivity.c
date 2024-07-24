#include "connectivity.h"

#include <metis.h>

#include "core/dict.h"
#include "core/memory.h"
#include "core/utils.h"

static void compute_periodic_cell_to_node(const Mesh *mesh, Dict *periodic);

static void find_inner_faces(Mesh *mesh, const Dict *periodic, Dict *f2n);

static void find_boundary_faces(Mesh *mesh, Dict *f2n);

static void find_sync_faces(Mesh *mesh, const Dict *periodic, Dict *f2n);

static void compute_face_connectivity(Mesh *mesh, const Dict *f2n);

static void compute_entity_to_face(Mesh *mesh);

void connectivity_cells(Mesh *mesh)
{
    idx_t ne = mesh->n_cells;
    idx_t nn = mesh->n_nodes;
    smart idx_t *eptr = memory_calloc(ne + 1, sizeof(*eptr));
    smart idx_t *eind = memory_calloc(ne * MAX_CELL_NODES, sizeof(*eind));
    for (long j = 0; j < mesh->n_cells; ++j) {
        eptr[j + 1] = eptr[j];
        for (long i = mesh->cell.i_node[j]; i < mesh->cell.i_node[j + 1]; ++i)
            eind[eptr[j + 1]++] = mesh->cell.node[i];

        if (j >= mesh->n_inner_cells)
            eind[eptr[j + 1]++] = nn++;  // add extra node to prevent unintentional connections
    }
    idx_t ncommon = N_DIMS;
    idx_t numflag = 0;
    idx_t *xadj, *adjncy;
    METIS_MeshToDual(&ne, &nn, eptr, eind, &ncommon, &numflag, &xadj, &adjncy);
    mesh->cell.i_cell = xadj;
    mesh->cell.cell = adjncy;
}

void connectivity_faces(Mesh *mesh)
{
    defer(dict_free) Dict periodic = dict_create(mesh->n_ghost_cells);
    compute_periodic_cell_to_node(mesh, &periodic);

    defer(dict_free) Dict f2n = dict_create(mesh->n_cells * MAX_CELL_FACES);
    find_inner_faces(mesh, &periodic, &f2n);
    find_boundary_faces(mesh, &f2n);
    find_sync_faces(mesh, &periodic, &f2n);

    compute_face_connectivity(mesh, &f2n);
    compute_entity_to_face(mesh);
}

long connectivity_nodes(const Mesh *mesh, long *node, long j0, long j1)
{
    // WARNING: only works for unique node indices
    ensure(j0 != j1);
    long n_nodes = 0;
    for (long i0 = mesh->cell.i_node[j0]; i0 < mesh->cell.i_node[j0 + 1]; ++i0) {
        for (long i1 = mesh->cell.i_node[j1]; i1 < mesh->cell.i_node[j1 + 1]; ++i1) {
            if (mesh->cell.node[i0] == mesh->cell.node[i1]) {
                if (node) node[n_nodes] = mesh->cell.node[i0];
                n_nodes += 1;
                break;
            }
        }
    }
    return n_nodes;
}

static void compute_periodic_cell_to_node(const Mesh *mesh, Dict *periodic)
{
    const ALIAS(i_cell, mesh->cell.i_cell);
    const ALIAS(cell, mesh->cell.cell);
    const ALIAS(i_node, mesh->cell.i_node);
    const ALIAS(node, mesh->cell.node);
    for (long j = mesh->n_inner_cells; j < mesh->n_inner_cells + mesh->n_ghost_cells; ++j) {
        if (i_cell[j + 1] - i_cell[j] != N_SIDES) continue;
        const long n_nodes = i_node[j + 1] - i_node[j];
        ensure(n_nodes >= N_DIMS);
        dict_insert(periodic, &cell[i_cell[j]], N_SIDES, &node[i_node[j]], n_nodes);
    }
}

static void find_inner_faces(Mesh *mesh, const Dict *periodic, Dict *f2n)
{
    const ALIAS(i_cell, mesh->cell.i_cell);
    for (long j = 0; j < mesh->n_inner_cells; ++j) {
        for (long i = i_cell[j]; i < i_cell[j + 1]; ++i) {
            if (mesh->cell.cell[i] >= mesh->n_inner_cells) continue;
            const long cell[N_SIDES] = {min(j, mesh->cell.cell[i]), max(j, mesh->cell.cell[i])};
            if (dict_lookup(f2n, cell, N_SIDES)) continue;

            long common[MAX_FACE_NODES];
            const long n_common = connectivity_nodes(mesh, common, cell[L], cell[R]);
            if (n_common >= N_DIMS) dict_insert(f2n, cell, N_SIDES, common, n_common);

            const DictItem *item = dict_lookup(periodic, cell, N_SIDES);
            if (item) dict_insert(f2n, (long[]){cell[R], cell[L]}, N_SIDES, item->val, item->nval);
        }
    }
    mesh->n_inner_faces = f2n->n_items;
}

static void find_boundary_faces(Mesh *mesh, Dict *f2n)
{
    const ALIAS(i_cell, mesh->cell.i_cell);
    for (long j = mesh->n_inner_cells; j < mesh->n_inner_cells + mesh->n_ghost_cells; ++j) {
        if (i_cell[j + 1] - i_cell[j] == N_SIDES) continue;
        for (long i = i_cell[j]; i < i_cell[j + 1]; ++i) {
            const long cell[N_SIDES] = {j, mesh->cell.cell[i]};

            long common[MAX_FACE_NODES];
            const long n_common = connectivity_nodes(mesh, common, cell[L], cell[R]);
            ensure(n_common >= N_DIMS);
            dict_insert(f2n, cell, N_SIDES, common, n_common);
        }
    }
    mesh->n_bound_faces = f2n->n_items - mesh->n_inner_faces;
}

static void find_sync_faces(Mesh *mesh, const Dict *periodic, Dict *f2n)
{
    const ALIAS(i_cell, mesh->cell.i_cell);
    for (long j = mesh->n_inner_cells + mesh->n_ghost_cells; j < mesh->n_cells; ++j) {
        for (long i = i_cell[j]; i < i_cell[j + 1]; ++i) {
            const long cell[N_SIDES] = {j, mesh->cell.cell[i]};

            long common[MAX_FACE_NODES];
            const long n_common = connectivity_nodes(mesh, common, cell[L], cell[R]);
            if (n_common >= N_DIMS) dict_insert(f2n, cell, N_SIDES, common, n_common);

            // if sync face is periodic, we have to look for the nodes in reverse order because
            // periodic dict only contains the connection from inner to sync cell
            const DictItem *item = dict_lookup(periodic, (long[]){cell[R], cell[L]}, N_SIDES);
            if (item) dict_insert(f2n, (long[]){cell[R], cell[L]}, N_SIDES, item->val, item->nval);
        }
    }
    mesh->n_faces = f2n->n_items;
}

static void compute_face_connectivity(Mesh *mesh, const Dict *f2n)
{
    long *i_node = memory_calloc(mesh->n_faces + 1, sizeof(*i_node));
    long *node = memory_calloc(mesh->n_faces * MAX_FACE_NODES, sizeof(*node));
    long(*cell)[N_SIDES] = memory_calloc(mesh->n_faces, sizeof(*cell));
    smart DictItem *item = dict_serialize_by_index(f2n);
    for (long j = 0; j < mesh->n_faces; ++j) {
        ensure(item[j].nval >= N_DIMS);
        i_node[j + 1] = i_node[j];
        for (long i = 0; i < item[j].nval; ++i) node[i_node[j + 1]++] = item[j].val[i];

        ensure(item[j].nkey == N_SIDES && item[j].key[L] != item[j].key[R]);
        cell[j][L] = min(item[j].key[L], item[j].key[R]);
        cell[j][R] = max(item[j].key[L], item[j].key[R]);
    }
    mesh->face.i_node = i_node;
    mesh->face.node = memory_realloc(node, i_node[mesh->n_faces], sizeof(*node));
    mesh->face.cell = cell;
}

static void compute_entity_to_face(Mesh *mesh)
{
    const ALIAS(j_cell, mesh->entity.j_cell);
    const ALIAS(i_cell, mesh->cell.i_cell);
    long *j_face = memory_calloc(mesh->n_entities + 1, sizeof(*j_face));
    j_face[0] = mesh->n_inner_faces;
    for (long e = 0; e < mesh->n_entities; ++e) {
        j_face[e + 1] += j_face[e];
        if (j_cell[e] < mesh->n_inner_cells) continue;
        if (i_cell[j_cell[e] + 1] - i_cell[j_cell[e]] == N_SIDES) continue;
        j_face[e + 1] += j_cell[e + 1] - j_cell[e];
    }
    mesh->entity.j_face = j_face;
}
