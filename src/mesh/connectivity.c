#include "connectivity.h"

#include <assert.h>
#include <metis.h>

#include "periodic.h"
#include "teal/dict.h"
#include "teal/memory.h"
#include "teal/utils.h"

static void find_inner_faces(Mesh *mesh, const Dict *periodic, Dict *f2n);

static void find_ghost_faces(Mesh *mesh, Dict *f2n);

static void find_sync_faces(Mesh *mesh, const Dict *periodic, Dict *f2n);

static void compute_face_connectivity(Mesh *mesh, const Dict *f2n);

static void compute_entity_to_face(Mesh *mesh);

long connectivity_nodes(const Mesh *mesh, long *node, long cell_a, long cell_b)
{
    assert(cell_a != cell_b);
    long n_nodes = 0;
    for (long i = mesh->cell.i_node[cell_a]; i < mesh->cell.i_node[cell_a + 1]; ++i) {
        for (long j = mesh->cell.i_node[cell_b]; j < mesh->cell.i_node[cell_b + 1]; ++j) {
            if (mesh->cell.node[i] == mesh->cell.node[j]) {
                if (node) node[n_nodes] = mesh->cell.node[i];
                n_nodes += 1;
                break;  // WARNING: only works for unique node indices
            }
        }
    }
    return n_nodes;
}

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
    mesh->cell.i_cell = memory_duplicate(xadj, ne + 1, sizeof(*xadj));
    mesh->cell.cell = memory_duplicate(adjncy, xadj[ne], sizeof(*adjncy));
    METIS_Free(xadj);
    METIS_Free(adjncy);
}

void connectivity_faces(Mesh *mesh)
{
    defer(dict_free) Dict *periodic = dict_create(mesh->n_ghost_cells);
    periodic_cell_to_node(mesh, periodic);

    defer(dict_free) Dict *f2n = dict_create(mesh->n_cells * MAX_CELL_FACES);
    find_inner_faces(mesh, periodic, f2n);
    find_ghost_faces(mesh, f2n);
    find_sync_faces(mesh, periodic, f2n);

    compute_face_connectivity(mesh, f2n);
    compute_entity_to_face(mesh);
}

static void find_inner_faces(Mesh *mesh, const Dict *periodic, Dict *f2n)
{
    for (long j = 0; j < mesh->n_inner_cells; ++j) {
        for (long i = mesh->cell.i_cell[j]; i < mesh->cell.i_cell[j + 1]; ++i) {
            if (mesh->cell.cell[i] >= mesh->n_inner_cells) continue;
            const Vector2l cell = {min(j, mesh->cell.cell[i]), max(j, mesh->cell.cell[i])};
            if (dict_lookup(f2n, cell, N_SIDES)) continue;

            long common[MAX_FACE_NODES];
            const long n_common = connectivity_nodes(mesh, common, cell[L], cell[R]);
            if (n_common >= N_DIMS) dict_insert(f2n, cell, common, N_SIDES, n_common);

            const DictItem *item = dict_lookup(periodic, cell, N_SIDES);
            if (item)
                dict_insert(f2n, (Vector2l){cell[R], cell[L]}, item->val, N_SIDES, item->nval);
        }
    }
    mesh->n_inner_faces = f2n->n_items;
}

static void find_ghost_faces(Mesh *mesh, Dict *f2n)
{
    for (long j = mesh->n_inner_cells; j < mesh->n_inner_cells + mesh->n_ghost_cells; ++j) {
        if (mesh->cell.i_cell[j + 1] - mesh->cell.i_cell[j] == N_SIDES) continue;
        for (long i = mesh->cell.i_cell[j]; i < mesh->cell.i_cell[j + 1]; ++i) {
            const Vector2l cell = {j, mesh->cell.cell[i]};

            long common[MAX_FACE_NODES];
            const long n_common = connectivity_nodes(mesh, common, cell[L], cell[R]);
            assert(n_common >= N_DIMS);
            dict_insert(f2n, cell, common, N_SIDES, n_common);
        }
    }
    mesh->n_ghost_faces = f2n->n_items - mesh->n_inner_faces;
}

static void find_sync_faces(Mesh *mesh, const Dict *periodic, Dict *f2n)
{
    for (long j = mesh->n_inner_cells + mesh->n_ghost_cells; j < mesh->n_cells; ++j) {
        for (long i = mesh->cell.i_cell[j]; i < mesh->cell.i_cell[j + 1]; ++i) {
            const Vector2l cell = {j, mesh->cell.cell[i]};

            long common[MAX_FACE_NODES];
            const long n_common = connectivity_nodes(mesh, common, cell[L], cell[R]);
            if (n_common >= N_DIMS) dict_insert(f2n, cell, common, N_SIDES, n_common);

            // if sync face is periodic, we have to look for the nodes in reverse order because
            // periodic dict only contains the connection from inner to sync cell
            const DictItem *item = dict_lookup(periodic, (Vector2l){cell[R], cell[L]}, N_SIDES);
            if (item)
                dict_insert(f2n, (Vector2l){cell[R], cell[L]}, item->val, N_SIDES, item->nval);
        }
    }
    mesh->n_faces = f2n->n_items;
}

static void compute_face_connectivity(Mesh *mesh, const Dict *f2n)
{
    long *i_node = memory_calloc(mesh->n_faces + 1, sizeof(*i_node));
    long *node = memory_calloc(mesh->n_faces * MAX_FACE_NODES, sizeof(*node));
    Vector2l *cell = memory_calloc(mesh->n_faces, sizeof(*cell));
    smart DictItem *item = dict_serialize_by_index(f2n);
    for (long j = 0; j < mesh->n_faces; ++j) {
        assert(item[j].nval >= N_DIMS);
        i_node[j + 1] = i_node[j];
        for (long i = 0; i < item[j].nval; ++i) node[i_node[j + 1]++] = item[j].val[i];

        assert(item[j].nkey == N_SIDES && item[j].key[L] != item[j].key[R]);
        cell[j][L] = min(item[j].key[L], item[j].key[R]);
        cell[j][R] = max(item[j].key[L], item[j].key[R]);
    }
    mesh->face.i_node = i_node;
    mesh->face.node = memory_realloc(node, i_node[mesh->n_faces], sizeof(*node));
    mesh->face.cell = cell;
}

static void compute_entity_to_face(Mesh *mesh)
{
    const alias(j_cell, mesh->entity.j_cell);
    const alias(i_cell, mesh->cell.i_cell);
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
