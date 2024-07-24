#include "read_gmsh.h"

#include <gmshc.h>
#include <stdio.h>
#include <string.h>

#include "core/memory.h"
#include "core/utils.h"

static void create_nodes(Mesh *mesh, long **node_idx);

static void create_cells(Mesh *mesh, const long *node_idx);

static void create_entities(Mesh *mesh);

static long *tag_2_idx(const size_t *tag, long n_tags);

void read_gmsh(Mesh *mesh, const char *fname)
{
    int ierr;

    gmshInitialize(0, 0, 0, 0, &ierr);
    gmshOpen(fname, &ierr);

    const int n_dims = gmshModelGetDimension(&ierr);
    if (n_dims != N_DIMS) error("unsupported mesh dimension '%d'", n_dims);

    const char *ext = strrchr(fname, '.');
    if (ext && !strcmp(ext, ".geo")) gmshModelMeshGenerate(N_DIMS, &ierr);

    smart long *node_idx;
    create_nodes(mesh, &node_idx);
    create_cells(mesh, node_idx);
    create_entities(mesh);

    gmshFinalize(&ierr);
}

static void create_nodes(Mesh *mesh, long **node_idx)
{
    int ierr;

    size_t *node_tag, n_node_tags;
    double *coord;
    gmshModelMeshGetNodes(&node_tag, &n_node_tags, &coord, 0, 0, 0, N_DIMS, -1, 1, 0, &ierr);

    *node_idx = tag_2_idx(node_tag, n_node_tags);
    gmshFree(node_tag);

    mesh->n_nodes = n_node_tags;
    mesh->node.coord = (void *)coord;
}

static void create_cells(Mesh *mesh, const long *node_idx)
{
    int ierr;

    long n_cells = 0, n_cell_nodes = 0;
    long *i_node = memory_calloc(n_cells + 1, sizeof(*i_node));
    long *node = 0;

    for (long n = 0, dim = N_DIMS; dim >= N_DIMS - 1; --dim) {
        int *ptag;
        size_t n_ptags;
        gmshModelGetPhysicalGroups(&ptag, &n_ptags, dim, &ierr);
        for (size_t p = 1; p < n_ptags; p += 2) {
            int *etag;
            size_t n_etags;
            gmshModelGetEntitiesForPhysicalGroup(dim, ptag[p], &etag, &n_etags, &ierr);
            for (size_t e = 0; e < n_etags; ++e) {
                int *element_type;
                size_t n_element_types, **element_tag, *n_element_tags, **node_tag, *n_node_tags;
                gmshModelMeshGetElements(&element_type, &n_element_types, &element_tag,
                                         &n_element_tags, 0, &node_tag, &n_node_tags, 0, dim,
                                         etag[e], &ierr);
                for (size_t t = 0; t < n_element_types; ++t) {
                    char *element_name;
                    int element_dim, element_order, n_element_nodes, n_element_primary_nodes;
                    gmshModelMeshGetElementProperties(element_type[t], &element_name, &element_dim,
                                                      &element_order, &n_element_nodes, 0, 0,
                                                      &n_element_primary_nodes, &ierr);
                    if (element_order != 1) error("unsupported element order '%ld'", element_order);
                    gmshFree(element_name);

                    n_cells += n_element_tags[t];
                    n_cell_nodes += n_node_tags[t];
                    i_node = memory_realloc(i_node, n_cells + 1, sizeof(*i_node));
                    node = memory_realloc(node, n_cell_nodes, sizeof(*node));

                    for (size_t j = 0; j < n_element_tags[t]; ++j) {
                        i_node[n + 1] = i_node[n];
                        for (long i = 0; i < n_element_nodes; ++i)
                            node[i_node[n + 1]++] = node_idx[node_tag[t][j * n_element_nodes + i]];
                        n += 1;
                    }
                    gmshFree(node_tag[t]);
                    gmshFree(element_tag[t]);
                }
                gmshFree(n_node_tags);
                gmshFree(node_tag);
                gmshFree(n_element_tags);
                gmshFree(element_tag);
                gmshFree(element_type);
            }
            gmshFree(etag);
        }
        if (dim == N_DIMS) mesh->n_inner_cells = n_cells;
        if (dim == N_DIMS - 1) mesh->n_ghost_cells = n_cells - mesh->n_inner_cells;
        gmshFree(ptag);
    }

    mesh->n_cells = n_cells;
    mesh->cell.i_node = i_node;
    mesh->cell.node = node;
}

static void create_entities(Mesh *mesh)
{
    int ierr;

    long n_entities = 0;
    char(*name)[NAMELEN] = 0;
    long *j_cell = memory_calloc(n_entities + 1, sizeof(*j_cell));

    for (long dim = N_DIMS; dim >= N_DIMS - 1; --dim) {
        int *ptag;
        size_t n_ptags;
        gmshModelGetPhysicalGroups(&ptag, &n_ptags, dim, &ierr);
        for (size_t p = 1; p < n_ptags; p += 2) {
            char *pname;
            gmshModelGetPhysicalName(dim, ptag[p], &pname, &ierr);

            name = memory_realloc(name, n_entities + 1, sizeof(*name));
            j_cell = memory_realloc(j_cell, n_entities + 2, sizeof(*j_cell));
            if (strlen(pname))
                strcpy(name[n_entities], pname);
            else
                sprintf(name[n_entities], "%d", ptag[p]);
            j_cell[n_entities + 1] = j_cell[n_entities];
            gmshFree(pname);

            int *etag;
            size_t n_etags;
            gmshModelGetEntitiesForPhysicalGroup(dim, ptag[p], &etag, &n_etags, &ierr);
            for (size_t e = 0; e < n_etags; ++e) {
                size_t *n_ghost_tags, nn_ghost_tags;
                gmshModelMeshGetElements(0, 0, 0, &n_ghost_tags, &nn_ghost_tags, 0, 0, 0, dim,
                                         etag[e], &ierr);
                for (size_t t = 0; t < nn_ghost_tags; ++t)
                    j_cell[n_entities + 1] += n_ghost_tags[t];
                gmshFree(n_ghost_tags);
            }
            gmshFree(etag);

            n_entities += 1;
        }
        gmshFree(ptag);
    }

    mesh->n_entities = n_entities;
    mesh->entity.name = name;
    mesh->entity.j_cell = j_cell;
    mesh->entity.offset = memory_calloc(mesh->n_entities, sizeof(*mesh->entity.offset));
}

static long *tag_2_idx(const size_t *tag, long n_tags)
{
    size_t max_tag = tag[0];
    for (long i = 1; i < n_tags; ++i)
        if (tag[i] > max_tag) max_tag = tag[i];

    long *idx = memory_calloc(max_tag + 1, sizeof(*idx));
    for (size_t i = 0; i < max_tag + 1; ++i) idx[i] = -1;
    for (long i = 0; i < n_tags; ++i) idx[tag[i]] = i;

    return idx;
}
