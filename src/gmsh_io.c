#include "gmsh_io.h"

#include <assert.h>
#include <gmshc.h>
#include <stdio.h>
#include <string.h>

#include "memory.h"
#include "mesh.h"
#include "utils.h"

static void create_nodes(Mesh *mesh, long *node_map);
static void create_cells(Mesh *mesh, const long *node_map);
static void create_entities(Mesh *mesh);

void gmsh_read(Mesh *mesh, const char *fname)
{
    int ierr;

    gmshInitialize(0, 0, 0, 0, &ierr);
    gmshOptionSetNumber("General.Verbosity", 3, &ierr);
    gmshOpen(fname, &ierr);

    const char *ext = utils_extension(fname);
    if (!strcmp(ext, "geo")) {
        gmshModelGeoSynchronize(&ierr);
        gmshModelOccSynchronize(&ierr);
        gmshModelMeshGenerate(N_DIMS, &ierr);
    }
    assert(gmshModelGetDimension(&ierr) == N_DIMS && "unsupported mesh dimension");

    gmsh_export(mesh);

    gmshFinalize(&ierr);
}

void gmsh_export(Mesh *mesh)
{
    int ierr;

    size_t max_tag;
    gmshModelMeshGetMaxNodeTag(&max_tag, &ierr);
    cleanup long *node_map = memory_calloc(max_tag + 1, sizeof(*node_map));

    create_nodes(mesh, node_map);
    create_cells(mesh, node_map);
    create_entities(mesh);
}

static void create_nodes(Mesh *mesh, long *node_map)
{
    int ierr;

    size_t n_nodeTags, *nodeTags;
    double *coord;
    gmshModelMeshGetNodes(&nodeTags, &n_nodeTags, &coord, 0, 0, 0, -1, -1, 0, 0, &ierr);

    for (size_t i = 0; i < n_nodeTags; ++i) node_map[nodeTags[i]] = i;
    gmshFree(nodeTags);

    double(*x)[3] = TCAST(x, coord);
    mesh->n_nodes = n_nodeTags;
    assert(mesh->n_nodes > 0 && "mesh contains no nodes");
    mesh->node.x = memory_calloc(mesh->n_nodes, sizeof(*mesh->node.x));
    for (long i = 0; i < mesh->n_nodes; ++i)
        for (long d = 0; d < N_DIMS; ++d) mesh->node.x[i][d] = x[i][d];
    gmshFree(coord);
}

static void create_cells(Mesh *mesh, const long *node_map)
{
    int ierr;

    // count number of elements and number of cell nodes
    long n_element_tags[2] = {}, n_node_tags[2] = {};
    for (long n = 0, d = N_DIMS; d >= N_DIMS - 1; --d) {
        size_t dimTags_n;
        int *dimTags;
        gmshModelGetPhysicalGroups(&dimTags, &dimTags_n, d, &ierr);
        for (size_t p = 0; p < dimTags_n; p += 2) {
            size_t tags_n;
            int *tags;
            gmshModelGetEntitiesForPhysicalGroup(dimTags[p], dimTags[p + 1], &tags, &tags_n, &ierr);
            for (size_t e = 0; e < tags_n; ++e) {
                size_t elementTypes_n, *elementTags_n, *nodeTags_n;
                gmshModelMeshGetElements(0, &elementTypes_n, 0, &elementTags_n, 0, 0, &nodeTags_n,
                                         0, d, tags[e], &ierr);
                for (size_t i = 0; i < elementTypes_n; ++i) {
                    n_element_tags[n] += elementTags_n[i];
                    n_node_tags[n] += nodeTags_n[i];
                }
                gmshFree(elementTags_n);
                gmshFree(nodeTags_n);
            }
            gmshFree(tags);
        }
        gmshFree(dimTags);
        n += 1;
    }

    mesh->n_inner_cells = n_element_tags[0];
    mesh->n_ghost_cells = n_element_tags[1];
    mesh->n_cells = mesh->n_inner_cells + mesh->n_ghost_cells;
    assert(mesh->n_cells > 0 && "mesh contains no cells");
    mesh->cell.i_node = memory_calloc(mesh->n_cells + 1, sizeof(*mesh->cell.i_node));
    mesh->cell.node = memory_calloc(n_node_tags[0] + n_node_tags[1], sizeof(*mesh->cell.node));

    // create cell to node map
    long n = 0;
    for (long d = N_DIMS; d >= N_DIMS - 1; --d) {
        size_t dimTags_n;
        int *dimTags;
        gmshModelGetPhysicalGroups(&dimTags, &dimTags_n, d, &ierr);
        for (size_t p = 0; p < dimTags_n; p += 2) {
            size_t tags_n;
            int *tags;
            gmshModelGetEntitiesForPhysicalGroup(dimTags[p], dimTags[p + 1], &tags, &tags_n, &ierr);
            for (size_t e = 0; e < tags_n; ++e) {
                size_t elementTypes_n, *elementTags_n, *nodeTags_n, **nodeTags;
                int *elementTypes;
                gmshModelMeshGetElements(&elementTypes, &elementTypes_n, 0, &elementTags_n, 0,
                                         &nodeTags, &nodeTags_n, 0, d, tags[e], &ierr);
                for (size_t k = 0; k < elementTypes_n; ++k) {
                    long n_nodes;
                    switch (elementTypes[k]) {
                        case 1:  // 2-node line
                            n_nodes = 2;
                            break;
                        case 2:  // 3-node triangle
                            n_nodes = 3;
                            break;
                        case 3:  // 4-node quadrangle
                            n_nodes = 4;
                            break;
                        default:
                            assert(!"unsupported mesh dimension");
                    }
                    for (size_t j = 0; j < elementTags_n[k]; ++j) {
                        mesh->cell.i_node[n + 1] = mesh->cell.i_node[n] + n_nodes;
                        for (long i = 0; i < n_nodes; ++i)
                            mesh->cell.node[mesh->cell.i_node[n] + i] =
                                node_map[nodeTags[k][j * n_nodes + i]];
                        n += 1;
                    }
                }
                gmshFree(elementTags_n);
                gmshFree(nodeTags_n);
                gmshFree(elementTypes);
                for (size_t i = 0; i < elementTypes_n; ++i) gmshFree(nodeTags[i]);
                gmshFree(nodeTags);
            }
            gmshFree(tags);
        }
        gmshFree(dimTags);
    }
    assert(n == mesh->n_cells && "cell creation error");
}

static void create_entities(Mesh *mesh)
{
    int ierr;

    // count number of mesh entites
    for (long d = N_DIMS; d >= N_DIMS - 1; --d) {
        size_t dimTags_n;
        gmshModelGetPhysicalGroups(0, &dimTags_n, d, &ierr);
        mesh->n_entities += dimTags_n / 2;
    }
    assert(mesh->n_entities > 0 && "mesh contains no entities");

    mesh->entity.name = memory_calloc(mesh->n_entities, sizeof(*mesh->entity.name));
    mesh->entity.j_cell = memory_calloc(mesh->n_entities + 1, sizeof(*mesh->entity.j_cell));

    // create mesh entities
    long n = 0;
    for (long d = N_DIMS; d >= N_DIMS - 1; --d) {
        size_t dimTags_n;
        int *dimTags;
        gmshModelGetPhysicalGroups(&dimTags, &dimTags_n, d, &ierr);
        for (size_t p = 0; p < dimTags_n; p += 2) {
            gmshModelGetPhysicalName(dimTags[p], dimTags[p + 1], &mesh->entity.name[n], &ierr);
            if (!strcmp(mesh->entity.name[n], "")) {
                char buf[128];
                snprintf(buf, sizeof(buf), "%d", dimTags[p + 1]);
                mesh->entity.name[n] = utils_strdup(buf);
            }
            mesh->entity.j_cell[n + 1] = mesh->entity.j_cell[n];

            size_t tags_n;
            int *tags;
            gmshModelGetEntitiesForPhysicalGroup(dimTags[p], dimTags[p + 1], &tags, &tags_n, &ierr);
            for (size_t e = 0; e < tags_n; ++e) {
                size_t elementTypes_n, *elementTags_n;
                gmshModelMeshGetElements(0, &elementTypes_n, 0, &elementTags_n, 0, 0, 0, 0, d,
                                         tags[e], &ierr);
                for (size_t i = 0; i < elementTypes_n; ++i) {
                    mesh->entity.j_cell[n + 1] += elementTags_n[i];
                }
                gmshFree(elementTags_n);
            }
            gmshFree(tags);
            n += 1;
        }
        gmshFree(dimTags);
    }
    assert(n == mesh->n_entities && "entity creation error");
}
