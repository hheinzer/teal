#include "mesh.h"

#include <assert.h>
#include <gmshc.h>
#include <hdf5.h>
#include <math.h>
#include <metis.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "array.h"
#include "dict.h"
#include "global.h"
#include "gmsh_io.h"
#include "hdf5_io.h"
#include "memory.h"
#include "partition.h"
#include "sync.h"
#include "utils.h"

static void compute_cell_to_cell(Mesh *mesh);
static void compute_face_connectivity(Mesh *mesh);
static void compute_face_areas(Mesh *mesh);
static void compute_cell_volumes(Mesh *mesh);
static void compute_face_centers(Mesh *mesh);
static void compute_cell_centers(Mesh *mesh);
static void compute_face_normals(Mesh *mesh);
static void compute_cell_projections(Mesh *mesh);
static void compute_face_gauss_weights(Mesh *mesh);
static void compute_cell_gauss_weights(Mesh *mesh);
static void compute_face_gradient_weights(Mesh *mesh);
static void compute_face_reconstruction(Mesh *mesh);
static void compute_cell_reconstruction(Mesh *mesh);
static void compute_entity_to_face(Mesh *mesh);
static long common_cell_nodes(const Mesh *mesh, const long cell1, const long cell2, long *node);
static double volume(const double (*x)[N_DIMS], const long *node, const long n_nodes,
                     const long dim);
static void center(double *center, const double (*x)[N_DIMS], const long *node, const long n_nodes,
                   const long dim);
static void normal(double *normal, const double (*x)[N_DIMS], const long *node, const long n_nodes);
static double *cross(const double *a, const double *b, double *c);

Mesh mesh_create(const double *x1, const double *x2, const long *n_cells) {
    for (long d = 0; d < N_DIMS; ++d) {
        assert(x1[d] < x2[d] && "illegal domain extent");
        assert(n_cells[d] > 0 && "illegal number of cells");
    }

    Mesh mesh = {};

    // allocate memory for nodes
    long n_nodes[N_DIMS];
    for (long d = 0; d < N_DIMS; ++d) n_nodes[d] = n_cells[d] + 1;
    mesh.n_nodes = array_product(n_nodes, N_DIMS);
    mesh.node.x = memory_calloc(mesh.n_nodes, sizeof(*mesh.node.x));

    // compute node coordinates
    for (long n = 0, j = 0; j < n_nodes[1]; ++j) {
        for (long i = 0; i < n_nodes[0]; ++i) {
            mesh.node.x[n][X] = x1[X] + (x2[X] - x1[X]) * i / n_cells[X];
            mesh.node.x[n][Y] = x1[Y] + (x2[Y] - x1[Y]) * j / n_cells[Y];
            n += 1;
        }
    }

    // allocate memory for cell to node info
    mesh.n_inner_cells = array_product(n_cells, N_DIMS);
    mesh.n_ghost_cells = 2 * array_sum(n_cells, N_DIMS);
    mesh.n_cells = mesh.n_inner_cells + mesh.n_ghost_cells;
    mesh.cell.i_node = memory_calloc(mesh.n_cells + 1, sizeof(*mesh.cell.i_node));
    mesh.cell.node =
        memory_calloc(4 * mesh.n_inner_cells + 2 * mesh.n_ghost_cells, sizeof(*mesh.cell.node));

    // compute inner cell to node info
    long n = 0;
    for (long j = 0; j < n_cells[Y]; ++j) {
        for (long i = 0; i < n_cells[X]; ++i) {
            mesh.cell.i_node[n + 1] = mesh.cell.i_node[n] + 4;
            mesh.cell.node[mesh.cell.i_node[n] + 0] = (i + 0) + (j + 0) * n_nodes[X];
            mesh.cell.node[mesh.cell.i_node[n] + 1] = (i + 1) + (j + 0) * n_nodes[X];
            mesh.cell.node[mesh.cell.i_node[n] + 2] = (i + 1) + (j + 1) * n_nodes[X];
            mesh.cell.node[mesh.cell.i_node[n] + 3] = (i + 0) + (j + 1) * n_nodes[X];
            n += 1;
        }
    }

    // compute ghost cell to node info
    for (long j = 0; j < 4; ++j) {  // bottom, right, top, left
        for (long i = 0; i < n_cells[j % 2]; ++i) {
            mesh.cell.i_node[n + 1] = mesh.cell.i_node[n] + 2;
            switch (j) {
                case 0:
                    mesh.cell.node[mesh.cell.i_node[n] + 0] = (i + 0) + (0 + 0) * n_nodes[X];
                    mesh.cell.node[mesh.cell.i_node[n] + 1] = (i + 1) + (0 + 0) * n_nodes[X];
                    break;
                case 1:
                    mesh.cell.node[mesh.cell.i_node[n] + 0] = n_cells[X] + (i + 0) * n_nodes[X];
                    mesh.cell.node[mesh.cell.i_node[n] + 1] = n_cells[X] + (i + 1) * n_nodes[X];
                    break;
                case 2:
                    mesh.cell.node[mesh.cell.i_node[n] + 0] = (i + 1) + n_cells[Y] * n_nodes[X];
                    mesh.cell.node[mesh.cell.i_node[n] + 1] = (i + 0) + n_cells[Y] * n_nodes[X];
                    break;
                case 3:
                    mesh.cell.node[mesh.cell.i_node[n] + 0] = (0 + 0) + (i + 1) * n_nodes[X];
                    mesh.cell.node[mesh.cell.i_node[n] + 1] = (0 + 0) + (i + 0) * n_nodes[X];
                    break;
                default:
                    break;
            }
            n += 1;
        }
    }

    mesh.n_entities = 5;
    mesh.entity.name = memory_calloc(mesh.n_entities, sizeof(*mesh.entity.name));
    mesh.entity.j_cell = memory_calloc(mesh.n_entities + 1, sizeof(*mesh.entity.j_cell));

    mesh.entity.name[0] = utils_strdup("domain");
    mesh.entity.name[1] = utils_strdup("bottom");
    mesh.entity.name[2] = utils_strdup("right");
    mesh.entity.name[3] = utils_strdup("top");
    mesh.entity.name[4] = utils_strdup("left");

    // compute entity cell ranges
    mesh.entity.j_cell[1] = mesh.n_inner_cells;
    for (long e = 1; e < mesh.n_entities; ++e)
        mesh.entity.j_cell[e + 1] = mesh.entity.j_cell[e] + n_cells[(e - 1) % 2];

    mesh_finalize(&mesh);
    return mesh;
}

Mesh mesh_read(const char *fname) {
    Mesh mesh = {};
    const char *ext = utils_extension(fname);
    if (!strcmp(ext, "geo") || !strcmp(ext, "msh")) {
        gmsh_read(&mesh, fname);
    } else {
        assert(!"unsupported mesh format");
    }
    mesh_finalize(&mesh);
    return mesh;
}

void mesh_finalize(Mesh *mesh) {
    compute_cell_to_cell(mesh);
    mesh->entity.bc.name = memory_calloc(mesh->n_entities, sizeof(*mesh->entity.bc.name));
    mesh->entity.bc.context = memory_calloc(mesh->n_entities, sizeof(*mesh->entity.bc.context));
    mesh->entity.bc.apply = memory_calloc(mesh->n_entities, sizeof(*mesh->entity.bc.apply));
}

void mesh_free(Mesh *mesh) {
    free(mesh->node.x);
    free(mesh->node.map);

    free(mesh->cell.i_node);
    free(mesh->cell.node);
    free(mesh->cell.i_cell);
    free(mesh->cell.cell);
    free(mesh->cell.volume);
    free(mesh->cell.center);
    free(mesh->cell.projection);
    free(mesh->cell.gauss_weight);
    free(mesh->cell.reconstruction);

    free(mesh->face.i_node);
    free(mesh->face.node);
    free(mesh->face.cell);
    free(mesh->face.area);
    free(mesh->face.center);
    free(mesh->face.normal);
    free(mesh->face.gauss_weight);
    free(mesh->face.gradient_weight);
    free(mesh->face.reconstruction);

    for (long e = 0; e < mesh->n_entities; ++e) free(mesh->entity.name[e]);
    free(mesh->entity.name);
    free(mesh->entity.j_cell);
    free(mesh->entity.j_face);
    for (long e = 0; e < mesh->n_entities; ++e) free(mesh->entity.bc.name[e]);
    free(mesh->entity.bc.name);
    free(mesh->entity.bc.context);
    free(mesh->entity.bc.apply);

    free(mesh->sync.i_recv);
    free(mesh->sync.i_send);
    free(mesh->sync.send);

    *mesh = (typeof(*mesh)){};
}

void mesh_modify(Mesh *mesh, Modify *modify) {
    for (long i = 0; i < mesh->n_nodes; ++i) modify(mesh->node.x[i]);
}

void mesh_split_entity(Mesh *mesh, const char *name, const long dim, const double x) {
    assert(dim < N_DIMS && "unsupported dimension");

    // find entity
    long E = -1;
    for (long e = 0; e < mesh->n_entities; ++e)
        if (!strcmp(mesh->entity.name[e], name)) E = e;
    assert(E != -1 && "entity not found");

    // compute cell split
    const long n_cells = mesh->entity.j_cell[E + 1] - mesh->entity.j_cell[E];
    cleanup long *split = memory_calloc(n_cells, sizeof(*split));
    const ALIAS(xn, mesh->node.x);
    for (long s = 0, j = mesh->entity.j_cell[E]; j < mesh->entity.j_cell[E + 1]; ++j) {
        const long n_nodes = mesh->cell.i_node[j + 1] - mesh->cell.i_node[j];
        double mean[N_DIMS] = {};
        for (long i = mesh->cell.i_node[j]; i < mesh->cell.i_node[j + 1]; ++i) {
            for (long d = 0; d < N_DIMS; ++d) mean[d] += xn[mesh->cell.node[i]][d] / n_nodes;
        }
        split[s++] = (mean[dim] > x);
    }

    // reorder cell to node
    long *i_node = memory_calloc(mesh->n_cells + 1, sizeof(*i_node));
    long *node = memory_calloc(mesh->cell.i_node[mesh->n_cells], sizeof(*node));
    for (long n = 0, e = 0; e < mesh->n_entities; ++e) {
        if (e != E) {
            for (long j = mesh->entity.j_cell[e]; j < mesh->entity.j_cell[e + 1]; ++j) {
                i_node[n + 1] = i_node[n];
                for (long i = mesh->cell.i_node[j]; i < mesh->cell.i_node[j + 1]; ++i)
                    node[i_node[n + 1]++] = mesh->cell.node[i];
                n += 1;
            }
        } else {
            for (long k = 0; k < 2; ++k) {
                for (long s = 0, j = mesh->entity.j_cell[e]; j < mesh->entity.j_cell[e + 1]; ++j) {
                    if (split[s++] != k) continue;
                    i_node[n + 1] = i_node[n];
                    for (long i = mesh->cell.i_node[j]; i < mesh->cell.i_node[j + 1]; ++i)
                        node[i_node[n + 1]++] = mesh->cell.node[i];
                    n += 1;
                }
            }
        }
    }
    free(mesh->cell.i_node);
    free(mesh->cell.node);
    mesh->cell.i_node = i_node;
    mesh->cell.node = node;

    // reorder cell to cell
    long *i_cell = memory_calloc(mesh->n_cells + 1, sizeof(*i_cell));
    long *cell = memory_calloc(mesh->cell.i_cell[mesh->n_cells], sizeof(*cell));
    for (long n = 0, e = 0; e < mesh->n_entities; ++e) {
        if (e != E) {
            for (long j = mesh->entity.j_cell[e]; j < mesh->entity.j_cell[e + 1]; ++j) {
                i_cell[n + 1] = i_cell[n];
                for (long i = mesh->cell.i_cell[j]; i < mesh->cell.i_cell[j + 1]; ++i)
                    cell[i_cell[n + 1]++] = mesh->cell.cell[i];
                n += 1;
            }
        } else {
            for (long k = 0; k < 2; ++k) {
                for (long s = 0, j = mesh->entity.j_cell[e]; j < mesh->entity.j_cell[e + 1]; ++j) {
                    if (split[s++] != k) continue;
                    i_cell[n + 1] = i_cell[n];
                    for (long i = mesh->cell.i_cell[j]; i < mesh->cell.i_cell[j + 1]; ++i)
                        cell[i_cell[n + 1]++] = mesh->cell.cell[i];
                    n += 1;
                }
            }
        }
    }
    free(mesh->cell.i_cell);
    free(mesh->cell.cell);
    mesh->cell.i_cell = i_cell;
    mesh->cell.cell = cell;

    // split entity
    const long n_entities = mesh->n_entities + 1;
    char **ename = memory_calloc(n_entities + 1, sizeof(*ename));
    long *j_cell = memory_calloc(n_entities + 1, sizeof(*j_cell));
    for (long n = 0, e = 0; e < mesh->n_entities; ++e) {
        if (e != E) {
            ename[n] = utils_strdup(mesh->entity.name[e]);
            j_cell[n + 1] = j_cell[n] + mesh->entity.j_cell[e + 1] - mesh->entity.j_cell[e];
            n += 1;
        } else {
            for (long k = 0; k < 2; ++k) {
                char buf[256];
                snprintf(buf, sizeof(buf), "%s%ld", mesh->entity.name[e], k);
                ename[n] = utils_strdup(buf);
                j_cell[n + 1] = j_cell[n];
                for (long s = 0; s < n_cells; ++s)
                    if (split[s] == k) j_cell[n + 1] += 1;
                n += 1;
            }
        }
    }
    for (long e = 0; e < mesh->n_entities; ++e) free(mesh->entity.name[e]);
    free(mesh->entity.name);
    free(mesh->entity.j_cell);
    mesh->entity.name = ename;
    mesh->entity.j_cell = j_cell;

    free(mesh->entity.bc.name);
    free(mesh->entity.bc.context);
    free(mesh->entity.bc.apply);
    mesh->entity.bc.name = memory_calloc(n_entities, sizeof(*mesh->entity.bc.name));
    mesh->entity.bc.context = memory_calloc(n_entities, sizeof(*mesh->entity.bc.context));
    mesh->entity.bc.apply = memory_calloc(n_entities, sizeof(*mesh->entity.bc.apply));

    mesh->n_entities = n_entities;
}

void mesh_set_periodic_condition(Mesh *mesh, const char *name0, const char *name1) {
    // find entities
    long E[2] = {-1, -1};
    for (long e = 0; e < mesh->n_entities; ++e) {
        if (!strcmp(mesh->entity.name[e], name0)) E[0] = e;
        if (!strcmp(mesh->entity.name[e], name1)) E[1] = e;
    }
    assert(E[0] != -1 && E[1] != -1 && "entity not found");

    // compute number of periodic cells
    const long np = mesh->entity.j_cell[E[0] + 1] - mesh->entity.j_cell[E[0]];
    const long n1 = mesh->entity.j_cell[E[1] + 1] - mesh->entity.j_cell[E[1]];
    assert(np == n1 && "inconsistent number of faces");

    // compute cell coordinates
    cleanup double(*x)[np][N_DIMS] = memory_calloc(2, sizeof(*x));
    for (long e = 0; e < 2; ++e) {
        for (long n = 0, j = mesh->entity.j_cell[E[e]]; j < mesh->entity.j_cell[E[e] + 1]; ++j) {
            long m = 0;
            for (long i = mesh->cell.i_node[j]; i < mesh->cell.i_node[j + 1]; ++i) {
                for (long d = 0; d < N_DIMS; ++d) {
                    x[e][n][d] += mesh->node.x[mesh->cell.node[i]][d];
                }
                m += 1;
            }
            for (long d = 0; d < N_DIMS; ++d) x[e][n][d] /= m;
            n += 1;
        }
    }

    // compute mean coordinates
    double mean[2][N_DIMS] = {};
    for (long e = 0; e < 2; ++e)
        for (long i = 0; i < np; ++i)
            for (long d = 0; d < N_DIMS; ++d) mean[e][d] += x[e][i][d] / np;

    // compute coordinate to cell dict
    fcleanup(dict_free) Dict x2c = dict_create(2 * np);
    for (long e = 0; e < 2; ++e) {
        for (long i = 0; i < np; ++i) {
            long key[N_DIMS], val = mesh->entity.j_cell[E[e]] + i;
            for (long d = 0; d < N_DIMS; ++d) key[d] = round(x[e][i][d] / EPS);
            dict_insert(&x2c, key, N_DIMS, &val, 1);
        }
    }
    assert(x2c.n_items == 2 * np && "coordinate resolution error");

    // compute cell to cell dict
    fcleanup(dict_free) Dict c2c = dict_create(2 * np);
    for (long e = 0; e < 2; ++e) {
        for (long i = 0; i < np; ++i) {
            long key[N_DIMS], *val;
            for (long d = 0; d < N_DIMS; ++d)
                key[d] = round((x[e][i][d] - mean[e][d] + mean[(e + 1) % 2][d]) / EPS);
            dict_lookup(&x2c, key, N_DIMS, &val);

            const long ghost = mesh->entity.j_cell[E[e]] + i;
            const long inner = mesh->cell.cell[mesh->cell.i_cell[*val]];
            dict_insert(&c2c, &ghost, 1, &inner, 1);
        }
    }
    assert(c2c.n_items == 2 * np && "cell to cell mapping error");

    // modify connectivity of inner cells
    for (long e = 0; e < 2; ++e) {
        for (long jj = mesh->entity.j_cell[E[e]]; jj < mesh->entity.j_cell[E[e] + 1]; ++jj) {
            const long j = mesh->cell.cell[mesh->cell.i_cell[jj]];
            for (long i = mesh->cell.i_cell[j]; i < mesh->cell.i_cell[j + 1]; ++i) {
                long *val;
                if (dict_lookup(&c2c, &mesh->cell.cell[i], 1, &val)) mesh->cell.cell[i] = *val;
            }
        }
    }

    // compute periodic cell to cell map
    long *i_cell = memory_calloc(mesh->n_cells + 1, sizeof(*i_cell));
    long *cell = memory_calloc(mesh->n_cells * MAX_CELL_FACES, sizeof(*cell));
    for (long n = 0, e = 0; e < mesh->n_entities; ++e) {
        for (long j = mesh->entity.j_cell[e]; j < mesh->entity.j_cell[e + 1]; ++j) {
            i_cell[n + 1] = i_cell[n];
            for (long i = mesh->cell.i_cell[j]; i < mesh->cell.i_cell[j + 1]; ++i) {
                if (j == mesh->cell.cell[i]) continue;  // filter out self connections
                cell[i_cell[n + 1]++] = mesh->cell.cell[i];
            }
            if (e == E[0] || e == E[1]) {
                long *val;
                dict_lookup(&c2c, &j, 1, &val);
                cell[i_cell[n + 1]++] = *val;
            }
            n += 1;
        }
    }
    free(mesh->cell.i_cell);
    free(mesh->cell.cell);
    mesh->cell.i_cell = i_cell;
    mesh->cell.cell = memory_realloc(cell, i_cell[mesh->n_cells], sizeof(*cell));

    // mark entities
    for (long e = 0; e < 2; ++e) mesh->entity.bc.name[E[e]] = utils_strdup("periodic");
}

void mesh_set_boundary_condition(Mesh *mesh, const char *name, const char *bc, const double *state,
                                 Function *custom) {
    for (long e = 0; e < mesh->n_entities; ++e) {
        if (!strcmp(mesh->entity.name[e], name)) {
            mesh->entity.bc.name[e] = utils_strdup(bc);
            mesh->entity.bc.context[e] = (BCContext){state, custom};
            return;
        }
    }
    assert(!"entity not found");
}

void mesh_generate(Mesh *mesh) {
    partition(mesh, true);
    compute_face_connectivity(mesh);
    compute_face_areas(mesh);
    compute_cell_volumes(mesh);
    compute_face_centers(mesh);
    compute_cell_centers(mesh);
    compute_face_normals(mesh);
    compute_cell_projections(mesh);
    compute_face_gauss_weights(mesh);
    compute_cell_gauss_weights(mesh);
    compute_face_gradient_weights(mesh);
    compute_face_reconstruction(mesh);
    compute_cell_reconstruction(mesh);
    compute_entity_to_face(mesh);
    mesh->volume = sync_sum(array_sum(mesh->cell.volume, mesh->n_inner_cells));
}

void mesh_print(const Mesh *mesh) {
    const long n_inner_nodes = sync_sum(mesh->n_inner_nodes);

    double node_min[N_DIMS], node_max[N_DIMS];
    sync_min(array_min_s(*mesh->node.x, mesh->n_inner_nodes, N_DIMS, node_min), N_DIMS);
    sync_max(array_max_s(*mesh->node.x, mesh->n_inner_nodes, N_DIMS, node_max), N_DIMS);

    const long n_inner_cells = sync_sum(mesh->n_inner_cells);
    const long n_ghost_cells = sync_sum(mesh->n_ghost_cells);

    const long n_inner_faces = sync_sum(mesh->n_inner_faces);
    const long n_bound_faces = sync_sum(mesh->n_bound_faces);
    const long n_faces = sync_sum(mesh->n_faces);

    long n_entity_cells[mesh->n_entities];
    for (long e = 0; e < mesh->n_entities; ++e)
        n_entity_cells[e] = mesh->entity.j_cell[e + 1] - mesh->entity.j_cell[e];
    sync_sum(n_entity_cells, mesh->n_entities);

    if (mesh->rank == 0) {
        printf("Mesh summary:\n");

        printf(" | " FMT_KEY ": %ld\n", "number of nodes", n_inner_nodes);

        printf(" | " FMT_KEY ": ", "min/max node");
        array_print(0, node_min, N_DIMS, " / ");
        array_print(0, node_max, N_DIMS, "\n");

        printf(" | " FMT_KEY ": %ld\n", "number of inner cells", n_inner_cells);
        printf(" | " FMT_KEY ": %ld\n", "number of ghost cells", n_ghost_cells);

        printf(" | " FMT_KEY ": %g\n", "total volume", mesh->volume);

        printf(" | " FMT_KEY ": %ld\n", "number of faces",
               (n_inner_faces + n_bound_faces + n_faces) / 2);

        char buf[256] = "";
        for (long e = 0; e < mesh->n_entities; ++e) {
            snprintf(buf, sizeof(buf), "entity %ld", e + 1);
            printf(" | " FMT_KEY ": %s(%ld) %s\n", buf, mesh->entity.name[e], n_entity_cells[e],
                   (mesh->entity.bc.name[e] ? mesh->entity.bc.name[e] : ""));
        }
    }
}

void mesh_write(const Mesh *mesh, const char *prefix) {
    char fname[256];
    snprintf(fname, sizeof(fname), "%s_mesh.vtkhdf", prefix);

    const long n_inner_nodes = mesh->n_inner_nodes;
    const long n_inner_cells = mesh->n_inner_cells;
    const int rank = mesh->rank;

    const long n_points = sync_sum(n_inner_nodes);
    const long n_local_conns = mesh->cell.i_node[n_inner_cells];
    const long n_conns = sync_sum(n_local_conns);
    const long conn_offset = sync_exscan_sum(n_local_conns);
    const long n_cells = sync_sum(n_inner_cells);

    // https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html#unstructured-grid
    hid_t file = hdf5_file_create(fname);
    hid_t vtkhdf = hdf5_group_create(file, "VTKHDF");

    const long version[] = {2, 0};
    hdf5_write_attribute(vtkhdf, "Version", version, 1, HDF5_DIMS(2));
    hdf5_write_attribute(vtkhdf, "Type", "UnstructuredGrid");

    hdf5_write_dataset(vtkhdf, "NumberOfPoints", &n_points, 1, HDF5_DIMS(rank == 0));
    hdf5_write_dataset(vtkhdf, "NumberOfConnectivityIds", &n_conns, 1, HDF5_DIMS(rank == 0));
    hdf5_write_dataset(vtkhdf, "NumberOfCells", &n_cells, 1, HDF5_DIMS(rank == 0));

    cleanup double(*points)[3] = memory_calloc(n_inner_nodes, sizeof(*points));
    for (long i = 0; i < n_inner_nodes; ++i)
        for (long d = 0; d < N_DIMS; ++d) points[i][d] = mesh->node.x[i][d];
    hdf5_write_dataset(vtkhdf, "Points", *points, 2, HDF5_DIMS(n_inner_nodes, 3));

    cleanup long *conn = memory_calloc(n_local_conns, sizeof(*conn));
    for (long i = 0; i < n_local_conns; ++i) conn[i] = mesh->node.map[mesh->cell.node[i]];
    hdf5_write_dataset(vtkhdf, "Connectivity", conn, 1, HDF5_DIMS(n_local_conns));

    const long n_offsets = n_inner_cells + (rank == 0);
    cleanup long *offsets = memory_calloc(n_offsets, sizeof(*offsets));
    for (long i = 0; i < n_offsets; ++i)
        offsets[i] = mesh->cell.i_node[(rank != 0) + i] + conn_offset;
    hdf5_write_dataset(vtkhdf, "Offsets", offsets, 1, HDF5_DIMS(n_offsets));

    cleanup unsigned char *types = memory_calloc(n_inner_cells, sizeof(*types));
    for (long i = 0; i < n_inner_cells; ++i) {
        switch (mesh->cell.i_node[i + 1] - mesh->cell.i_node[i]) {
            case 3:
                types[i] = 5;  // VTK_TRIANGLE
                break;
            case 4:
                types[i] = 9;  // VTK_QUAD
                break;
            default:
                assert(!"unsupported cell type");
        }
    }
    hdf5_write_dataset(vtkhdf, "Types", types, 1, HDF5_DIMS(n_inner_cells));

    hid_t cell_data = hdf5_group_create(vtkhdf, "CellData");
    cleanup long *buf = memory_calloc(n_inner_cells, sizeof(*buf));

    for (long i = 0; i < n_inner_cells; ++i) buf[i] = i;
    hdf5_write_dataset(cell_data, "index", buf, 1, HDF5_DIMS(n_inner_cells));

    for (long i = 0; i < n_inner_cells; ++i) buf[i] = rank;
    hdf5_write_dataset(cell_data, "rank", buf, 1, HDF5_DIMS(n_inner_cells));

    hdf5_write_dataset(cell_data, "volume", mesh->cell.volume, 1, HDF5_DIMS(n_inner_cells));

    hdf5_write_dataset(cell_data, "center", *mesh->cell.center, 2,
                       HDF5_DIMS(n_inner_cells, N_DIMS));

    hdf5_group_close(cell_data);
    hdf5_group_close(vtkhdf);

    hdf5_file_close(file);
}

static void compute_cell_to_cell(Mesh *mesh) {
    // create node to cell dict
    fcleanup(dict_free) Dict n2c = dict_create(mesh->n_nodes);
    for (long j = 0; j < mesh->n_cells; ++j)
        for (long i = mesh->cell.i_node[j]; i < mesh->cell.i_node[j + 1]; ++i)
            dict_append(&n2c, &mesh->cell.node[i], 1, &j, 1);

    // create cell to cell dict
    fcleanup(dict_free) Dict c2c_dict = dict_create(mesh->n_cells);
    for (long j1 = 0; j1 < mesh->n_cells; ++j1) {
        for (long i = mesh->cell.i_node[j1]; i < mesh->cell.i_node[j1 + 1]; ++i) {
            long *c2, nc2 = dict_lookup(&n2c, &mesh->cell.node[i], 1, &c2);
            for (const long *j2 = c2; j2 < c2 + nc2; ++j2) {
                if (*j2 == j1) continue;
                if (common_cell_nodes(mesh, j1, *j2, 0) >= N_DIMS) {
                    long *c1, nc1 = dict_lookup(&c2c_dict, &j1, 1, &c1);
                    if (nc1 && array_contains(c1, nc1, *j2)) continue;
                    dict_append(&c2c_dict, &j1, 1, j2, 1);
                }
            }
        }
    }
    cleanup DictItem *c2c = dict_serialize(&c2c_dict);

    // create cell to cell map
    long *i_cell = memory_calloc(mesh->n_cells + 1, sizeof(*i_cell));
    long *cell = memory_calloc(mesh->n_cells * MAX_CELL_FACES, sizeof(*cell));
    for (long j = 0; j < mesh->n_cells; ++j) {
        const long *c = c2c[j].val;
        const long nc = c2c[j].n_val;
        i_cell[j + 1] = i_cell[j] + nc;
        for (long i = 0; i < nc; ++i) cell[i_cell[j] + i] = c[i];
    }
    mesh->cell.i_cell = i_cell;
    mesh->cell.cell = memory_realloc(cell, i_cell[mesh->n_cells], sizeof(*cell));
}

static void compute_face_connectivity(Mesh *mesh) {
    // find periodic connections
    fcleanup(dict_free) Dict periodic = dict_create(mesh->n_ghost_cells);
    for (long e = 0; e < mesh->n_entities; ++e) {
        if (!mesh->entity.bc.name[e] || strcmp(mesh->entity.bc.name[e], "periodic")) continue;
        for (long j = mesh->entity.j_cell[e]; j < mesh->entity.j_cell[e + 1]; ++j) {
            const long *cell = &mesh->cell.cell[mesh->cell.i_cell[j]];
            const long *node = &mesh->cell.node[mesh->cell.i_node[j]];
            const long n_nodes = mesh->cell.i_node[j + 1] - mesh->cell.i_node[j];
            dict_insert(&periodic, cell, 2, node, n_nodes);
        }
    }

    // find all inner faces
    fcleanup(dict_free) Dict f2n_dict = dict_create(mesh->n_cells * MAX_CELL_FACES);
    for (long j = 0; j < mesh->n_inner_cells; ++j) {
        for (long i = mesh->cell.i_cell[j]; i < mesh->cell.i_cell[j + 1]; ++i) {
            if (mesh->cell.cell[i] >= mesh->n_inner_cells) continue;
            const long key[2] = {MIN(j, mesh->cell.cell[i]), MAX(j, mesh->cell.cell[i])};
            long *node;
            if (!dict_lookup(&f2n_dict, key, 2, &node)) {
                long cnode[MAX_FACE_NODES];
                const long n_cnodes = common_cell_nodes(mesh, j, mesh->cell.cell[i], cnode);
                if (n_cnodes >= N_DIMS) {
                    dict_insert(&f2n_dict, key, 2, cnode, n_cnodes);
                } else {
                    long *pnode, n_pnodes = dict_lookup(&periodic, key, 2, &pnode);
                    assert(n_pnodes >= N_DIMS && "illegal number of face nodes");
                    dict_insert(&f2n_dict, key, 2, pnode, n_pnodes);
                }
            }
        }
    }
    mesh->n_inner_faces = f2n_dict.n_items;

    // find all boundary faces
    for (long e = 0; e < mesh->n_entities; ++e) {
        if (!mesh->entity.bc.name[e] || !strcmp(mesh->entity.bc.name[e], "periodic")) continue;
        for (long j = mesh->entity.j_cell[e]; j < mesh->entity.j_cell[e + 1]; ++j) {
            for (long i = mesh->cell.i_cell[j]; i < mesh->cell.i_cell[j + 1]; ++i) {
                const long key[2] = {j, mesh->cell.cell[i]};
                long node[MAX_FACE_NODES];
                const long n_nodes = common_cell_nodes(mesh, key[0], key[1], node);
                assert(n_nodes >= N_DIMS && "illegal number of face nodes");
                dict_insert(&f2n_dict, key, 2, node, n_nodes);
            }
        }
    }
    mesh->n_bound_faces = f2n_dict.n_items - mesh->n_inner_faces;

    // find all mpi faces
    for (long j = mesh->n_inner_cells + mesh->n_ghost_cells; j < mesh->n_cells; ++j) {
        for (long i = mesh->cell.i_cell[j]; i < mesh->cell.i_cell[j + 1]; ++i) {
            const long key[2] = {j, mesh->cell.cell[i]};
            long node[MAX_FACE_NODES];
            const long n_nodes = common_cell_nodes(mesh, key[0], key[1], node);
            if (n_nodes >= N_DIMS) {
                dict_insert(&f2n_dict, key, 2, node, n_nodes);
            } else {
                long *pnode, n_pnodes = dict_lookup(&periodic, key, 2, &pnode);
                assert(n_pnodes >= N_DIMS && "illegal number of face nodes");
                dict_insert(&f2n_dict, key, 2, pnode, n_pnodes);
            }
        }
    }
    mesh->n_faces = f2n_dict.n_items;

    cleanup DictItem *f2n = dict_serialize(&f2n_dict);
    long *i_node = memory_calloc(mesh->n_faces + 1, sizeof(*i_node));
    long *node = memory_calloc(mesh->n_faces * MAX_FACE_NODES, sizeof(*node));
    long(*cell)[2] = memory_calloc(mesh->n_faces, sizeof(*cell));
    for (long i = 0; i < mesh->n_faces; ++i) {
        i_node[i + 1] = i_node[i];
        for (long j = 0; j < f2n[i].n_val; ++j) node[i_node[i + 1]++] = f2n[i].val[j];
        cell[i][0] = MIN(f2n[i].key[0], f2n[i].key[1]);
        cell[i][1] = MAX(f2n[i].key[0], f2n[i].key[1]);
    }
    mesh->face.i_node = i_node;
    mesh->face.node = memory_realloc(node, i_node[mesh->n_faces], sizeof(*node));
    mesh->face.cell = cell;
}

static void compute_face_areas(Mesh *mesh) {
    const ALIAS(x, mesh->node.x);
    const ALIAS(node, mesh->face.node);
    const ALIAS(i_node, mesh->face.i_node);
    mesh->face.area = memory_calloc(mesh->n_faces, sizeof(*mesh->face.area));
    for (long i = 0; i < mesh->n_faces; ++i)
        mesh->face.area[i] = volume(x, &node[i_node[i]], i_node[i + 1] - i_node[i], N_DIMS - 1);
}

static void compute_cell_volumes(Mesh *mesh) {
    const ALIAS(x, mesh->node.x);
    const ALIAS(node, mesh->cell.node);
    const ALIAS(i_node, mesh->cell.i_node);
    mesh->cell.volume = memory_calloc(mesh->n_inner_cells, sizeof(*mesh->cell.volume));
    for (long i = 0; i < mesh->n_inner_cells; ++i)
        mesh->cell.volume[i] = volume(x, &node[i_node[i]], i_node[i + 1] - i_node[i], N_DIMS);
}

static void compute_face_centers(Mesh *mesh) {
    const ALIAS(x, mesh->node.x);
    const ALIAS(node, mesh->face.node);
    const ALIAS(i_node, mesh->face.i_node);
    mesh->face.center = memory_calloc(mesh->n_faces, sizeof(*mesh->face.center));
    for (long i = 0; i < mesh->n_faces; ++i)
        center(mesh->face.center[i], x, &node[i_node[i]], i_node[i + 1] - i_node[i], N_DIMS - 1);
}

static void compute_cell_centers(Mesh *mesh) {
    const ALIAS(x, mesh->node.x);
    const ALIAS(node, mesh->cell.node);
    const ALIAS(i_node, mesh->cell.i_node);
    mesh->cell.center = memory_calloc(mesh->n_cells, sizeof(*mesh->cell.center));

    for (long i = 0; i < mesh->n_inner_cells; ++i)
        center(mesh->cell.center[i], x, &node[i_node[i]], i_node[i + 1] - i_node[i], N_DIMS);

    // mirror inner cell center at face to get ghost cell center
    for (long i = mesh->n_inner_cells; i < mesh->n_inner_cells + mesh->n_ghost_cells; ++i) {
        double xf[N_DIMS], nf[N_DIMS];
        center(xf, x, &node[i_node[i]], i_node[i + 1] - i_node[i], N_DIMS - 1);
        normal(nf, x, &node[i_node[i]], i_node[i + 1] - i_node[i]);

        const long inner = mesh->cell.cell[mesh->cell.i_cell[i]];
        const double *xi = mesh->cell.center[inner];

        double dx[N_DIMS];
        for (long d = 0; d < N_DIMS; ++d) dx[d] = xf[d] - xi[d];

        const double dot = array_dot(dx, nf, N_DIMS);
        for (long d = 0; d < N_DIMS; ++d) mesh->cell.center[i][d] = xi[d] + 2 * dot * nf[d];
    }

    sync_all(mesh, N_DIMS, mesh->cell.center);
}

static void compute_face_normals(Mesh *mesh) {
    const ALIAS(x, mesh->node.x);
    const ALIAS(node, mesh->face.node);
    const ALIAS(i_node, mesh->face.i_node);
    mesh->face.normal = memory_calloc(mesh->n_faces, sizeof(*mesh->face.normal));
    for (long i = 0; i < mesh->n_faces; ++i) {
        // compute normal
        normal(mesh->face.normal[i], x, &node[i_node[i]], i_node[i + 1] - i_node[i]);

        // compute vector from face center to (inner) cell center
        double FC[N_DIMS];
        for (long d = 0; d < N_DIMS; ++d)
            FC[d] = mesh->cell.center[mesh->face.cell[i][0]][d] - mesh->face.center[i][d];

        // make sure that normal points outwards
        double dot = array_dot(FC, mesh->face.normal[i], N_DIMS);
        assert(fabs(dot) > EPS && "implausible zero dot product");
        if (dot > 0)
            for (long d = 0; d < N_DIMS; ++d) mesh->face.normal[i][d] *= -1;
    }
}

static void compute_cell_projections(Mesh *mesh) {
    const ALIAS(cell, mesh->face.cell);
    const ALIAS(normal, mesh->face.normal);
    const ALIAS(area, mesh->face.area);
    mesh->cell.projection = memory_calloc(mesh->n_inner_cells, sizeof(*mesh->cell.projection));
    for (long i = 0; i < mesh->n_faces; ++i)
        for (long j = 0; j < 2 && cell[i][j] < mesh->n_inner_cells; ++j)
            for (long d = 0; d < N_DIMS; ++d)
                mesh->cell.projection[cell[i][j]][d] += fabs(normal[i][d]) * area[i] / 2;
}

static void compute_face_gauss_weights(Mesh *mesh) {
    const long n_faces = mesh->n_faces;
    const ALIAS(cell, mesh->face.cell);
    const ALIAS(i_node, mesh->cell.i_node);
    const ALIAS(volume, mesh->cell.volume);
    mesh->face.gauss_weight = memory_calloc(n_faces, sizeof(*mesh->face.gauss_weight));
    for (long i = 0; i < n_faces; ++i) {
        for (long j = 0; j < 2 && cell[i][j] < mesh->n_inner_cells; ++j) {
            const long n_nodes = i_node[cell[i][j] + 1] - i_node[cell[i][j]];
            switch (N_DIMS) {
                case 2:
                    switch (n_nodes) {
                        case 3:
                            mesh->face.gauss_weight[i][j] = volume[cell[i][j]] / 3;
                            break;
                        case 4:
                            mesh->face.gauss_weight[i][j] = volume[cell[i][j]] / 6;
                            break;
                        default:
                            assert(!"unsupported element");
                    }
                    break;
                default:
                    assert(!"unsupported dimension");
            }
        }
    }
}

static void compute_cell_gauss_weights(Mesh *mesh) {
    const long n_inner_cells = mesh->n_inner_cells;
    const ALIAS(i_node, mesh->cell.i_node);
    const ALIAS(volume, mesh->cell.volume);
    mesh->cell.gauss_weight = memory_calloc(n_inner_cells, sizeof(*mesh->cell.gauss_weight));
    for (long i = 0; i < n_inner_cells; ++i) {
        const long n_nodes = i_node[i + 1] - i_node[i];
        switch (N_DIMS) {
            case 2:
                switch (n_nodes) {
                    case 3:
                        mesh->cell.gauss_weight[i] = 0;
                        break;
                    case 4:
                        mesh->cell.gauss_weight[i] = volume[i] / 3;
                        break;
                    default:
                        assert(!"unsupported element");
                }
                break;
            default:
                assert(!"unsupported dimension");
        }
    }
}

static void compute_face_gradient_weights(Mesh *mesh) {
    // build periodic cell dict
    const ALIAS(bc, mesh->entity.bc.name);
    fcleanup(dict_free) Dict periodic = dict_create(mesh->n_cells - mesh->n_inner_cells);
    for (long e = 0; e < mesh->n_entities; ++e) {
        if (!bc[e] || strcmp(bc[e], "periodic")) continue;
        for (long ghost = mesh->entity.j_cell[e]; ghost < mesh->entity.j_cell[e + 1]; ++ghost) {
            assert(mesh->cell.i_cell[ghost + 1] - mesh->cell.i_cell[ghost] == 2);
            const long *cell = &mesh->cell.cell[mesh->cell.i_cell[ghost]];  // inner -> outer
            dict_insert(&periodic, cell, 2, &ghost, 1);
        }
    }

    // compute entries of upper triangular matrix (Haselbacher 2000, eq. (13))
    const ALIAS(cell, mesh->face.cell);
    const ALIAS(center, mesh->cell.center);
    cleanup double *r11 = memory_calloc(mesh->n_inner_cells, sizeof(*r11));
    cleanup double *r12 = memory_calloc(mesh->n_inner_cells, sizeof(*r12));
    cleanup double *r22 = memory_calloc(mesh->n_inner_cells, sizeof(*r22));
    for (long i = 0; i < mesh->n_faces; ++i) {
        double dx[N_DIMS] = {};
        long *outer;
        if (dict_lookup(&periodic, cell[i], 2, &outer)) {
            for (long d = 0; d < N_DIMS; ++d) dx[d] = center[*outer][d] - center[cell[i][0]][d];
        } else {
            for (long d = 0; d < N_DIMS; ++d) dx[d] = center[cell[i][1]][d] - center[cell[i][0]][d];
        }

        for (long j = 0; j < 2 && cell[i][j] < mesh->n_inner_cells; ++j) {
            r11[cell[i][j]] += dx[0] * dx[0];
            r12[cell[i][j]] += dx[0] * dx[1];
            r22[cell[i][j]] += dx[1] * dx[1];
            for (long d = 0; d < N_DIMS; ++d) dx[d] *= -1;
        }
    }
    for (long i = 0; i < mesh->n_inner_cells; ++i) {
        r11[i] = sqrt(r11[i]);
        r12[i] = r12[i] / r11[i];
        r22[i] = sqrt(r22[i] - r12[i] * r12[i]);
    }

    // compute gradient weights (Haselbacher 2000, eq. (12))
    mesh->face.gradient_weight = memory_calloc(mesh->n_faces, sizeof(*mesh->face.gradient_weight));
    for (long i = 0; i < mesh->n_faces; ++i) {
        double dx[N_DIMS] = {};
        long *outer;
        if (dict_lookup(&periodic, cell[i], 2, &outer)) {
            for (long d = 0; d < N_DIMS; ++d) dx[d] = center[*outer][d] - center[cell[i][0]][d];
        } else {
            for (long d = 0; d < N_DIMS; ++d) dx[d] = center[cell[i][1]][d] - center[cell[i][0]][d];
        }

        for (long j = 0; j < 2 && cell[i][j] < mesh->n_inner_cells; ++j) {
            double a1 = dx[0] / (r11[cell[i][j]] * r11[cell[i][j]]);
            double a2 = (dx[1] - r12[cell[i][j]] / r11[cell[i][j]] * dx[0]) /
                        (r22[cell[i][j]] * r22[cell[i][j]]);
            if (isnan(a1)) a1 = 0;
            if (isnan(a2)) a2 = 0;
            mesh->face.gradient_weight[i][j][X] = a1 - r12[cell[i][j]] / r11[cell[i][j]] * a2;
            mesh->face.gradient_weight[i][j][Y] = a2;
            for (long d = 0; d < N_DIMS; ++d) dx[d] *= -1;
        }
    }
}

static void compute_face_reconstruction(Mesh *mesh) {
    // build periodic cell dict
    const ALIAS(bc, mesh->entity.bc.name);
    fcleanup(dict_free) Dict periodic = dict_create(mesh->n_cells - mesh->n_inner_cells);
    for (long e = 0; e < mesh->n_entities; ++e) {
        if (!bc[e] || strcmp(bc[e], "periodic")) continue;
        for (long ghost = mesh->entity.j_cell[e]; ghost < mesh->entity.j_cell[e + 1]; ++ghost) {
            assert(mesh->cell.i_cell[ghost + 1] - mesh->cell.i_cell[ghost] == 2);
            const long *cell = &mesh->cell.cell[mesh->cell.i_cell[ghost]];  // inner -> outer
            dict_insert(&periodic, cell, 2, &ghost, 1);
        }
    }

    const long n_faces = mesh->n_faces;
    const ALIAS(xc, mesh->cell.center);
    const ALIAS(xf, mesh->face.center);
    double(*reconstruction)[2][N_DIMS] = memory_calloc(n_faces, sizeof(*reconstruction));
    for (long i = 0; i < n_faces; ++i) {
        long cell[2] = {mesh->face.cell[i][0], mesh->face.cell[i][1]};
        long *outer;
        if (dict_lookup(&periodic, cell, 2, &outer)) cell[1] = *outer;

        for (long j = 0; j < 2; ++j)
            for (long d = 0; d < N_DIMS; ++d) reconstruction[i][j][d] = xf[i][d] - xc[cell[j]][d];
    }
    mesh->face.reconstruction = reconstruction;
}

static void compute_cell_reconstruction(Mesh *mesh) {
    // build cell to face dict
    fcleanup(dict_free) Dict c2f = dict_create(mesh->n_faces);
    for (long i = 0; i < mesh->n_faces; ++i) {
        const long *cell = mesh->face.cell[i];
        const long key[2] = {MIN(cell[0], cell[1]), MAX(cell[0], cell[1])};
        dict_insert(&c2f, key, 2, &i, 1);
    }

    const long n_inner_cells = mesh->n_inner_cells;
    const ALIAS(i_cell, mesh->cell.i_cell);
    const ALIAS(cell, mesh->cell.cell);
    const ALIAS(xc, mesh->cell.center);
    const ALIAS(xf, mesh->face.center);
    double(*reconstruction)[N_DIMS] = memory_calloc(i_cell[n_inner_cells], sizeof(*reconstruction));
    for (long j = 0; j < n_inner_cells; ++j) {
        for (long i = i_cell[j]; i < i_cell[j + 1]; ++i) {
            const long key[2] = {MIN(j, cell[i]), MAX(j, cell[i])};
            long *face, n_faces = dict_lookup(&c2f, key, 2, &face);
            assert(n_faces == 1 && "cell reconstruction error");
            for (long d = 0; d < N_DIMS; ++d) reconstruction[i][d] = xf[*face][d] - xc[j][d];
        }
    }
    mesh->cell.reconstruction = reconstruction;
}

static void compute_entity_to_face(Mesh *mesh) {
    mesh->entity.j_face = memory_calloc(mesh->n_entities + 1, sizeof(*mesh->entity.j_face));
    for (long e = 0; e < mesh->n_entities; ++e) {
        if (mesh->entity.j_cell[e] < mesh->n_inner_cells) continue;
        if (!mesh->entity.bc.name[e] || !strcmp(mesh->entity.bc.name[e], "periodic")) continue;
        mesh->entity.j_face[e + 1] = mesh->entity.j_cell[e + 1] - mesh->entity.j_cell[e];
    }
    for (long e = 0; e < mesh->n_entities; ++e)
        mesh->entity.j_face[e + 1] += mesh->entity.j_face[e];
    for (long e = 0; e < mesh->n_entities + 1; ++e) mesh->entity.j_face[e] += mesh->n_inner_faces;
}

static long common_cell_nodes(const Mesh *mesh, const long cell1, const long cell2, long *node) {
    // WARNING: only works for unique node indices
    long n_common = 0;
    for (long i1 = mesh->cell.i_node[cell1]; i1 < mesh->cell.i_node[cell1 + 1]; ++i1) {
        for (long i2 = mesh->cell.i_node[cell2]; i2 < mesh->cell.i_node[cell2 + 1]; ++i2) {
            if (mesh->cell.node[i1] == mesh->cell.node[i2]) {
                if (node) node[n_common] = mesh->cell.node[i1];
                n_common += 1;
                break;
            }
        }
    }
    return n_common;
}

static double volume(const double (*x)[N_DIMS], const long *node, const long n_nodes,
                     const long dim) {
    double v = 0, AB[3] = {}, AC[3] = {}, AD[3] = {}, c[3] = {};
    switch (dim) {
        case 1:
            switch (n_nodes) {
                case 2:  // line
                    for (long d = 0; d < N_DIMS; ++d) AB[d] = x[node[1]][d] - x[node[0]][d];
                    v += array_norm2(AB, 3);
                    break;
                default:
                    printf("%ld\n", n_nodes);
                    assert(!"unsupported element");
            }
            break;
        case 2:
            switch (n_nodes) {
                case 3:  // triangle
                    for (long d = 0; d < N_DIMS; ++d) AB[d] = x[node[1]][d] - x[node[0]][d];
                    for (long d = 0; d < N_DIMS; ++d) AC[d] = x[node[2]][d] - x[node[0]][d];
                    v += array_norm2(cross(AB, AC, c), 3) / 2;
                    break;
                case 4:  // quadrangle
                    for (long d = 0; d < N_DIMS; ++d) AB[d] = x[node[1]][d] - x[node[0]][d];
                    for (long d = 0; d < N_DIMS; ++d) AC[d] = x[node[2]][d] - x[node[0]][d];
                    for (long d = 0; d < N_DIMS; ++d) AD[d] = x[node[3]][d] - x[node[0]][d];
                    v += array_norm2(cross(AB, AC, c), 3) / 2;
                    v += array_norm2(cross(AC, AD, c), 3) / 2;
                    break;
                default:
                    printf("%ld\n", n_nodes);
                    assert(!"unsupported element");
            }
            break;
        default:
            assert(!"unsupported dimension");
    }
    return v;
}

static void center(double *center, const double (*x)[N_DIMS], const long *node, const long n_nodes,
                   const long dim) {
    double c1[3] = {}, c2[3] = {}, v1, v2, AB[3] = {}, AC[3] = {}, AD[3] = {}, c[3] = {};
    switch (dim) {
        case 1:
            switch (n_nodes) {
                case 2:  // line
                    for (long d = 0; d < N_DIMS; ++d)
                        center[d] = (x[node[0]][d] + x[node[1]][d]) / 2;
                    break;
                default:
                    assert(!"unsupported element");
            }
            break;
        case 2:
            switch (n_nodes) {
                case 3:  // triangle
                    for (long d = 0; d < N_DIMS; ++d)
                        center[d] = (x[node[0]][d] + x[node[1]][d] + x[node[2]][d]) / 3;
                    break;
                case 4:  // quadrangle
                    for (long d = 0; d < N_DIMS; ++d) AB[d] = x[node[1]][d] - x[node[0]][d];
                    for (long d = 0; d < N_DIMS; ++d) AC[d] = x[node[2]][d] - x[node[0]][d];
                    for (long d = 0; d < N_DIMS; ++d) AD[d] = x[node[3]][d] - x[node[0]][d];
                    v1 = array_norm2(cross(AB, AC, c), 3) / 2;
                    v2 = array_norm2(cross(AC, AD, c), 3) / 2;

                    for (long d = 0; d < N_DIMS; ++d)
                        c1[d] = (x[node[0]][d] + x[node[1]][d] + x[node[2]][d]) / 3;
                    for (long d = 0; d < N_DIMS; ++d)
                        c2[d] = (x[node[0]][d] + x[node[2]][d] + x[node[3]][d]) / 3;

                    for (long d = 0; d < N_DIMS; ++d)
                        center[d] = (v1 * c1[d] + v2 * c2[d]) / (v1 + v2);
                    break;
                default:
                    assert(!"unsupported element");
            }
            break;
        default:
            assert(!"unsupported dimension");
    }
}

static void normal(double *normal, const double (*x)[N_DIMS], const long *node,
                   const long n_nodes) {
    double AB[3] = {}, c[3] = {}, norm2;
    switch (N_DIMS) {
        case 2:
            switch (n_nodes) {
                case 2:  // line
                    for (long d = 0; d < N_DIMS; ++d) AB[d] = x[node[1]][d] - x[node[0]][d];
                    norm2 = array_norm2(cross(AB, (double[]){0, 0, 1}, c), 3);
                    for (long d = 0; d < N_DIMS; ++d) normal[d] = c[d] / norm2;
                    break;
                default:
                    assert(!"unsupported element");
            }
            break;
        default:
            assert(!"unsupported dimension");
    }
}

static double *cross(const double *a, const double *b, double *c) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    return c;
}
