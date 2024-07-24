#ifndef MESH_H
#define MESH_H

#define NAMELEN 128

#define MAX_CELL_NODES 8
#define MAX_CELL_FACES 6
#define MAX_FACE_NODES 4

enum Dimension { X, Y, Z, N_DIMS };

enum Side { L, R, N_SIDES };

typedef struct Mesh {
    long n_nodes, n_inner_nodes, n_sync_nodes;
    long n_cells, n_inner_cells, n_ghost_cells, n_sync_cells;
    long n_faces, n_inner_faces, n_bound_faces, n_sync_faces;
    long n_entities;
    double volume;

    struct {
        long *idx;
        double (*coord)[N_DIMS];
    } node;

    struct {
        long *i_node, *node;
        long *i_cell, *cell;
        double *volume;
        double (*center)[N_DIMS];
        double (*projection)[N_DIMS];
        double (*reconstruction)[N_DIMS];
    } cell;

    struct {
        long *i_node, *node;
        long (*cell)[N_SIDES];
        double *area;
        double (*center)[N_DIMS];
        double (*normal)[N_DIMS * N_DIMS];
        double (*gradient_weight)[N_DIMS];
        double (*reconstruction)[N_SIDES][N_DIMS];
        double (*gradient_correction)[N_DIMS + 1];
    } face;

    struct {
        char (*name)[NAMELEN];
        long *j_cell, *j_face;
        double (*offset)[N_DIMS];
    } entity;

    struct {
        long *j_recv, *i_send, *send;
    } sync;
} Mesh;

typedef void Modify(double *x);

Mesh mesh_create(const double *x0, const double *x1, const long *n_cells);

Mesh mesh_read(const char *fname);

void mesh_free(Mesh *mesh);

void mesh_modify(Mesh *mesh, Modify *modify);

void mesh_split(Mesh *mesh, const char *entity, const double *root, const double *direction);

void mesh_periodic(Mesh *mesh, const char *entity0, const char *entity1);

void mesh_generate(Mesh *mesh);

void mesh_print(const Mesh *mesh);

void mesh_write(const Mesh *mesh, const char *prefix);

#endif
