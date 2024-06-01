#ifndef MESH_H
#define MESH_H

typedef enum Dimension { X, Y, N_DIMS } Dimension;
#define MAX_CELL_NODES 4
#define MAX_CELL_FACES 4
#define MAX_FACE_NODES 2

#define ALIAS(a, m) typeof(m) a = m

typedef void Modify(double *x);
typedef void Function(const double *x, const double time, const double *u, double *f);

typedef struct BCContext {
    const double *state;
    Function *custom;
} BCContext;
typedef void ApplyBC(const BCContext context, const double *n, const double *ui, double *ug);

typedef struct Mesh {
    long n_nodes, n_inner_nodes;
    long n_cells, n_inner_cells, n_ghost_cells;
    long n_faces, n_inner_faces, n_bound_faces;
    long n_entities;
    int rank, size;
    double volume;

    struct {
        double (*x)[N_DIMS];
        long *map;
    } node;

    struct {
        long *i_node, *node;
        long *i_cell, *cell;
        double *volume;
        double (*center)[N_DIMS];
        double (*projection)[N_DIMS];
        double *gauss_weight;
        double (*reconstruction)[N_DIMS];
    } cell;

    struct {
        long *i_node, *node;
        long (*cell)[2];  // inner, outer
        double *area;
        double (*center)[N_DIMS];
        double (*normal)[N_DIMS];              // points from inner to outer cell
        double (*gauss_weight)[2];             // inner, outer
        double (*gradient_weight)[2][N_DIMS];  // inner, outer
        double (*reconstruction)[2][N_DIMS];   // inner, outer
    } face;

    struct {
        char **name;
        long *j_cell, *j_face;
        struct {
            char **name;
            BCContext *context;
            ApplyBC **apply;
        } bc;
    } entity;

    struct {
        long *i_recv, *i_send, *send;
    } sync;
} Mesh;

Mesh mesh_create(const double *x0, const double *x1, const long *n_cells);

Mesh mesh_read(const char *fname);

void mesh_finalize(Mesh *mesh);

void mesh_free(Mesh *mesh);

void mesh_modify(Mesh *mesh, Modify *modify);

void mesh_split_entity(Mesh *mesh, const char *name, const long dim, const double x);

void mesh_set_periodic_condition(Mesh *mesh, const char *name0, const char *name1);

void mesh_set_boundary_condition(Mesh *mesh, const char *name, const char *bc, const double *state,
                                 Function *custom);

void mesh_generate(Mesh *mesh);

void mesh_print(const Mesh *mesh);

void mesh_write(const Mesh *mesh, const char *prefix);

#endif
