#pragma once

#define NAMELEN 128

#define MAX_CELL_NODES 8
#define MAX_CELL_FACES 6
#define MAX_FACE_NODES 4

enum Dimension { X, Y, Z, N_DIMS };

enum Side { L, R, N_SIDES };

/* Modify the mesh coordinates 'x'. */
typedef void Modify(double *x);

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

/* Create a pseudo-structured mesh, where 'x0' is the minimum coordinates, 'x1' is the maximum
 * coordinates, and 'n_cells' is the number of cells per spatial dimension. */
Mesh mesh_create(const double *x0, const double *x1, const long *n_cells);

/* Read a mesh from the file 'fname'. */
Mesh mesh_read(const char *fname);

/* Deallocate memory of 'mesh' and set all fields to 0. */
void mesh_free(Mesh *mesh);

/* Modify the nodes of 'mesh' using a 'modify' function. */
void mesh_modify(Mesh *mesh, Modify *modify);

/* Split an 'entity' of 'mesh' using a splitting plane defined by a 'root' coordinate and a
 * 'normal' vector. Cells located on the positive side of the plane (w.r.t the 'normal') are
 * assigned to the first new entity and cells on the negative side to the second. The names of the
 * new entities are the same as 'entity' but with "-0" and "-1" appended. */
void mesh_split(Mesh *mesh, const char *entity, const double *root, const double *normal);

/* Make entities 'entity0' and 'entity1' of 'mesh' periodic. */
void mesh_periodic(Mesh *mesh, const char *entity0, const char *entity1);

/* Generate the 'mesh', partition it, and compute all cell and face properties. After calling this
 * function, the mesh may no longer be modified. */
void mesh_generate(Mesh *mesh);

/* Print a summary of 'mesh'. */
void mesh_print(const Mesh *mesh);

/* Write the 'mesh' to an HDF5 file in the VTKHDF format. The file name is based on 'prefix'. */
void mesh_write(const Mesh *mesh, const char *prefix);
