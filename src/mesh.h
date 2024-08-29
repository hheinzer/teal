#pragma once

#include "teal.h"

#define MAX_CELL_NODES 8
#define MAX_CELL_FACES 6
#define MAX_FACE_NODES 4

enum Dimension { X, Y, Z, N_DIMS };

enum Side { L, R, N_SIDES };

typedef struct Mesh Mesh;

typedef void Modify(Vector3d x);

struct Mesh {
    long n_inner_nodes, n_sync_nodes, n_nodes;
    long n_inner_cells, n_ghost_cells, n_sync_cells, n_cells;
    long n_inner_faces, n_ghost_faces, n_sync_faces, n_faces;
    long n_entities;
    double volume;

    struct {
        long *global;
        Vector3d *coord;
    } node;

    struct {
        long *i_node, *node;
        long *i_cell, *cell;
        double *volume;
        Vector3d *center;
        Vector3d *projection;
        Vector3d *to_cell;
    } cell;

    struct {
        long *i_node, *node;
        Vector2l *cell;
        double *area;
        Vector3d *center;
        Matrix3d *basis;
        Vector3d (*to_cell)[N_SIDES];
        Vector3d *weight;
        Vector4d *correction;
    } face;

    struct {
        String *name;
        long *j_cell, *j_face;
        Vector3d *offset;
    } entity;

    struct {
        long *j_recv, *i_send, *send;
    } sync;
};

Mesh *mesh_create(const Vector3d a, const Vector3d b, const Vector3l n_cells);

Mesh *mesh_read(const char *fname);

void mesh_modify(Mesh *mesh, Modify *modify);

void mesh_split(Mesh *mesh, const char *entity, const Vector3d root, const Vector3d normal);

void mesh_periodic(Mesh *mesh, const char *entity_a, const char *entity_b);

void mesh_generate(Mesh **mesh);

void mesh_print(const Mesh *mesh);

void mesh_write(const Mesh *mesh, const char *prefix);

void mesh_free(Mesh **mesh);
