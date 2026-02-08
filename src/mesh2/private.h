#pragma once

#include <stdint.h>

#include "matrix2.h"
#include "mesh2.h"
#include "vector2.h"

enum MeshLimits {
    MAX_CELL_NODES = 8,
    MAX_CELL_FACES = 6,
    MAX_FACE_NODES = 4,
};

typedef struct {
    int num;
    int num_inner;
    int64_t *global;
    Vector *coord;
} MeshNodes;

typedef struct {
    int *off, *idx;  // compressed sparse row
} Graph;

typedef struct {
    int num;
    int num_inner;
    int off_boundary;
    int off_periodic;
    Graph node;
    Graph cell;
    double *volume;
    Vector *center;
    Vector *projection;
    Vector *offset;
} MeshCells;

typedef struct {
    int left, right;
} IntPair;

typedef struct {
    Vector left, right;
} VectorPair;

typedef struct {
    Vector unit;
    double norm;
} UnitNorm;

typedef struct {
    int num;
    int num_inner;
    int off_boundary;
    Graph node;
    IntPair *cell_idx;
    double *area;
    Vector *center;
    Matrix *basis;
    VectorPair *weight;
    VectorPair *offset;
    UnitNorm *correction;
} MeshFaces;

typedef char Name[128];

typedef struct {
    int num;
    int num_inner;
    int off_boundary;
    Name *name;
    int *cell_off;
    int *face_off;
    Matrix *rotation;
    Vector *translation;
} MeshEntities;

typedef struct {
    int num;
    int (*tag)[2];
    int *rank;
    int *recv_off;
    Graph send;
} MeshNeighbors;

struct mesh {
    int generated;
    MeshNodes nodes;
    MeshCells cells;
    MeshFaces faces;
    MeshEntities entities;
    MeshNeighbors neighbors;
};

void mesh2_reorder_nodes(Mesh2 *mesh, const int *map, int beg, int end);

void mesh2_reorder_cells(Mesh2 *mesh, const int *map, int beg, int end);
