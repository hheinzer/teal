#pragma once

#include "vector.h"

enum { MAX_CELL_NODES = 8 };
enum { MAX_CELL_FACES = 6 };
enum { MAX_FACE_NODES = 4 };

// Node types (in order):
// - inner: rank owns this node
// - outer: other rank owns this node; used to define cells
typedef struct {
    int num;
    int num_inner;
    long *global;
    vector *coord;
} MeshNodes;

// Compressed sparse row graph
typedef struct {
    int num;
    int *off;   // [num + 1] offsets into idx
    long *idx;  // [off[num]] indices
} Graph;

// Cell types (in order):
// - inner: rank owns this cell
// - boundary: rank owns this cell; used for boundary conditions
// - periodic: rank owns this cell; used for periodic boundaries
// - neighbor: other rank owns this cell; used for halo exchange
typedef struct {
    int num;
    int num_inner;
    int off_boundary;  // num_inner + num_boundary
    int off_periodic;  // num_inner + num_boundary + num_periodic
    Graph node;        // cell-to-node connectivity
    Graph cell;        // cell-to-cell connectivity
    double *volume, sum_volume;
    vector *center;
    vector *projection;  // axis-aligned projection of cell volume
    vector *offset;      // cell-to-face offset
} MeshCells;

typedef struct {
    long left, right;  // left is always inner; right can be inner or outer
} Adjacent;

typedef struct {
    vector normal, tangent1, tangent2;
} Basis;

typedef struct {
    vector left, right;
} Weight;

typedef struct {
    vector left, right;
} Offset;

typedef struct {
    vector unit;
    double norm;
} Correction;

// Face types (in order):
// - inner: face between two inner cells
// - boundary: face between an inner and a boundary cell
// - periodic: face between an inner and a periodic cell
// - neighbor: face between an inner and a neighbor cell
typedef struct {
    int num;
    int num_inner;
    int off_boundary;  // num_inner + num_boundary
    int off_periodic;  // num_inner + num_boundary + num_periodic
    Graph node;        // face-to-node connectivity
    Adjacent *cell;    // face-to-cell connectivity
    double *area;
    vector *center;
    Basis *basis;            // local orthonormal basis
    Weight *weight;          // least-squares weights for gradient reconstruction
    Offset *offset;          // cell-to-face offset
    Correction *correction;  // correction vector for gradient averaging
} MeshFaces;

typedef char Name[128];

// Entity types (in order):
// - inner: group of inner cells
// - boundary: group of boundary cells
typedef struct {
    int num;
    int num_inner;
    Name *name;
    int *cell_off;  // entity-to-cell offsets
    int *face_off;  // entity-to-face offsets
} MeshEntities;

typedef struct {
    int num;
    int *cell_off;  // periodic-to-cell offsets
    int *face_off;  // periodic-to-face offsets
    vector *translation;
} MeshPeriodics;

typedef struct {
    int num;
    int *rank;
    int *recv_off;  // neighbor-to-cell offsets for receiving
    Graph send;     // neighbor-to-cell graph for sending
} MeshNeighbors;

typedef struct {
    MeshNodes nodes;
    MeshCells cells;
    MeshFaces faces;
    MeshEntities entities;
    MeshPeriodics periodics;
    MeshNeighbors neighbors;
} Mesh;

// Create a Cartesian mesh with optional per-axis periodicity.
Mesh *mesh_create(vector min_coord, vector max_coord, const long *num_cells, const int *periodic,
                  int ndims);

// Read a mesh from a file.
Mesh *mesh_read(const char *fname);

// Split a mesh entity by a plane into below (<entity>-a) and above the plane (<entity>-b).
void mesh_split(Mesh *mesh, const char *entity, vector root, vector normal, int ndims);

// Generate derived mesh data.
void mesh_generate(Mesh *mesh);

// Print a brief mesh summary.
void mesh_summary(const Mesh *mesh);

// Write a mesh to a file.
void mesh_write(const Mesh *mesh, const char *fname);

// Release all mesh allocations.
void mesh_free(Mesh *mesh);
