#pragma once

enum { MAX_CELL_NODES = 8 };
enum { MAX_CELL_FACES = 6 };
enum { MAX_FACE_NODES = 4 };

typedef struct {
    double *x, *y, *z;
} Vector;

// node order:
// - inner: rank owns this node
// - outer: other owns this node, is used to specify cells
typedef struct {
    long num;
    long num_inner;
    long *global;
    Vector coord;
} MeshNodes;

// compressed sparse row graph
typedef struct {
    long num;
    long *off;  // [num + 1] offsets into idx
    long *idx;  // [off[num]] indices
} Graph;

// cell order:
// - inner: rank owns this cell
// - boundary: rank owns this cell, is used to specify boundary conditions
// - periodic: rank owns this cell, is used to specify a periodic boundary condition; the connected
//   cell can be on this rank or on another rank
// - neighbor: other rank owns this cell, used to receive cell values from neighbor ranks
typedef struct {
    long num;
    long num_inner;
    long off_boundary;  // num_inner + num_boundary
    long off_periodic;  // num_inner + num_boundary + num_periodic
    Graph node;         // cell-to-node connectivity
    Graph cell;         // cell-to-cell connectivity
    double *volume, sum_volume;
    Vector center;
    Vector projection;  // axis-aligned projection of cell volume
    Vector offset;      // cell-to-face offset
} MeshCells;

typedef struct {
    long *left, *right;  // left is always inner
} Adjacent;

typedef struct {
    Vector normal, tangent1, tangent2;
} Basis;

typedef struct {
    Vector left, right;
} Weight;

typedef struct {
    Vector left, right;
} Offset;

typedef struct {
    Vector unit;
    double *norm;
} Correction;

// face order:
// - inner: face between two inner cells
// - boundary: face between an inner and a boundary cell
// - periodic: face between and inner and a periodic cell
// - neighbor: face between and inner and a neighbor cell
typedef struct {
    long num;
    long num_inner;
    long off_boundary;  // num_inner + num_boundary
    long off_periodic;  // num_inner + num_boundary + num_periodic
    Graph node;         // face-to-node connectivity
    Adjacent cell;      // face-to-cell connectivity
    double *area;
    Vector center;
    Basis basis;            // local orthonormal basis
    Weight weight;          // least-squares weights for gradient reconstruction
    Offset offset;          // cell-to-face offset
    Correction correction;  // correction vector for gradient averaging
} MeshFaces;

// entity order:
// - inner: group of inner cells
// - boundary: group of boundary cells
typedef struct {
    int num;
    int num_inner;
    char (*name)[128];
    long *cell;  // entity-to-cell offsets
    long *face;  // entity-to-face offsets
} MeshEntities;

typedef struct {
    Vector x, y, z;
} Matrix;

// WARN: not the final design
typedef struct {
    int num;
    long *cell;  // periodic-to-cell offsets
    long *face;  // periodic-to-face offsets
    Vector translation;
    Matrix rotation;
} MeshPeriodics;

typedef struct {
    int num;
    int *rank;
    long *recv;  // neighbor-to-cell offsets for receiving
    Graph send;  // neighbor-to-cell graph for sending
} MeshNeighbors;

typedef struct {
    MeshNodes nodes;
    MeshCells cells;
    MeshFaces faces;
    MeshEntities entities;
    MeshPeriodics periodics;
    MeshNeighbors neighbors;
} Mesh;

Mesh *mesh_create(const double *min_coord, const double *max_coord, const long *num_cells,
                  const int *periodic, int ndims);

Mesh *mesh_read(const char *fname);

void mesh_split(Mesh *mesh, const char *entity, const double *root, const double *normal,
                int ndims);

void mesh_generate(Mesh *mesh);

void mesh_summary(const Mesh *mesh);

void mesh_write(const Mesh *mesh, const char *prefix);

void mesh_free(Mesh *mesh);
