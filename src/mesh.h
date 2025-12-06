// The mesh is composed of nodes, cells, faces, entities, and neighbor communication meta data.
//
// Mesh components are ordered as follows:
// - Nodes: inner nodes first, then neighbor nodes
// - Cells: inner cells, then ghost cells (for boundary conditions), then periodic cells, then
//   neighbor cells (used for communication only)
// - Faces: inner faces, then ghost faces, then neighbor faces
// - Entities: inner entities (subsets of the interior mesh), then ghost entities, then periodic
//   entities
// - Neighbors: periodic neighbors first, then MPI neighbors
// This ordering guarantees deterministic ordering and enables slicing by category.
//
// Connectivity is stored in compressed sparse row (CSR) form:
// - For a graph with N rows, off has length N+1, idx has length off[N]
// - Row r spans idx[ off[r] ... off[r+1]-1 ]
//
// Entities represent groups of cells/faces. Each entity has a name, a contiguous range of
// cells/faces, and (for periodic boundaries) a translation. Face ranges are only defined for
// non-inner entities, while inner entities always refer to subsets of interior mesh (useful for
// specifying initial conditions).
//
// Neighbor meta data describes MPI ranks adjacent to the own rank and contains send/recv graphs
// that define which cells are exchanged during communication.
//
// All geometry information is computed by `mesh_generate()` after mesh creation or reading. The
// mesh may be modified before, but not after, as derived data could be invalidated.
//
// Supported cell types are:
// - 2D: triangle (3 nodes), quadrangle (4 nodes)
// - 3D: tetrahedron (4 nodes), pyramid (5 nodes), wedge (6 nodes), hexahedron (8 nodes)
//
// Periodic boundaries must be representable through translation; rotations are not supported.
#pragma once

#include "teal.h"

enum { MAX_CELL_NODES = 8 };
enum { MAX_CELL_FACES = 6 };
enum { MAX_FACE_NODES = 4 };

typedef struct {
    long num;
    long num_inner;
    vector *coord;
} MeshNodes;

typedef struct {
    long *off, *idx;  // CSR graph
} Graph;

typedef struct {
    long num;
    long num_inner;
    long off_ghost;     // num_inner + num_ghost
    long off_periodic;  // num_inner + num_ghost + num_periodic
    Graph node;         // cell-to-node connectivity
    Graph cell;         // cell-to-cell connectivity
    scalar *volume;
    vector *center;
    vector *projection;  // axis-aligned projection of cell volume
    vector *offset;      // cell-to-face offset
    scalar sum_volume;
} MeshCells;

typedef struct {
    long left, right;  // adjacent face cells (left: always inner; right: inner or outer)
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
    scalar norm;
} Correction;

typedef struct {
    long num;
    long num_inner;
    long off_ghost;  // num_inner + num_ghost
    Graph node;      // face-to-node connectivity
    Adjacent *cell;  // face-to-cell connectivity
    scalar *area;
    vector *center;
    Basis *basis;            // local orthonormal basis
    Weight *weight;          // least-squares weights for gradient reconstruction
    Offset *offset;          // cell-to-face offset
    Correction *correction;  // correction vector for gradient averaging
} MeshFaces;

typedef char Name[128];

typedef struct {
    long num;
    long num_inner;
    long off_ghost;  // num_inner + num_ghost
    Name *name;
    long *cell_off;       // entity-to-cell offsets
    long *face_off;       // entity-to-face offsets
    vector *translation;  // entity translation (zero for non-periodic)
} MeshEntities;

typedef struct {
    long num;
    long *rank;
    long *recv_off;  // neighbor-to-cell offsets for receiving
    Graph send;      // neighbor-to-cell graph for sending
} MeshNeighbors;

typedef struct {
    MeshNodes nodes;
    MeshCells cells;
    MeshFaces faces;
    MeshEntities entities;
    MeshNeighbors neighbors;
} Mesh;

// Create a distributed Cartesian mesh with optional per-axis periodicity.
Mesh *mesh_create(vector min_coord, vector max_coord, tuple num_cells, const bool *periodic);

// Read a mesh from disk.
Mesh *mesh_read(const char *fname);

// Split an entity's cells by a plane through `root` with normal `normal`.
void mesh_split(Mesh *mesh, const char *entity, vector root, vector normal);

// Build connectivity, faces, neighbor graphs, geometry, and reconstruction weights.
void mesh_generate(Mesh *mesh);

// Check the mesh integrity and report violations.
void mesh_check(const Mesh *mesh);

// Print a global mesh summary.
void mesh_summary(const Mesh *mesh);

// Write mesh to a VTKHDF file.
void mesh_write(const Mesh *mesh, const char *prefix);
