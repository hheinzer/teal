/*
 * The mesh is composed of nodes, cells, faces, entities, and neighbor communication meta data.
 *
 * Mesh components are ordered as follows:
 * - Nodes: inner nodes first, then neighbor nodes
 * - Cells: inner cells, then ghost cells (for boundary conditions), then periodic cells, then
 *   neighbor cells (used for communication only)
 * - Faces: inner faces, then ghost faces, then neighbor faces
 * - Entities: inner entities (subsets of the interior mesh), then ghost entities, then periodic
 *   entities
 * - Neighbors: periodic neighbors first, then MPI neighbors
 * This ordering guarantees deterministic global indices and enables slicing by category.
 *
 * Connectivity is stored in compressed sparse row (CSR) form:
 * - For a graph with N rows, off has length N+1, idx has length off[N]
 * - Row r spans idx[ off[r] ... off[r+1]-1 ]
 *
 * Entities represent groups of cells/faces. Each entity has a name, a contiguous range of
 * cells/faces, and (for periodic boundaries) a geometric offset. Face ranges are only defined for
 * non-inner entities, while inner entities always refer to subsets of interior mesh (useful for
 * specifying initial conditions).
 *
 * Neighbor meta data describes MPI ranks adjacent to the own rank and contains send/recv graphs
 * that define which cells are exchanged during communication.
 *
 * All geometry information is computed lazily by `mesh_build()` after mesh creation or reading. The
 * mesh may be modified before, but not after, as derived data could be invalidated.
 *
 * Supported cell types are:
 * - 2D: triangle (3 nodes), quadrangle (4 nodes)
 * - 3D: tetrahedron (4 nodes), pyramid (5 nodes), wedge (6 nodes), hexahedron (8 nodes)
 *
 * Periodic boundaries must be representable through translation; rotations is not supported.
 */
#pragma once

#include "teal.h"

enum { MAX_CELL_NODES = 8 };
enum { MAX_CELL_FACES = 6 };
enum { MAX_FACE_NODES = 4 };

typedef struct {
    long num;
    long num_inner;
    long *global;
    vector *coord;
} MeshNodes;

typedef struct {
    long *off, *idx;  // CSR graph
} MeshGraph;

typedef struct {
    long num;
    long num_inner;
    long num_ghost;
    long num_periodic;
    MeshGraph node;  // cell-to-node connectivity
    MeshGraph cell;  // cell-to-cell connectivity
    double *volume;
    vector *center;
    vector *projection;  // axis-aligned projection of cell volume
} MeshCells;

typedef struct {
    long left, right;  // adjacent face cells (left: always inner; right: inner or outer)
} MeshFaceCell;

typedef struct {
    long num;
    long num_inner;
    long num_ghost;
    MeshGraph node;      // face-to-node connectivity
    MeshFaceCell *cell;  // face-to-cell connectivity
    double *area;
    vector *center;
    matrix *basis;   // local orthonormal basis (x: normal; y,z: tangents)
    vector *weight;  // least-squares weights for gradient reconstruction
} MeshFaces;

typedef struct {
    long num;
    long num_inner;
    long num_ghost;
    string *name;
    long *cell_off;  // entity-to-cell offsets
    long *face_off;  // entity-to-face offsets
    vector *offset;  // entity translation (zero for non-periodic)
} MeshEntities;

typedef struct {
    long num;
    long *rank;      // neighbor ranks
    long *recv_off;  // neighbor-to-cell offsets for receiving
    MeshGraph send;  // neighbor-to-cell graph for sending
} MeshNeighbors;

typedef struct {
    MeshNodes nodes;
    MeshCells cells;
    MeshFaces faces;
    MeshEntities entities;
    MeshNeighbors neighbors;
} Mesh;

/* Create a distributed Cartesian mesh with optional per-axis periodicity. */
Mesh *mesh_create(vector min_coord, vector max_coord, tuple num_cells, flags periodic);

/* Read a mesh from disk. */
Mesh *mesh_read(const char *fname);

/* Split an entity's cells by a plane through `root` with unit-normal `normal`. */
void mesh_split(Mesh *mesh, const char *entity, vector root, vector normal);

/* Build connectivity, faces, neighbor graphs, geometry, and reconstruction weights. */
void mesh_build(Mesh *mesh);

/* Print a per-rank, human-readable mesh dump. */
void mesh_print(const Mesh *mesh);

void mesh_test(const Mesh *mesh);

/* Print a global mesh summary. */
void mesh_summary(const Mesh *mesh);

/* Create HDF5 file `fname` and write nodes, cells, and entities groups. */
void mesh_write(const Mesh *mesh, const char *fname);
