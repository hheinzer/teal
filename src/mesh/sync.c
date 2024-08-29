#include "sync.h"

#include "teal/memory.h"
#include "teal/sync.h"
#include "teal/utils.h"

#define size(buf) \
    _Generic(buf, long *: sizeof(long), Vector3d *: sizeof(double), String *: sizeof(char))
#define alloc(buf, cnt) buf = memory_calloc(cnt, size(buf))

#define type(buf) _Generic(buf, long *: MPI_LONG, Vector3d *: MPI_DOUBLE, String *: MPI_CHAR)
#define send(buf, cnt) MPI_Send(buf, cnt, type(buf), rank, tag++, sync.comm)
#define recv(buf, cnt) MPI_Recv(buf, cnt, type(buf), 0, tag++, sync.comm, MPI_STATUS_IGNORE)

#define alloc_recv(buf, cnt) (alloc(buf, cnt), recv(buf, cnt))

static void send_mesh(Mesh *mesh, int rank);

static void recv_mesh(Mesh *mesh);

void sync_all(const Mesh *mesh, double *u, long ldu)
{
    const alias(j_recv, mesh->sync.j_recv);
    const alias(i_send, mesh->sync.i_send);
    const alias(send, mesh->sync.send);

    smart MPI_Request *req_recv = memory_calloc(sync.size, sizeof(*req_recv));  // NOLINT
    smart MPI_Request *req_send = memory_calloc(sync.size, sizeof(*req_send));  // NOLINT
    smart double *buf = memory_calloc(i_send[sync.size] * ldu, sizeof(*buf));

    sync_irecv(j_recv, req_recv, u, ldu);
    sync_isend(i_send, send, u, req_send, buf, ldu);
    sync_waitall(req_recv);
    sync_waitall(req_send);
}

void sync_mesh(Mesh *mesh, int rank)
{
    if (rank > 0)
        send_mesh(mesh, rank);
    else
        recv_mesh(mesh);
}

static void send_mesh(Mesh *mesh, int rank)
{
    int tag = 0;
    send(&mesh->n_nodes, 1);
    send(&mesh->n_inner_nodes, 1);
    send(&mesh->n_sync_nodes, 1);
    send(&mesh->n_cells, 1);
    send(&mesh->n_inner_cells, 1);
    send(&mesh->n_ghost_cells, 1);
    send(&mesh->n_sync_cells, 1);
    send(&mesh->n_entities, 1);

    send(mesh->node.global, mesh->n_nodes);
    send(mesh->node.coord, mesh->n_nodes * N_DIMS);

    send(mesh->cell.i_node, mesh->n_cells + 1);
    send(mesh->cell.node, mesh->cell.i_node[mesh->n_cells]);
    send(mesh->cell.i_cell, mesh->n_cells + 1);
    send(mesh->cell.cell, mesh->cell.i_cell[mesh->n_cells]);

    send(mesh->entity.name, mesh->n_entities * sizeof(String));
    send(mesh->entity.j_cell, mesh->n_entities + 1);
    send(mesh->entity.offset, mesh->n_entities * N_DIMS);

    send(mesh->sync.j_recv, sync.size + 1);
    send(mesh->sync.i_send, sync.size + 1);
    send(mesh->sync.send, mesh->sync.i_send[sync.size]);
}

static void recv_mesh(Mesh *mesh)
{
    int tag = 0;
    recv(&mesh->n_nodes, 1);
    recv(&mesh->n_inner_nodes, 1);
    recv(&mesh->n_sync_nodes, 1);
    recv(&mesh->n_cells, 1);
    recv(&mesh->n_inner_cells, 1);
    recv(&mesh->n_ghost_cells, 1);
    recv(&mesh->n_sync_cells, 1);
    recv(&mesh->n_entities, 1);

    alloc_recv(mesh->node.global, mesh->n_nodes);
    alloc_recv(mesh->node.coord, mesh->n_nodes * N_DIMS);

    alloc_recv(mesh->cell.i_node, mesh->n_cells + 1);
    alloc_recv(mesh->cell.node, mesh->cell.i_node[mesh->n_cells]);
    alloc_recv(mesh->cell.i_cell, mesh->n_cells + 1);
    alloc_recv(mesh->cell.cell, mesh->cell.i_cell[mesh->n_cells]);

    alloc_recv(mesh->entity.name, mesh->n_entities * sizeof(String));
    alloc_recv(mesh->entity.j_cell, mesh->n_entities + 1);
    alloc_recv(mesh->entity.offset, mesh->n_entities * N_DIMS);

    alloc_recv(mesh->sync.j_recv, sync.size + 1);
    alloc_recv(mesh->sync.i_send, sync.size + 1);
    alloc_recv(mesh->sync.send, mesh->sync.i_send[sync.size]);
}
