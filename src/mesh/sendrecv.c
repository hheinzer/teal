#include "sendrecv.h"

#include <mpi.h>
#include <unistd.h>

#include "core/memory.h"
#include "teal.h"

#define SIZE(buf)                          \
    _Generic(buf,                          \
        long *: sizeof(long),              \
        double(*)[N_DIMS]: sizeof(double), \
        char(*)[NAMELEN]: sizeof(char))
#define ALLOC(buf, cnt) buf = memory_calloc(cnt, SIZE(buf))

#define TYPE(buf) \
    _Generic(buf, long *: MPI_LONG, double(*)[N_DIMS]: MPI_DOUBLE, char(*)[NAMELEN]: MPI_CHAR)
#define SEND(buf, cnt) MPI_Send(buf, cnt, TYPE(buf), rank, tag++, teal.comm)
#define RECV(buf, cnt) MPI_Recv(buf, cnt, TYPE(buf), root, tag++, teal.comm, MPI_STATUS_IGNORE)
#define ALLOC_RECV(buf, cnt) \
    ALLOC(buf, cnt);         \
    RECV(buf, cnt)

static int ready = 0;

static void wake_ranks(void);
static void suspend_rank(void);

void send_mesh(Mesh *mesh, long rank)
{
    if (!ready) wake_ranks();

    int tag = 0;
    SEND(&mesh->n_nodes, 1);
    SEND(&mesh->n_inner_nodes, 1);
    SEND(&mesh->n_sync_nodes, 1);
    SEND(&mesh->n_cells, 1);
    SEND(&mesh->n_inner_cells, 1);
    SEND(&mesh->n_ghost_cells, 1);
    SEND(&mesh->n_sync_cells, 1);
    SEND(&mesh->n_entities, 1);

    SEND(mesh->node.idx, mesh->n_nodes);
    SEND(mesh->node.coord, mesh->n_nodes * N_DIMS);

    SEND(mesh->cell.i_node, mesh->n_cells + 1);
    SEND(mesh->cell.node, mesh->cell.i_node[mesh->n_cells]);
    SEND(mesh->cell.i_cell, mesh->n_cells + 1);
    SEND(mesh->cell.cell, mesh->cell.i_cell[mesh->n_cells]);

    SEND(mesh->entity.name, mesh->n_entities * NAMELEN);
    SEND(mesh->entity.j_cell, mesh->n_entities + 1);
    SEND(mesh->entity.offset, mesh->n_entities * N_DIMS);

    SEND(mesh->sync.j_recv, teal.size + 1);
    SEND(mesh->sync.i_send, teal.size + 1);
    SEND(mesh->sync.send, mesh->sync.i_send[teal.size]);
}

void recv_mesh(Mesh *mesh, long root)
{
    suspend_rank();

    int tag = 0;
    RECV(&mesh->n_nodes, 1);
    RECV(&mesh->n_inner_nodes, 1);
    RECV(&mesh->n_sync_nodes, 1);
    RECV(&mesh->n_cells, 1);
    RECV(&mesh->n_inner_cells, 1);
    RECV(&mesh->n_ghost_cells, 1);
    RECV(&mesh->n_sync_cells, 1);
    RECV(&mesh->n_entities, 1);

    ALLOC_RECV(mesh->node.idx, mesh->n_nodes);
    ALLOC_RECV(mesh->node.coord, mesh->n_nodes * N_DIMS);

    ALLOC_RECV(mesh->cell.i_node, mesh->n_cells + 1);
    ALLOC_RECV(mesh->cell.node, mesh->cell.i_node[mesh->n_cells]);
    ALLOC_RECV(mesh->cell.i_cell, mesh->n_cells + 1);
    ALLOC_RECV(mesh->cell.cell, mesh->cell.i_cell[mesh->n_cells]);

    ALLOC_RECV(mesh->entity.name, mesh->n_entities * NAMELEN);
    ALLOC_RECV(mesh->entity.j_cell, mesh->n_entities + 1);
    ALLOC_RECV(mesh->entity.offset, mesh->n_entities * N_DIMS);

    ALLOC_RECV(mesh->sync.j_recv, teal.size + 1);
    ALLOC_RECV(mesh->sync.i_send, teal.size + 1);
    ALLOC_RECV(mesh->sync.send, mesh->sync.i_send[teal.size]);
}

static void wake_ranks(void)
{
    ready = 1;
    MPI_Request req;
    MPI_Ibcast(&ready, 1, MPI_INT, 0, teal.comm, &req);
    MPI_Wait(&req, MPI_STATUS_IGNORE);
}

static void suspend_rank(void)
{
    MPI_Request req;
    MPI_Ibcast(&ready, 1, MPI_INT, 0, teal.comm, &req);
    while (1) {
        int flag;
        MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
        if (flag) break;
        sleep(1);
    }
}
