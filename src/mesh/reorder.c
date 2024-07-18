#include "reorder.h"

#include <stdlib.h>

#include "core/memory.h"

void reorder_cells(Mesh *mesh, const long *old2new, const long *new2old)
{
    long *i_node = memory_calloc(mesh->n_cells + 1, sizeof(*i_node));
    long *i_cell = memory_calloc(mesh->n_cells + 1, sizeof(*i_cell));
    long *node = memory_calloc(mesh->cell.i_node[mesh->n_cells], sizeof(*node));
    long *cell = memory_calloc(mesh->cell.i_cell[mesh->n_cells], sizeof(*cell));
    for (long jnew = 0; jnew < mesh->n_cells; ++jnew) {
        const long jold = new2old[jnew];

        i_node[jnew + 1] = i_node[jnew];
        for (long i = mesh->cell.i_node[jold]; i < mesh->cell.i_node[jold + 1]; ++i)
            node[i_node[jnew + 1]++] = mesh->cell.node[i];

        i_cell[jnew + 1] = i_cell[jnew];
        for (long i = mesh->cell.i_cell[jold]; i < mesh->cell.i_cell[jold + 1]; ++i)
            cell[i_cell[jnew + 1]++] = old2new[mesh->cell.cell[i]];
    }
    free(mesh->cell.i_node);
    free(mesh->cell.i_cell);
    free(mesh->cell.node);
    free(mesh->cell.cell);
    mesh->cell.i_node = i_node;
    mesh->cell.i_cell = i_cell;
    mesh->cell.node = node;
    mesh->cell.cell = cell;
}
