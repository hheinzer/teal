#include <assert.h>
#include <limits.h>

#include "mesh2.h"
#include "teal2.h"

static int is_valid(const int *map, int num)
{
    int *count = teal2_calloc(num, sizeof(*count));
    for (int i = 0; i < num; i++) {
        if (0 <= map[i] && map[i] < num) {
            count[map[i]] += 1;
        }
    }

    int min = 1;
    int max = 1;
    int sum = 0;
    for (int i = 0; i < num; i++) {
        min = (count[i] < min) ? count[i] : min;
        max = (count[i] > max) ? count[i] : max;
        sum += count[i];
    }

    teal2_free(count);
    return min == 1 && max == 1 && sum == num;
}

typedef struct {
    long global;
    Vector coord;
} Node;

void mesh2_reorder_nodes(Mesh *mesh, const int *map, int beg, int end)
{
    assert(mesh && is_valid(map, end - beg) && 0 <= beg && beg <= end && end <= mesh->nodes.num);

    int num_nodes = end - beg;
    Node *node = teal2_calloc(num_nodes, sizeof(*node));

    int num = 0;
    for (int i = beg; i < end; i++) {
        node[map[num]].global = mesh->nodes.global[i];
        node[map[num]].coord = mesh->nodes.coord[i];
        num += 1;
    }
    assert(num == num_nodes);

    num = 0;
    for (int i = beg; i < end; i++) {
        mesh->nodes.global[i] = node[num].global;
        mesh->nodes.coord[i] = node[num].coord;
        num += 1;
    }
    assert(num == num_nodes);

    if (mesh->cells.node.off && mesh->cells.node.idx) {
        for (int i = 0; i < mesh->cells.node.off[mesh->cells.num]; i++) {
            int idx = mesh->cells.node.idx[i];
            if (beg <= idx && idx < end) {
                mesh->cells.node.idx[i] = beg + map[idx - beg];
            }
        }
    }

    teal2_free(node);
}

typedef struct {
    int num_nodes;
    int node[MAX_CELL_NODES];
    int num_cells;
    int cell[MAX_CELL_FACES];
} Cell;

void mesh2_reorder_cells(Mesh *mesh, const int *map, int beg, int end)
{
    assert(mesh && is_valid(map, end - beg) && 0 <= beg && beg <= end && end <= mesh->cells.num);

    int num_cells = end - beg;
    Cell *cell = teal2_calloc(num_cells, sizeof(*cell));

    int num = 0;
    for (int i = beg; i < end; i++) {
        cell[map[num]].num_nodes = 0;
        for (int j = mesh->cells.node.off[i]; j < mesh->cells.node.off[i + 1]; j++) {
            assert(cell[map[num]].num_nodes < MAX_CELL_NODES);
            cell[map[num]].node[cell[map[num]].num_nodes++] = mesh->cells.node.idx[j];
        }
        if (mesh->cells.cell.off && mesh->cells.cell.idx) {
            for (int j = mesh->cells.cell.off[i]; j < mesh->cells.cell.off[i + 1]; j++) {
                assert(cell[map[num]].num_cells < MAX_CELL_FACES);
                cell[map[num]].cell[cell[map[num]].num_cells++] = mesh->cells.cell.idx[j];
            }
        }
        num += 1;
    }
    assert(num == num_cells);

    num = 0;
    for (int i = beg; i < end; i++) {
        mesh->cells.node.off[i + 1] = mesh->cells.node.off[i];
        for (int j = 0; j < cell[num].num_nodes; j++) {
            mesh->cells.node.idx[mesh->cells.node.off[i + 1]++] = cell[num].node[j];
        }
        if (mesh->cells.cell.off && mesh->cells.cell.idx) {
            mesh->cells.cell.off[i + 1] = mesh->cells.cell.off[i];
            for (int j = 0; j < cell[num].num_cells; j++) {
                mesh->cells.cell.idx[mesh->cells.cell.off[i + 1]++] = cell[num].cell[j];
            }
        }
        num += 1;
    }
    assert(num == num_cells);

    if (mesh->cells.cell.off && mesh->cells.cell.idx) {
        for (int i = 0; i < mesh->cells.cell.off[mesh->cells.num]; i++) {
            int idx = mesh->cells.cell.idx[i];
            if (beg <= idx && idx < end) {
                mesh->cells.cell.idx[i] = beg + map[idx - beg];
            }
        }
    }

    if (mesh->neighbors.send.off && mesh->neighbors.send.idx) {
        for (int i = 0; i < mesh->neighbors.send.off[mesh->neighbors.num]; i++) {
            int idx = mesh->neighbors.send.idx[i];
            if (beg <= idx && idx < end) {
                mesh->neighbors.send.idx[i] = beg + map[idx - beg];
            }
        }
    }

    teal2_free(cell);
}
