#include <assert.h>
#include <limits.h>
#include <string.h>

#include "mesh2.h"
#include "reorder.h"
#include "teal2.h"
#include "utils2.h"

static int entity_index(const Mesh2 *mesh, const char *entity)
{
    for (int i = 0; i < mesh->entities.off_boundary; i++) {
        if (!strcmp(mesh->entities.name[i], entity)) {
            return i;
        }
    }
    teal2_error("invalid entity (%s)", entity);
}

static int reorder_cells(Mesh2 *mesh, Vector root, Vector normal, int index)
{
    int beg = mesh->entities.cell_off[index];
    int end = mesh->entities.cell_off[index + 1];
    int num_cells = end - beg;

    int *map = teal2_calloc(num_cells, sizeof(*map));
    int num = 0;
    for (int i = beg; i < end; i++) {
        Vector mean = {0};
        for (int j = mesh->cells.node.off[i]; j < mesh->cells.node.off[i + 1]; j++) {
            vector2_iadd(&mean, mesh->nodes.coord[mesh->cells.node.idx[j]]);
        }
        vector2_idiv(&mean, mesh->cells.node.off[i + 1] - mesh->cells.node.off[i]);
        map[num++] = (vector2_dot(vector2_sub(mean, root), normal) <= 0);
    }
    assert(num == num_cells);

    int offset = 0;
    for (int i = 0; i < num_cells; i++) {
        offset += map[i];
    }

    int lhs = 0;
    int rhs = 0;
    for (int i = 0; i < num_cells; i++) {
        map[i] = map[i] ? lhs++ : (offset + rhs++);
    }
    assert(lhs == offset && rhs == num_cells - offset);

    mesh2_reorder_cells(mesh, map, beg, end);

    teal2_free(map);
    return offset;
}

static void split_entity(Mesh2 *mesh, int index, int offset)
{
    int num_entities = mesh->entities.num + 1;
    int num_move = mesh->entities.num - index;

    String *name = teal2_realloc(mesh->entities.name, num_entities, sizeof(*name));
    move(&name[index + 1], &name[index], num_move, sizeof(*name));
    strcat(name[index], "-a");
    strcat(name[index + 1], "-b");

    int *cell_off = teal2_realloc(mesh->entities.cell_off, num_entities + 1, sizeof(*cell_off));
    move(&cell_off[index + 2], &cell_off[index + 1], num_move, sizeof(*cell_off));
    cell_off[index + 1] = cell_off[index] + offset;

    Matrix *rotation = teal2_realloc(mesh->entities.rotation, num_entities, sizeof(*rotation));
    move(&rotation[index + 1], &rotation[index], num_move, sizeof(*rotation));

    Vector *translation =
        teal2_realloc(mesh->entities.translation, num_entities, sizeof(*translation));
    move(&translation[index + 1], &translation[index], num_move, sizeof(*translation));

    mesh->entities.num = num_entities;
    if (index < mesh->entities.num_inner) {
        mesh->entities.num_inner += 1;
    }
    mesh->entities.off_boundary += 1;
    mesh->entities.name = name;
    mesh->entities.cell_off = cell_off;
    mesh->entities.rotation = rotation;
    mesh->entities.translation = translation;
}

void mesh2_split(Mesh2 *mesh, const char *entity, Vector root, Vector normal)
{
    assert(mesh && entity);
    int index = entity_index(mesh, entity);
    int offset = reorder_cells(mesh, root, normal, index);
    split_entity(mesh, index, offset);
}
