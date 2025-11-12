#include "equations.h"
#include "sync.h"
#include "teal/arena.h"
#include "teal/assert.h"
#include "teal/vector.h"

void *equations_gradient(const Equations *eqns, void *variable_)
{
    assert(eqns && variable_);

    number num = eqns->mesh->cells.num;

    number num_faces = eqns->mesh->faces.num;
    number num_inner = eqns->mesh->faces.num_inner;
    number off_ghost = eqns->mesh->faces.off_ghost;
    Adjacent *cell = eqns->mesh->faces.cell;
    vector *weight = eqns->mesh->faces.weight;

    number stride = eqns->variables.stride;
    scalar(*variable)[stride] = variable_;
    vector(*gradient)[stride] = arena_calloc(num, sizeof(*gradient));

    Request req = sync_variables(eqns, variable, stride);

    for (number i = 0; i < num_inner; i++) {
        number left = cell[i].left;
        number right = cell[i].right;
        for (number j = 0; j < stride; j++) {
            scalar diff = variable[right][j] - variable[left][j];
            vector_inc(&gradient[left][j], vector_mul(diff, weight[i]));
            vector_inc(&gradient[right][j], vector_mul(diff, weight[i]));
        }
    }
    for (number i = num_inner; i < off_ghost; i++) {
        number left = cell[i].left;
        number right = cell[i].right;
        for (number j = 0; j < stride; j++) {
            scalar diff = variable[right][j] - variable[left][j];
            vector_inc(&gradient[left][j], vector_mul(diff, weight[i]));
        }
    }

    sync_wait(eqns, req.recv);

    for (number i = off_ghost; i < num_faces; i++) {
        number left = cell[i].left;
        number right = cell[i].right;
        for (number j = 0; j < stride; j++) {
            scalar diff = variable[right][j] - variable[left][j];
            vector_inc(&gradient[left][j], vector_mul(diff, weight[i]));
            // vector_inc(&gradient[right][j], vector_mul(diff, weight[i])); // TODO: needed?
        }
    }

    sync_wait(eqns, req.send);

    return gradient;
}
