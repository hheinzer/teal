#include "equations.h"
#include "sync.h"
#include "teal/arena.h"
#include "teal/assert.h"

void *equations_gradient(const Equations *eqns, void *variable_)
{
    assert(eqns && variable_);

    int num = eqns->mesh->cells.num;

    int num_faces = eqns->mesh->faces.num;
    int num_inner = eqns->mesh->faces.num_inner;
    int off_ghost = eqns->mesh->faces.off_ghost;
    Adjacent *cell = eqns->mesh->faces.cell;
    vector *weight = eqns->mesh->faces.weight;

    int stride = eqns->variables.stride;
    scalar(*variable)[stride] = variable_;
    vector(*gradient)[stride] = arena_calloc(num, sizeof(*gradient));

    Request req = sync_variables(eqns, variable, stride);

    for (int i = 0; i < num_inner; i++) {
        int left = cell[i].left;
        int right = cell[i].right;
        for (int j = 0; j < stride; j++) {
            scalar diff = variable[right][j] - variable[left][j];
            gradient[left][j].x += diff * weight[i].x;
            gradient[left][j].y += diff * weight[i].y;
            gradient[left][j].z += diff * weight[i].z;
            gradient[right][j].x += diff * weight[i].x;
            gradient[right][j].y += diff * weight[i].y;
            gradient[right][j].z += diff * weight[i].z;
        }
    }
    for (int i = num_inner; i < off_ghost; i++) {
        int left = cell[i].left;
        int right = cell[i].right;
        for (int j = 0; j < stride; j++) {
            scalar diff = variable[right][j] - variable[left][j];
            gradient[left][j].x += diff * weight[i].x;
            gradient[left][j].y += diff * weight[i].y;
            gradient[left][j].z += diff * weight[i].z;
        }
    }

    sync_wait(eqns, req.recv);

    for (int i = off_ghost; i < num_faces; i++) {
        int left = cell[i].left;
        int right = cell[i].right;
        for (int j = 0; j < stride; j++) {
            scalar diff = variable[right][j] - variable[left][j];
            gradient[left][j].x += diff * weight[i].x;
            gradient[left][j].y += diff * weight[i].y;
            gradient[left][j].z += diff * weight[i].z;
        }
    }

    sync_wait(eqns, req.send);

    return gradient;
}
