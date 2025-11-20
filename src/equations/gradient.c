#include <assert.h>

#include "equations.h"
#include "sync.h"
#include "teal/arena.h"

void *equations_gradient(const Equations *eqns, void *variable_)
{
    assert(eqns && variable_);

    long num = eqns->mesh->cells.num;

    long num_faces = eqns->mesh->faces.num;
    long num_inner = eqns->mesh->faces.num_inner;
    long off_ghost = eqns->mesh->faces.off_ghost;
    Adjacent *cell = eqns->mesh->faces.cell;
    vector *weight = eqns->mesh->faces.weight;

    long stride = eqns->variables.stride;
    scalar(*variable)[stride] = variable_;
    vector(*gradient)[stride] = arena_calloc(num, sizeof(*gradient));

    Request req = sync_variables(eqns, variable, stride);

    for (long i = 0; i < num_inner; i++) {
        long left = cell[i].left;
        long right = cell[i].right;
        for (long j = 0; j < stride; j++) {
            scalar diff = variable[right][j] - variable[left][j];
            gradient[left][j].x += diff * weight[i].x;
            gradient[left][j].y += diff * weight[i].y;
            gradient[left][j].z += diff * weight[i].z;
            gradient[right][j].x += diff * weight[i].x;
            gradient[right][j].y += diff * weight[i].y;
            gradient[right][j].z += diff * weight[i].z;
        }
    }
    for (long i = num_inner; i < off_ghost; i++) {
        long left = cell[i].left;
        long right = cell[i].right;
        for (long j = 0; j < stride; j++) {
            scalar diff = variable[right][j] - variable[left][j];
            gradient[left][j].x += diff * weight[i].x;
            gradient[left][j].y += diff * weight[i].y;
            gradient[left][j].z += diff * weight[i].z;
        }
    }

    sync_wait(eqns, req.recv);

    for (long i = off_ghost; i < num_faces; i++) {
        long left = cell[i].left;
        long right = cell[i].right;
        for (long j = 0; j < stride; j++) {
            scalar diff = variable[right][j] - variable[left][j];
            gradient[left][j].x += diff * weight[i].x;
            gradient[left][j].y += diff * weight[i].y;
            gradient[left][j].z += diff * weight[i].z;
        }
    }

    sync_wait(eqns, req.send);

    return gradient;
}
