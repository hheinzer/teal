#include <assert.h>
#include <string.h>

#include "equations.h"
#include "exchange.h"

void equations_gradient(const Equations *eqns, void *primitive_, void *gradient_)
{
    assert(eqns && primitive_ && gradient_);

    int num_cells = eqns->mesh->cells.num;

    int num_faces = eqns->mesh->faces.num;
    int num_inner = eqns->mesh->faces.num_inner;
    int off_boundary = eqns->mesh->faces.off_boundary;
    Pair *cell_idx = eqns->mesh->faces.cell_idx;
    VectorPair *weight = eqns->mesh->faces.weight;

    int stride = eqns->primitive.stride;
    double (*primitive)[stride] = primitive_;
    Vector(*gradient)[stride] = gradient_;

    Exchange exchange = equations_exchange(eqns, primitive, stride);

    memset(gradient, 0, num_cells * sizeof(*gradient));

    for (int i = 0; i < num_inner; i++) {
        int left = cell_idx[i].left;
        int right = cell_idx[i].right;
        for (int j = 0; j < stride; j++) {
            double diff = primitive[right][j] - primitive[left][j];
            vector_iadd(&gradient[left][j], vector_mul(diff, weight[i].left));
            vector_iadd(&gradient[right][j], vector_mul(-diff, weight[i].right));
        }
    }

    for (int i = num_inner; i < off_boundary; i++) {
        int left = cell_idx[i].left;
        int right = cell_idx[i].right;
        for (int j = 0; j < stride; j++) {
            double diff = primitive[right][j] - primitive[left][j];
            vector_iadd(&gradient[left][j], vector_mul(diff, weight[i].left));
        }
    }

    equations_exchange_wait_recv(eqns, exchange);

    for (int i = off_boundary; i < num_faces; i++) {
        int left = cell_idx[i].left;
        int right = cell_idx[i].right;
        for (int j = 0; j < stride; j++) {
            double diff = primitive[right][j] - primitive[left][j];
            vector_iadd(&gradient[left][j], vector_mul(diff, weight[i].left));
        }
    }

    equations_exchange_wait_send(eqns, exchange);
}
