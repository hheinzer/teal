#include <math.h>
#include <string.h>

#include "../test.h"
#include "equations2.h"
#include "mesh2.h"
#include "sync2.h"
#include "teal2.h"
#include "utils2.h"

static double timestep(const void *primitive, const double *property, double volume,
                       Vector projection)
{
    return 1;
}

static void compute(void *primitive_, const double *property, Vector center, double time,
                    const void *context)
{
    double fac = 2 * M_PI;
    Vector *velocity = primitive_;
    velocity->x = sin(fac * center.x) * cos(fac * center.y) * cos(fac * center.z);
    velocity->y = cos(fac * center.x) * sin(fac * center.y) * cos(fac * center.z);
    velocity->z = 0;
}

static double test_gradient(int num)
{
    Vector min_coord = {0, 0, 0};
    Vector max_coord = {1, 1, 1};
    Triple num_cells = {num, num, num};
    Triple periodic = {1, 1, 1};
    Mesh *mesh = mesh2_create(min_coord, max_coord, num_cells, periodic);
    mesh2_generate(mesh);

    Equations *eqns = equations2_create(mesh, "test", timestep, 0, 0, 0);
    equations2_create_primitive(eqns, (const char *[]){"velocity"}, (int[]){3}, 0, 1);
    equations2_set_initial_condition(eqns, "domain", compute, 0);

    int num_inner = eqns->mesh->cells.num_inner;
    int num_cells_ = eqns->mesh->cells.num;
    double *volume = eqns->mesh->cells.volume;
    double sum_volume = eqns->mesh->cells.sum_volume;
    Vector *center = eqns->mesh->cells.center;

    int stride = eqns->primitive.stride;
    double (*primitive)[stride] = eqns->primitive.data;
    Vector(*gradient)[stride] = teal2_calloc(num_cells_, (int)sizeof(*gradient));

    equations2_gradient(eqns, primitive, gradient);

    double error[stride];
    memset(error, 0, sizeof(error));

    double fac = 2 * M_PI;
    for (int i = 0; i < num_inner; i++) {
        double x = center[i].x, y = center[i].y, z = center[i].z;

        Vector exact[stride];
        exact[0].x = fac * cos(fac * x) * cos(fac * y) * cos(fac * z);
        exact[0].y = -fac * sin(fac * x) * sin(fac * y) * cos(fac * z);
        exact[0].z = -fac * sin(fac * x) * cos(fac * y) * sin(fac * z);
        exact[1].x = -fac * sin(fac * x) * sin(fac * y) * cos(fac * z);
        exact[1].y = fac * cos(fac * x) * cos(fac * y) * cos(fac * z);
        exact[1].z = -fac * cos(fac * x) * sin(fac * y) * sin(fac * z);
        exact[2] = (Vector){0, 0, 0};

        for (int j = 0; j < stride; j++) {
            error[j] += volume[i] * vector2_norm2(vector2_sub(gradient[i][j], exact[j]));
        }
    }

    sync2_sum(error, stride, MPI_DOUBLE);

    double max_error = 0;
    for (int i = 0; i < stride; i++) {
        error[i] = sqrt(error[i] / sum_volume);
        if (error[i] > max_error) {
            max_error = error[i];
        }
    }

    mesh2_destroy(mesh);
    equations2_destroy(eqns);

    teal2_free(gradient);
    return max_error;
}

int main(int argc, char **argv)
{
    teal2_init(&argc, &argv);

    double error16 = test_gradient(16);
    double error32 = test_gradient(32);
    double error64 = test_gradient(64);

    double rate1 = log(error16 / error32) / log(2);
    double rate2 = log(error32 / error64) / log(2);
    test(rate1 > 1.99);
    test(rate2 > 1.99);

    teal2_deinit();
}
