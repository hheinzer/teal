#include <string.h>

#include "euler.h"
#include "riemann.h"
#include "teal/matrix.h"
#include "teal/utils.h"
#include "teal/vector.h"

static Euler global_to_local(const Euler *global, matrix basis)
{
    Euler local;
    local.density = global->density;
    local.momentum = matrix_matvec(basis, global->momentum);
    local.energy = global->energy;
    local.velocity = matrix_matvec(basis, global->velocity);
    local.pressure = global->pressure;
    return local;
}

typedef struct {
    scalar density;
    vector momentum;
    scalar energy;
} Flux;

static void physical(Flux *flux, const Euler *face)
{
    flux->density = face->momentum.x;
    flux->momentum.x = (face->momentum.x * face->velocity.x) + face->pressure;
    flux->momentum.y = (face->momentum.x * face->velocity.y);
    flux->momentum.z = (face->momentum.x * face->velocity.z);
    flux->energy = (face->energy + face->pressure) * face->velocity.x;
}

static Flux local_to_global(const Flux *local, matrix basis)
{
    Flux global;
    global.density = local->density;
    global.momentum = matrix_matvec(matrix_transpose(basis), local->momentum);
    global.energy = local->energy;
    return global;
}

static void godunov(void *flux_, const void *left_, const void *right_, const scalar *property,
                    matrix basis)
{
    Flux *flux = flux_;
    Euler left = global_to_local(left_, basis);
    Euler right = global_to_local(right_, basis);
    scalar gamma = property[0];

    Euler face = riemann(&left, &right, gamma, 0);
    face.momentum = vector_mul(face.density, face.velocity);
    face.energy = (face.pressure / (gamma - 1)) + (vector_dot(face.momentum, face.velocity) / 2);

    physical(flux, &face);
    *flux = local_to_global(flux, basis);
}

// static void roe(void *flux_, const void *left_, const void *right_, const scalar *property,
//                 matrix basis)
//{
// }
//
// static void hll(void *flux_, const void *left_, const void *right_, const scalar *property,
//                 matrix basis)
//{
// }
//
// static void hllc(void *flux_, const void *left_, const void *right_, const scalar *property,
//                  matrix basis)
//{
// }
//
// static void hlle(void *flux_, const void *left_, const void *right_, const scalar *property,
//                  matrix basis)
//{
// }
//
// static void lxf(void *flux_, const void *left_, const void *right_, const scalar *property,
//                 matrix basis)
//{
// }

Convective *euler_convective(const char *name)
{
    if (!strcmp(name, "godunov")) {
        return godunov;
    }
    // if (!strcmp(name, "roe")) {
    //     return roe;
    // }
    // if (!strcmp(name, "hll")) {
    //     return hll;
    // }
    // if (!strcmp(name, "hllc")) {
    //     return hllc;
    // }
    // if (!strcmp(name, "hlle")) {
    //     return hlle;
    // }
    // if (!strcmp(name, "lxf")) {
    //     return lxf;
    // }
    error("invalid convective flux -- '%s'", name);
}
