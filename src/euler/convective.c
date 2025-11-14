#include <math.h>
#include <string.h>

#include "euler.h"
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
} Conserved;

static Conserved compute_flux(const Euler *local)
{
    Conserved flux;
    flux.density = local->momentum.x;
    flux.momentum.x = (local->momentum.x * local->velocity.x) + local->pressure;
    flux.momentum.y = (local->momentum.x * local->velocity.y);
    flux.momentum.z = (local->momentum.x * local->velocity.z);
    flux.energy = (local->energy + local->pressure) * local->velocity.x;
    return flux;
}

static void local_to_global(Conserved *flux, matrix basis)
{
    flux->momentum = matrix_matvec(matrix_transpose(basis), flux->momentum);
}

static void godunov(void *flux_, const void *left_, const void *right_, const scalar *property,
                    matrix basis)
{
    Conserved *flux = flux_;
    Euler left = global_to_local(left_, basis);
    Euler right = global_to_local(right_, basis);
    scalar gamma = property[0];

    Euler face = euler_riemann(&left, &right, gamma, 0);
    euler_conserved(&face, property);

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
