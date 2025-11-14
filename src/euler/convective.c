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

    *flux = compute_flux(&face);
    local_to_global(flux, basis);
}

static void entropy_fix(scalar eigenvalue[3], const Euler *left, const Euler *right, scalar gamma)
{
    scalar speed_of_sound_l = sqrt(gamma * left->pressure / left->density);
    scalar speed_of_sound_r = sqrt(gamma * right->pressure / right->density);
    scalar eigenvalue_l[3] = {
        left->velocity.x - speed_of_sound_l,
        left->velocity.x,
        left->velocity.x + speed_of_sound_l,
    };
    scalar eigenvalue_r[3] = {
        right->velocity.x - speed_of_sound_r,
        right->velocity.x,
        right->velocity.x + speed_of_sound_r,
    };
    for (number i = 0; i < 3; i++) {
        scalar delta =
            fmax(0, fmax(eigenvalue[i] - eigenvalue_l[i], eigenvalue_r[i] - eigenvalue[i]));
        if (fabs(eigenvalue[i]) < delta) {
            eigenvalue[i] = ((pow2(eigenvalue[i]) / delta) + delta) / 2;
        }
        else {
            eigenvalue[i] = fabs(eigenvalue[i]);
        }
    }
}

static Conserved compute_jump(const Euler *left, const Euler *right)
{
    Conserved jump;
    jump.density = right->density - left->density;
    jump.momentum = vector_sub(right->momentum, left->momentum);
    jump.energy = right->energy - left->energy;
    return jump;
}

static void matvec(scalar res[5], const scalar mat[5][5], const scalar vec[5])
{
    for (number i = 0; i < 5; i++) {
        res[i] = 0;
        for (number j = 0; j < 5; j++) {
            res[i] += mat[i][j] * vec[j];
        }
    }
}

static void roe(void *flux_, const void *left_, const void *right_, const scalar *property,
                matrix basis)
{
    Conserved *flux = flux_;
    Euler left = global_to_local(left_, basis);
    Euler right = global_to_local(right_, basis);
    scalar gamma = property[0];
    scalar gamma_m1 = gamma - 1;

    scalar enthalpy_l = (left.energy + left.pressure) / left.density;
    scalar enthalpy_r = (right.energy + right.pressure) / right.density;

    scalar sqrt_density_l = sqrt(left.density);
    scalar sqrt_density_r = sqrt(right.density);
    scalar factor = 1 / (sqrt_density_l + sqrt_density_r);
    vector velocity = vector_mul(factor, vector_add(vector_mul(sqrt_density_l, left.velocity),
                                                    vector_mul(sqrt_density_r, right.velocity)));
    scalar enthalpy = factor * (sqrt_density_l * enthalpy_l + enthalpy_r * sqrt_density_r);
    scalar velocity2 = vector_dot(velocity, velocity);
    scalar speed_of_sound = sqrt(gamma_m1 * (enthalpy - velocity2 / 2));

    scalar eigenvalue[3] = {
        velocity.x - speed_of_sound,
        velocity.x,
        velocity.x + speed_of_sound,
    };
    entropy_fix(eigenvalue, &left, &right, gamma);

    Conserved jump = compute_jump(&left, &right);
    scalar wave_strength[5];
    wave_strength[2] = jump.momentum.y - (velocity.y * jump.density);
    wave_strength[3] = jump.momentum.z - (velocity.z * jump.density);
    wave_strength[1] =
        gamma_m1 / pow2(speed_of_sound) *
        (jump.density * (enthalpy - pow2(velocity.x)) + velocity.x * jump.momentum.x - jump.energy +
         wave_strength[2] * velocity.y + wave_strength[3] * velocity.z);
    wave_strength[0] = (jump.density * (velocity.x + speed_of_sound) - jump.momentum.x -
                        speed_of_sound * wave_strength[1]) /
                       (2 * speed_of_sound);
    wave_strength[4] = jump.density - (wave_strength[0] + wave_strength[1]);

    scalar dissipation[5];
    scalar eigenvector[5][5] = {
        {
            1,
            1,
            0,
            0,
            1,
        },
        {
            velocity.x - speed_of_sound,
            velocity.x,
            0,
            0,
            velocity.x + speed_of_sound,
        },
        {
            velocity.y,
            velocity.y,
            1,
            0,
            velocity.y,
        },
        {
            velocity.z,
            velocity.z,
            0,
            1,
            velocity.z,
        },
        {
            enthalpy - (velocity.x * speed_of_sound),
            velocity2 / 2,
            velocity.y,
            velocity.z,
            enthalpy + (velocity.x * speed_of_sound),
        },
    };
    scalar characteristic[5] = {
        wave_strength[0] * eigenvalue[0], wave_strength[1] * eigenvalue[1],
        wave_strength[2] * eigenvalue[1], wave_strength[3] * eigenvalue[1],
        wave_strength[4] * eigenvalue[2],
    };
    matvec(dissipation, eigenvector, characteristic);

    Conserved flux_l = compute_flux(&left);
    Conserved flux_r = compute_flux(&right);
    flux->density = (flux_l.density + flux_r.density - dissipation[0]) / 2;
    flux->momentum.x = (flux_l.momentum.x + flux_r.momentum.x - dissipation[1]) / 2;
    flux->momentum.y = (flux_l.momentum.y + flux_r.momentum.y - dissipation[2]) / 2;
    flux->momentum.z = (flux_l.momentum.z + flux_r.momentum.z - dissipation[3]) / 2;
    flux->energy = (flux_l.energy + flux_r.energy - dissipation[4]) / 2;
    local_to_global(flux, basis);
}

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
    if (!strcmp(name, "roe")) {
        return roe;
    }
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
