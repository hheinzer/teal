#include <math.h>
#include <string.h>

#include "euler.h"
#include "teal/matrix.h"
#include "teal/utils.h"
#include "teal/vector.h"

static void symmetry(void *ghost_, const void *inner_, const void *reference_,
                     const scalar *property, matrix basis)
{
    unused(reference_);
    unused(property);

    Euler *ghost = ghost_;
    const Euler *inner = inner_;
    vector normal = basis.x;

    scalar velocity = vector_dot(inner->velocity, normal);
    ghost->density = inner->density;
    ghost->velocity = vector_sub(inner->velocity, vector_mul(2 * velocity, normal));
    ghost->pressure = inner->pressure;
}

static void supersonic_inflow(void *ghost_, const void *inner_, const void *reference_,
                              const scalar *property, matrix basis)
{
    unused(inner_);
    unused(property);
    unused(basis);

    Euler *ghost = ghost_;
    const Euler *reference = reference_;

    ghost->density = reference->density;
    ghost->velocity = reference->velocity;
    ghost->pressure = reference->pressure;
}

static void supersonic_outflow(void *ghost_, const void *inner_, const void *reference_,
                               const scalar *property, matrix basis)
{
    unused(reference_);
    unused(property);
    unused(basis);

    Euler *ghost = ghost_;
    const Euler *inner = inner_;

    ghost->density = inner->density;
    ghost->velocity = inner->velocity;
    ghost->pressure = inner->pressure;
}

static void subsonic_inflow(void *ghost_, const void *inner_, const void *reference_,
                            const scalar *property, matrix basis)
{
    Euler *ghost = ghost_;
    const Euler *inner = inner_;
    const Euler *reference = reference_;
    scalar gamma = property[0];
    vector normal = basis.x;

    scalar density = inner->density;
    scalar speed_of_sound = sqrt(gamma * inner->pressure / inner->density);
    scalar velocity_diff = vector_dot(vector_sub(reference->velocity, inner->velocity), normal);
    ghost->pressure =
        (reference->pressure + inner->pressure - density * speed_of_sound * velocity_diff) / 2;
    ghost->density =
        reference->density + ((ghost->pressure - reference->pressure) / pow2(speed_of_sound));
    ghost->velocity = vector_sub(
        reference->velocity,
        vector_mul((reference->pressure - ghost->pressure) / (density * speed_of_sound), normal));
}

static void subsonic_outflow(void *ghost_, const void *inner_, const void *reference_,
                             const scalar *property, matrix basis)
{
    Euler *ghost = ghost_;
    const Euler *inner = inner_;
    const Euler *reference = reference_;
    scalar gamma = property[0];
    vector normal = basis.x;

    scalar density = inner->density;
    scalar speed_of_sound = sqrt(gamma * inner->pressure / inner->density);
    ghost->pressure = reference->pressure;
    ghost->density = inner->density + ((ghost->pressure - inner->pressure) / pow2(speed_of_sound));
    ghost->velocity = vector_sub(
        inner->velocity,
        vector_mul((inner->pressure - ghost->pressure) / (density * speed_of_sound), normal));
}

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

static void matvec(scalar res[5], const scalar mat[5][5], const scalar vec[5], scalar factor)
{
    for (number i = 0; i < 5; i++) {
        res[i] = 0;
        for (number j = 0; j < 5; j++) {
            res[i] += mat[i][j] * vec[j];
        }
        res[i] *= factor;
    }
}

static void local_to_global(Euler *ghost, matrix basis)
{
    ghost->velocity = matrix_matvec(matrix_transpose(basis), ghost->velocity);
}

static void farfield(void *ghost_, const void *inner_, const void *reference_,
                     const scalar *property, matrix basis)
{
    Euler *ghost = ghost_;
    scalar gamma = property[0];
    scalar gamma_m1 = gamma - 1;

    Euler inner = global_to_local(inner_, basis);
    *ghost = global_to_local(reference_, basis);

    scalar velocity2_i = vector_dot(inner.velocity, inner.velocity);
    scalar speed_of_sound2_i = gamma * inner.pressure / inner.density;
    scalar speed_of_sound_i = sqrt(speed_of_sound2_i);
    scalar enthalpy_i = (velocity2_i / 2) + (speed_of_sound2_i / gamma_m1);
    scalar eigenvector_inv[5][5] = {
        {
            enthalpy_i + (speed_of_sound_i / gamma_m1 * (inner.velocity.x - speed_of_sound_i)),
            -(inner.velocity.x + (speed_of_sound_i / gamma_m1)),
            -inner.velocity.y,
            -inner.velocity.z,
            1,
        },
        {
            (-2 * enthalpy_i) + (4 / gamma_m1 * speed_of_sound2_i),
            2 * inner.velocity.x,
            2 * inner.velocity.y,
            2 * inner.velocity.z,
            -2,
        },
        {
            -2 * inner.velocity.y * speed_of_sound2_i / gamma_m1,
            0,
            2 * speed_of_sound2_i / gamma_m1,
            0,
            0,
        },
        {
            -2 * inner.velocity.z * speed_of_sound2_i / gamma_m1,
            0,
            0,
            2 * speed_of_sound2_i / gamma_m1,
            0,
        },
        {
            enthalpy_i - (speed_of_sound_i / gamma_m1 * (inner.velocity.x + speed_of_sound_i)),
            -inner.velocity.x + (speed_of_sound_i / gamma_m1),
            -inner.velocity.y,
            -inner.velocity.z,
            1,
        },
    };

    scalar factor = gamma_m1 / (2 * speed_of_sound2_i);
    scalar conserved_i[5] = {
        inner.density, inner.momentum.x, inner.momentum.y, inner.momentum.z, inner.energy,
    };
    scalar conserved_g[5] = {
        ghost->density, ghost->momentum.x, ghost->momentum.y, ghost->momentum.z, ghost->energy,
    };
    scalar characteristic_i[5];
    scalar characteristic_g[5];
    matvec(characteristic_i, eigenvector_inv, conserved_i, factor);
    matvec(characteristic_g, eigenvector_inv, conserved_g, factor);

    scalar speed_of_sound_g = sqrt(gamma * ghost->pressure / ghost->density);
    if (ghost->velocity.x - speed_of_sound_g > 0) {
        characteristic_g[0] = characteristic_i[0];
    }
    if (ghost->velocity.x > 0) {
        characteristic_g[1] = characteristic_i[1];
        characteristic_g[2] = characteristic_i[2];
        characteristic_g[3] = characteristic_i[3];
    }
    if (ghost->velocity.x + speed_of_sound_g > 0) {
        characteristic_g[4] = characteristic_i[4];
    }

    scalar eigenvector[5][5] = {
        {
            1,
            1,
            0,
            0,
            1,
        },
        {
            inner.velocity.x - speed_of_sound_i,
            inner.velocity.x,
            0,
            0,
            inner.velocity.x + speed_of_sound_i,
        },
        {
            inner.velocity.y,
            inner.velocity.y,
            1,
            0,
            inner.velocity.y,
        },
        {
            inner.velocity.z,
            inner.velocity.z,
            0,
            1,
            inner.velocity.z,
        },
        {
            enthalpy_i - (inner.velocity.x * speed_of_sound_i),
            velocity2_i / 2,
            inner.velocity.y,
            inner.velocity.z,
            enthalpy_i + (inner.velocity.x * speed_of_sound_i),

        },
    };
    matvec(conserved_g, eigenvector, characteristic_g, 1);

    ghost->density = conserved_g[0];
    ghost->momentum.x = conserved_g[1];
    ghost->momentum.y = conserved_g[2];
    ghost->momentum.z = conserved_g[3];
    ghost->energy = conserved_g[4];
    euler_primitive(ghost, property);

    local_to_global(ghost, basis);
}

Boundary *euler_boundary(const char *name)
{
    if (!strcmp(name, "symmetry") || !strcmp(name, "slipwall")) {
        return symmetry;
    }
    if (!strcmp(name, "supersonic inflow")) {
        return supersonic_inflow;
    }
    if (!strcmp(name, "supersonic outflow")) {
        return supersonic_outflow;
    }
    if (!strcmp(name, "subsonic inflow")) {
        return subsonic_inflow;
    }
    if (!strcmp(name, "subsonic outflow")) {
        return subsonic_outflow;
    }
    if (!strcmp(name, "farfield")) {
        return farfield;
    }
    error("invalid boundary condition -- '%s'", name);
}
