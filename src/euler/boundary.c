#include <math.h>
#include <string.h>

#include "euler.h"
#include "teal/utils.h"

static void symmetry(void *ghost_, const void *inner_, const void *reference_,
                     const scalar *property, const matrix *basis)
{
    unused(reference_);
    unused(property);

    Euler *ghost = ghost_;
    const Euler *inner = inner_;
    vector normal = basis->x;

    scalar velocity = (inner->velocity.x * normal.x) + (inner->velocity.y * normal.y) +
                      (inner->velocity.z * normal.z);

    ghost->density = inner->density;
    ghost->velocity.x = inner->velocity.x - (2 * velocity * normal.x);
    ghost->velocity.y = inner->velocity.y - (2 * velocity * normal.y);
    ghost->velocity.z = inner->velocity.z - (2 * velocity * normal.z);
    ghost->pressure = inner->pressure;
}

static void supersonic_inflow(void *ghost_, const void *inner_, const void *reference_,
                              const scalar *property, const matrix *basis)
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
                               const scalar *property, const matrix *basis)
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
                            const scalar *property, const matrix *basis)
{
    Euler *ghost = ghost_;
    const Euler *inner = inner_;
    const Euler *reference = reference_;
    scalar gamma = property[0];
    vector normal = basis->x;

    scalar density = inner->density;
    scalar speed_of_sound2 = gamma * inner->pressure / inner->density;
    scalar speed_of_sound = sqrt(speed_of_sound2);

    scalar velocity = ((reference->velocity.x - inner->velocity.x) * normal.x) +
                      ((reference->velocity.y - inner->velocity.y) * normal.y) +
                      ((reference->velocity.z - inner->velocity.z) * normal.z);

    scalar factor1 = density * speed_of_sound;
    ghost->pressure = (reference->pressure + inner->pressure - factor1 * velocity) / 2;
    ghost->density =
        reference->density + ((ghost->pressure - reference->pressure) / speed_of_sound2);

    scalar factor2 = (reference->pressure - ghost->pressure) / factor1;
    ghost->velocity.x = reference->velocity.x - (factor2 * normal.x);
    ghost->velocity.y = reference->velocity.y - (factor2 * normal.y);
    ghost->velocity.z = reference->velocity.z - (factor2 * normal.z);
}

static void subsonic_outflow(void *ghost_, const void *inner_, const void *reference_,
                             const scalar *property, const matrix *basis)
{
    Euler *ghost = ghost_;
    const Euler *inner = inner_;
    const Euler *reference = reference_;
    scalar gamma = property[0];
    vector normal = basis->x;

    scalar density = inner->density;
    scalar speed_of_sound2 = gamma * inner->pressure / inner->density;
    scalar speed_of_sound = sqrt(speed_of_sound2);

    ghost->pressure = reference->pressure;
    ghost->density = inner->density + ((ghost->pressure - inner->pressure) / speed_of_sound2);

    scalar factor = (inner->pressure - ghost->pressure) / (density * speed_of_sound);
    ghost->velocity.x = inner->velocity.x - (factor * normal.x);
    ghost->velocity.y = inner->velocity.y - (factor * normal.y);
    ghost->velocity.z = inner->velocity.z - (factor * normal.z);
}

static void matvec(vector *res, const matrix *mat, const vector *vec)
{
    res->x = (mat->x.x * vec->x) + (mat->x.y * vec->y) + (mat->x.z * vec->z);
    res->y = (mat->y.x * vec->x) + (mat->y.y * vec->y) + (mat->y.z * vec->z);
    res->z = (mat->z.x * vec->x) + (mat->z.y * vec->y) + (mat->z.z * vec->z);
}

static Euler global_to_local(const Euler *global, const matrix *basis)
{
    Euler local;
    local.density = global->density;
    matvec(&local.momentum, basis, &global->momentum);
    local.energy = global->energy;
    matvec(&local.velocity, basis, &global->velocity);
    local.pressure = global->pressure;
    return local;
}

static void matvec_t(vector *res, const matrix *mat, const vector *vec)
{
    res->x = (mat->x.x * vec->x) + (mat->y.x * vec->y) + (mat->z.x * vec->z);
    res->y = (mat->x.y * vec->x) + (mat->y.y * vec->y) + (mat->z.y * vec->z);
    res->z = (mat->x.z * vec->x) + (mat->y.z * vec->y) + (mat->z.z * vec->z);
}

static void local_to_global(Euler *ghost, const matrix *basis)
{
    vector velocity;
    matvec_t(&velocity, basis, &ghost->velocity);
    ghost->velocity = velocity;
}

static void farfield(void *ghost_, const void *inner_, const void *reference_,
                     const scalar *property, const matrix *basis)
{
    Euler *ghost = ghost_;
    scalar gamma = property[0];
    scalar gamma_m1 = gamma - 1;

    Euler inner = global_to_local(inner_, basis);
    *ghost = global_to_local(reference_, basis);

    scalar velocity2_i = pow2(inner.velocity.x) + pow2(inner.velocity.y) + pow2(inner.velocity.z);

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

    for (int i = 0; i < 5; i++) {
        characteristic_i[i] = 0;
        characteristic_g[i] = 0;
        for (int j = 0; j < 5; j++) {
            characteristic_i[i] += eigenvector_inv[i][j] * conserved_i[j];
            characteristic_g[i] += eigenvector_inv[i][j] * conserved_g[j];
        }
        characteristic_i[i] *= factor;
        characteristic_g[i] *= factor;
    }

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

    for (int i = 0; i < 5; i++) {
        conserved_g[i] = 0;
        for (int j = 0; j < 5; j++) {
            conserved_g[i] += eigenvector[i][j] * characteristic_g[j];
        }
    }

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
