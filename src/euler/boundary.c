#include <math.h>
#include <string.h>

#include "euler.h"
#include "teal.h"

static void symmetry(void *outer_, const void *inner_, const double *property, Matrix basis,
                     const void *context_)
{
    (void)property;
    (void)context_;
    EulerPrimitive *outer = outer_;
    const EulerPrimitive *inner = inner_;

    double velocity = vector_dot(inner->velocity, basis.x);

    outer->density = inner->density;
    outer->velocity = vector_sub(inner->velocity, vector_mul(2 * velocity, basis.x));
    outer->pressure = inner->pressure;
}

static void supersonic_inflow(void *outer_, const void *inner_, const double *property,
                              Matrix basis, const void *context_)
{
    (void)inner_;
    (void)property;
    (void)basis;
    EulerPrimitive *outer = outer_;
    const EulerPrimitive *context = context_;

    *outer = *context;
}

static void supersonic_outflow(void *outer_, const void *inner_, const double *property,
                               Matrix basis, const void *context_)
{
    (void)property;
    (void)basis;
    (void)context_;
    EulerPrimitive *outer = outer_;
    const EulerPrimitive *inner = inner_;

    *outer = *inner;
}

static void subsonic_inflow(void *outer_, const void *inner_, const double *property, Matrix basis,
                            const void *context_)
{
    EulerPrimitive *outer = outer_;
    const EulerPrimitive *inner = inner_;
    const EulerPrimitive *context = context_;
    double gamma = property[EULER_HEAT_CAPACITY_RATIO];

    double density = inner->density;
    double speed_of_sound2 = gamma * inner->pressure / inner->density;
    double speed_of_sound = sqrt(speed_of_sound2);
    double velocity = vector_dot(vector_sub(context->velocity, inner->velocity), basis.x);

    double factor1 = density * speed_of_sound;
    outer->pressure = (context->pressure + inner->pressure - (factor1 * velocity)) / 2;
    outer->density = context->density + ((outer->pressure - context->pressure) / speed_of_sound2);

    double factor2 = (context->pressure - outer->pressure) / factor1;
    outer->velocity = vector_sub(context->velocity, vector_mul(factor2, basis.x));
}

static void subsonic_outflow(void *outer_, const void *inner_, const double *property, Matrix basis,
                             const void *context_)
{
    EulerPrimitive *outer = outer_;
    const EulerPrimitive *inner = inner_;
    const EulerPrimitive *context = context_;
    double gamma = property[EULER_HEAT_CAPACITY_RATIO];

    double density = inner->density;
    double speed_of_sound2 = gamma * inner->pressure / inner->density;
    double speed_of_sound = sqrt(speed_of_sound2);

    outer->pressure = context->pressure;
    outer->density = inner->density + ((outer->pressure - inner->pressure) / speed_of_sound2);

    double factor = (inner->pressure - outer->pressure) / (density * speed_of_sound);
    outer->velocity = vector_sub(inner->velocity, vector_mul(factor, basis.x));
}

static void pressure_outflow(void *outer_, const void *inner_, const double *property, Matrix basis,
                             const void *context_)
{
    EulerPrimitive *outer = outer_;
    const EulerPrimitive *inner = inner_;
    const EulerPrimitive *context = context_;
    double gamma = property[EULER_HEAT_CAPACITY_RATIO];

    double speed_of_sound = sqrt(gamma * inner->pressure / inner->density);
    double velocity = vector_dot(inner->velocity, basis.x);
    double pressure = (velocity / speed_of_sound < 1) ? context->pressure : inner->pressure;

    outer->density = inner->density * pressure / inner->pressure;
    outer->velocity = inner->velocity;
    outer->pressure = pressure;
}

typedef struct {
    double density;
    Vector velocity;
    double pressure;
    Vector momentum;
    double energy;
} Euler;

static Euler global_to_local(const EulerPrimitive *global, const double *property, Matrix basis)
{
    Euler local;
    local.density = global->density;
    local.velocity = matrix_vector(basis, global->velocity);
    local.pressure = global->pressure;
    local.momentum = vector_mul(local.density, local.velocity);
    local.energy = (local.pressure / (property[EULER_HEAT_CAPACITY_RATIO] - 1)) +
                   (vector_dot(local.momentum, local.velocity) / 2);
    return local;
}

static void local_to_global(EulerPrimitive *global, Euler local, const double *property,
                            Matrix basis)
{
    local.velocity = vector_div(local.momentum, local.density);
    local.pressure = (property[EULER_HEAT_CAPACITY_RATIO] - 1) *
                     (local.energy - (vector_dot(local.momentum, local.velocity) / 2));
    global->density = fmax(1e-8, local.density);
    global->velocity = matrix_vector(matrix_transpose(basis), local.velocity);
    global->pressure = fmax(1e-8, local.pressure);
}

static void farfield(void *outer_, const void *inner_, const double *property, Matrix basis,
                     const void *context_)
{
    Euler outer = global_to_local(context_, property, basis);
    Euler inner = global_to_local(inner_, property, basis);
    double gamma = property[EULER_HEAT_CAPACITY_RATIO];
    double gamma_m1 = gamma - 1;

    double velocity2_i = vector_norm2(inner.velocity);
    double speed_of_sound2_i = gamma * inner.pressure / inner.density;
    double speed_of_sound_i = sqrt(speed_of_sound2_i);
    double enthalpy_i = (velocity2_i / 2) + (speed_of_sound2_i / gamma_m1);

    double eigenvector_inv[5][5] = {
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

    double factor = gamma_m1 / (2 * speed_of_sound2_i);
    double conserved_i[5] = {
        inner.density, inner.momentum.x, inner.momentum.y, inner.momentum.z, inner.energy,
    };
    double conserved_g[5] = {
        outer.density, outer.momentum.x, outer.momentum.y, outer.momentum.z, outer.energy,
    };
    double characteristic_i[5];
    double characteristic_g[5];

    for (long i = 0; i < 5; i++) {
        characteristic_i[i] = 0;
        characteristic_g[i] = 0;
        for (long j = 0; j < 5; j++) {
            characteristic_i[i] += eigenvector_inv[i][j] * conserved_i[j];
            characteristic_g[i] += eigenvector_inv[i][j] * conserved_g[j];
        }
        characteristic_i[i] *= factor;
        characteristic_g[i] *= factor;
    }

    double speed_of_sound_g = sqrt(gamma * outer.pressure / outer.density);
    if (outer.velocity.x - speed_of_sound_g > 0) {
        characteristic_g[0] = characteristic_i[0];
    }
    if (outer.velocity.x > 0) {
        characteristic_g[1] = characteristic_i[1];
        characteristic_g[2] = characteristic_i[2];
        characteristic_g[3] = characteristic_i[3];
    }
    if (outer.velocity.x + speed_of_sound_g > 0) {
        characteristic_g[4] = characteristic_i[4];
    }

    double eigenvector[5][5] = {
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

    for (long i = 0; i < 5; i++) {
        conserved_g[i] = 0;
        for (long j = 0; j < 5; j++) {
            conserved_g[i] += eigenvector[i][j] * characteristic_g[j];
        }
    }

    outer.density = conserved_g[0];
    outer.momentum.x = conserved_g[1];
    outer.momentum.y = conserved_g[2];
    outer.momentum.z = conserved_g[3];
    outer.energy = conserved_g[4];
    local_to_global(outer_, outer, property, basis);
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
    if (!strcmp(name, "pressure outflow")) {
        return pressure_outflow;
    }
    if (!strcmp(name, "farfield")) {
        return farfield;
    }
    teal_error("invalid boundary condition (%s)", name);
}
