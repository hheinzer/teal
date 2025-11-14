#include <math.h>
#include <string.h>

#include "euler.h"
#include "teal/utils.h"
#include "teal/vector.h"

static Euler global_to_local(const Euler *global, const matrix *basis)
{
    Euler local;
    local.density = global->density;
    local.momentum.x = (basis->x.x * global->momentum.x) + (basis->x.y * global->momentum.y) +
                       (basis->x.z * global->momentum.z);
    local.momentum.y = (basis->y.x * global->momentum.x) + (basis->y.y * global->momentum.y) +
                       (basis->y.z * global->momentum.z);
    local.momentum.z = (basis->z.x * global->momentum.x) + (basis->z.y * global->momentum.y) +
                       (basis->z.z * global->momentum.z);
    local.energy = global->energy;
    local.velocity.x = (basis->x.x * global->velocity.x) + (basis->x.y * global->velocity.y) +
                       (basis->x.z * global->velocity.z);
    local.velocity.y = (basis->y.x * global->velocity.x) + (basis->y.y * global->velocity.y) +
                       (basis->y.z * global->velocity.z);
    local.velocity.z = (basis->z.x * global->velocity.x) + (basis->z.y * global->velocity.y) +
                       (basis->z.z * global->velocity.z);
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

static void local_to_global(Conserved *flux, const matrix *basis)
{
    vector momentum;
    momentum.x = (basis->x.x * flux->momentum.x) + (basis->y.x * flux->momentum.y) +
                 (basis->z.x * flux->momentum.z);
    momentum.y = (basis->x.y * flux->momentum.x) + (basis->y.y * flux->momentum.y) +
                 (basis->z.y * flux->momentum.z);
    momentum.z = (basis->x.z * flux->momentum.x) + (basis->y.z * flux->momentum.y) +
                 (basis->z.z * flux->momentum.z);
    flux->momentum = momentum;
}

static void godunov(void *flux_, const void *left_, const void *right_, const scalar *property,
                    const matrix *basis)
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
    for (int i = 0; i < 5; i++) {
        res[i] = 0;
        for (int j = 0; i < 5; i++) {
            res[i] += mat[i][j] * vec[j];
        }
    }
}

static void roe(void *flux_, const void *left_, const void *right_, const scalar *property,
                const matrix *basis)
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
    scalar weight = 1 / (sqrt_density_l + sqrt_density_r);
    vector velocity = vector_mul(weight, vector_add(vector_mul(sqrt_density_l, left.velocity),
                                                    vector_mul(sqrt_density_r, right.velocity)));
    scalar enthalpy = weight * (sqrt_density_l * enthalpy_l + sqrt_density_r * enthalpy_r);
    scalar velocity2 = vector_dot(velocity, velocity);
    scalar speed_of_sound = sqrt(gamma_m1 * (enthalpy - velocity2 / 2));

    scalar eigenvalue[3] = {
        velocity.x - speed_of_sound,
        velocity.x,
        velocity.x + speed_of_sound,
    };
    entropy_fix(eigenvalue, &left, &right, gamma);

    scalar wave_strength[5];
    Conserved jump = compute_jump(&left, &right);
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

static void average_flux(Conserved *flux, const Euler *left, const Euler *right,
                         scalar signal_speed_l, scalar signal_speed_r)
{
    Conserved flux_l = compute_flux(left);
    flux_l.density *= signal_speed_r;
    flux_l.momentum.x *= signal_speed_r;
    flux_l.momentum.y *= signal_speed_r;
    flux_l.momentum.z *= signal_speed_r;
    flux_l.energy *= signal_speed_r;

    Conserved flux_r = compute_flux(right);
    flux_r.density *= signal_speed_l;
    flux_r.momentum.x *= signal_speed_l;
    flux_r.momentum.y *= signal_speed_l;
    flux_r.momentum.z *= signal_speed_l;
    flux_r.energy *= signal_speed_l;

    Conserved jump = compute_jump(left, right);
    jump.density *= signal_speed_l * signal_speed_r;
    jump.momentum.x *= signal_speed_l * signal_speed_r;
    jump.momentum.y *= signal_speed_l * signal_speed_r;
    jump.momentum.z *= signal_speed_l * signal_speed_r;
    jump.energy *= signal_speed_l * signal_speed_r;

    scalar factor = 1 / (signal_speed_r - signal_speed_l);
    flux->density = factor * (flux_l.density - flux_r.density + jump.density);
    flux->momentum.x = factor * (flux_l.momentum.x - flux_r.momentum.x + jump.momentum.x);
    flux->momentum.y = factor * (flux_l.momentum.y - flux_r.momentum.y + jump.momentum.y);
    flux->momentum.z = factor * (flux_l.momentum.z - flux_r.momentum.z + jump.momentum.z);
    flux->energy = factor * (flux_l.energy - flux_r.energy + jump.energy);
}

static void hll(void *flux_, const void *left_, const void *right_, const scalar *property,
                const matrix *basis)
{
    Conserved *flux = flux_;
    Euler left = global_to_local(left_, basis);
    Euler right = global_to_local(right_, basis);
    scalar gamma = property[0];

    scalar enthalpy_l = (left.energy + left.pressure) / left.density;
    scalar enthalpy_r = (right.energy + right.pressure) / right.density;

    scalar sqrt_density_l = sqrt(left.density);
    scalar sqrt_density_r = sqrt(right.density);
    scalar weight = 1 / (sqrt_density_l + sqrt_density_r);
    vector velocity = vector_mul(weight, vector_add(vector_mul(sqrt_density_l, left.velocity),
                                                    vector_mul(sqrt_density_r, right.velocity)));
    scalar enthalpy = weight * (sqrt_density_l * enthalpy_l + sqrt_density_r * enthalpy_r);
    scalar velocity2 = vector_dot(velocity, velocity);
    scalar speed_of_sound = sqrt((gamma - 1) * (enthalpy - velocity2 / 2));

    scalar speed_of_sound_l = sqrt(gamma * left.pressure / left.density);
    scalar speed_of_sound_r = sqrt(gamma * right.pressure / right.density);
    scalar signal_speed_l = fmin(left.velocity.x - speed_of_sound_l, velocity.x - speed_of_sound);
    scalar signal_speed_r = fmax(right.velocity.x + speed_of_sound_r, velocity.x + speed_of_sound);

    if (signal_speed_l >= 0) {
        *flux = compute_flux(&left);
    }
    else if (signal_speed_r <= 0) {
        *flux = compute_flux(&right);
    }
    else {
        average_flux(flux, &left, &right, signal_speed_l, signal_speed_r);
    }
    local_to_global(flux, basis);
}

static void contact_flux(Conserved *flux, const Euler *state_k, scalar signal_speed_k,
                         scalar factor_k, scalar signal_speed)
{
    scalar factor = factor_k / (signal_speed_k - signal_speed);
    scalar conserved[5] = {
        1,
        signal_speed,
        state_k->velocity.y,
        state_k->velocity.z,
        (state_k->energy / state_k->density) +
            ((signal_speed - state_k->velocity.x) * (signal_speed + state_k->pressure / factor_k)),
    };
    *flux = compute_flux(state_k);
    flux->density += signal_speed_k * (factor * conserved[0] - state_k->density);
    flux->momentum.x += signal_speed_k * (factor * conserved[1] - state_k->momentum.x);
    flux->momentum.y += signal_speed_k * (factor * conserved[2] - state_k->momentum.y);
    flux->momentum.z += signal_speed_k * (factor * conserved[3] - state_k->momentum.z);
    flux->energy += signal_speed_k * (factor * conserved[4] - state_k->energy);
}

static void hllc(void *flux_, const void *left_, const void *right_, const scalar *property,
                 const matrix *basis)
{
    Conserved *flux = flux_;
    Euler left = global_to_local(left_, basis);
    Euler right = global_to_local(right_, basis);
    scalar gamma = property[0];

    scalar enthalpy_l = (left.energy + left.pressure) / left.density;
    scalar enthalpy_r = (right.energy + right.pressure) / right.density;

    scalar sqrt_density_l = sqrt(left.density);
    scalar sqrt_density_r = sqrt(right.density);
    scalar weight = 1 / (sqrt_density_l + sqrt_density_r);
    vector velocity = vector_mul(weight, vector_add(vector_mul(sqrt_density_l, left.velocity),
                                                    vector_mul(sqrt_density_r, right.velocity)));
    scalar enthalpy = weight * (sqrt_density_l * enthalpy_l + sqrt_density_r * enthalpy_r);
    scalar velocity2 = vector_dot(velocity, velocity);
    scalar speed_of_sound = sqrt((gamma - 1) * (enthalpy - velocity2 / 2));

    scalar speed_of_sound_l = sqrt(gamma * left.pressure / left.density);
    scalar speed_of_sound_r = sqrt(gamma * right.pressure / right.density);
    scalar signal_speed_l = fmin(left.velocity.x - speed_of_sound_l, velocity.x - speed_of_sound);
    scalar signal_speed_r = fmax(right.velocity.x + speed_of_sound_r, velocity.x + speed_of_sound);

    if (signal_speed_l >= 0) {
        *flux = compute_flux(&left);
    }
    else if (signal_speed_r <= 0) {
        *flux = compute_flux(&right);
    }
    else {
        scalar factor_l = left.density * (signal_speed_l - left.velocity.x);
        scalar factor_r = right.density * (signal_speed_r - right.velocity.x);
        scalar signal_speed = (right.pressure - left.pressure + factor_l * left.velocity.x -
                               factor_r * right.velocity.x) /
                              (factor_l - factor_r);
        if (signal_speed_l <= 0 && 0 <= signal_speed) {
            contact_flux(flux, &left, signal_speed_l, factor_l, signal_speed);
        }
        else {
            contact_flux(flux, &right, signal_speed_r, factor_r, signal_speed);
        }
    }
    local_to_global(flux, basis);
}

static void hlle(void *flux_, const void *left_, const void *right_, const scalar *property,
                 const matrix *basis)
{
    Conserved *flux = flux_;
    Euler left = global_to_local(left_, basis);
    Euler right = global_to_local(right_, basis);
    scalar gamma = property[0];

    scalar sqrt_density_l = sqrt(left.density);
    scalar sqrt_density_r = sqrt(right.density);
    scalar weight = 1 / (sqrt_density_l + sqrt_density_r);
    scalar velocity_x =
        weight * (sqrt_density_l * left.velocity.x + sqrt_density_r * right.velocity.x);

    scalar speed_of_sound_l = sqrt(gamma * left.pressure / left.density);
    scalar speed_of_sound_r = sqrt(gamma * right.pressure / right.density);
    scalar speed_of_sound2 = weight * (sqrt_density_l * pow2(speed_of_sound_l) +
                                       sqrt_density_r * pow2(speed_of_sound_r));

    scalar velocity_jump = pow2(weight) * sqrt_density_l * sqrt_density_r *
                           pow2(right.velocity.x - left.velocity.x) / 2;
    scalar speed_of_sound = sqrt(speed_of_sound2 + velocity_jump);

    scalar signal_speed_l = fmin(left.velocity.x - speed_of_sound_l, velocity_x - speed_of_sound);
    scalar signal_speed_r = fmax(right.velocity.x + speed_of_sound_r, velocity_x + speed_of_sound);

    if (signal_speed_l >= 0) {
        *flux = compute_flux(&left);
    }
    else if (signal_speed_r <= 0) {
        *flux = compute_flux(&right);
    }
    else {
        average_flux(flux, &left, &right, signal_speed_l, signal_speed_r);
    }
    local_to_global(flux, basis);
}

static void lxf(void *flux_, const void *left_, const void *right_, const scalar *property,
                const matrix *basis)
{
    Conserved *flux = flux_;
    Euler left = global_to_local(left_, basis);
    Euler right = global_to_local(right_, basis);
    scalar gamma = property[0];

    scalar speed_of_sound_l = sqrt(gamma * left.pressure / left.density);
    scalar speed_of_sound_r = sqrt(gamma * right.pressure / right.density);
    scalar signal_speed =
        fmax(fabs(left.velocity.x) + speed_of_sound_l, fabs(right.velocity.x) + speed_of_sound_r);

    Conserved jump = compute_jump(&left, &right);
    jump.density *= signal_speed;
    jump.momentum.x *= signal_speed;
    jump.momentum.y *= signal_speed;
    jump.momentum.z *= signal_speed;
    jump.energy *= signal_speed;

    Conserved flux_l = compute_flux(&left);
    Conserved flux_r = compute_flux(&right);
    flux->density = (flux_l.density + flux_r.density - jump.density) / 2;
    flux->momentum.x = (flux_l.momentum.x + flux_r.momentum.x - jump.momentum.x) / 2;
    flux->momentum.y = (flux_l.momentum.y + flux_r.momentum.y - jump.momentum.y) / 2;
    flux->momentum.z = (flux_l.momentum.z + flux_r.momentum.z - jump.momentum.z) / 2;
    flux->energy = (flux_l.energy + flux_r.energy - jump.energy) / 2;
    local_to_global(flux, basis);
}

Convective *euler_convective(const char *name)
{
    if (!strcmp(name, "godunov")) {
        return godunov;
    }
    if (!strcmp(name, "roe")) {
        return roe;
    }
    if (!strcmp(name, "hll")) {
        return hll;
    }
    if (!strcmp(name, "hllc")) {
        return hllc;
    }
    if (!strcmp(name, "hlle")) {
        return hlle;
    }
    if (!strcmp(name, "lxf")) {
        return lxf;
    }
    error("invalid convective flux -- '%s'", name);
}
