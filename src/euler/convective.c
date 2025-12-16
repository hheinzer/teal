#include <math.h>
#include <string.h>

#include "euler.h"
#include "mesh/transform.h"
#include "teal/utils.h"
#include "teal/vector.h"

// Rotate global conserved/primitive into a face-aligned basis for 1D Riemann solves.
static Euler global_to_local(const Euler *global, const scalar *property, const Basis *basis)
{
    Euler local;
    local.density = global->density;
    transform_to_local(&local.velocity, basis, &global->velocity);
    local.pressure = global->pressure;
    euler_conserved(&local, property);
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

static void local_to_global(Conserved *flux, const Basis *basis)
{
    vector momentum;
    transform_to_global(&momentum, basis, &flux->momentum);
    flux->momentum = momentum;
}

// Exact Godunov flux via the exact Riemann solver.
static void godunov(void *flux_, const void *left_, const void *right_, const scalar *property,
                    const Basis *basis)
{
    Conserved *flux = flux_;
    Euler left = global_to_local(left_, property, basis);
    Euler right = global_to_local(right_, property, basis);
    scalar gamma = property[EULER_HEAT_CAPACITY_RATIO];

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

    for (long i = 0; i < 3; i++) {
        scalar lambda = eigenvalue[i];
        scalar delta = fmax(0, fmax(lambda - eigenvalue_l[i], eigenvalue_r[i] - lambda));
        if (fabs(lambda) < delta) {
            // eigenvalue[i] = delta;
            eigenvalue[i] = ((sq(lambda) / delta) + delta) / 2;
        }
        else {
            eigenvalue[i] = fabs(lambda);
        }
    }
}

static Conserved compute_jump(const Euler *left, const Euler *right)
{
    Conserved jump;
    jump.density = right->density - left->density;
    jump.momentum.x = right->momentum.x - left->momentum.x;
    jump.momentum.y = right->momentum.y - left->momentum.y;
    jump.momentum.z = right->momentum.z - left->momentum.z;
    jump.energy = right->energy - left->energy;
    return jump;
}

static vector roe_velocity(scalar weight, scalar sqrt_density_l, scalar sqrt_density_r,
                           const Euler *left, const Euler *right)
{
    vector velocity;
    velocity.x = weight * (sqrt_density_l * left->velocity.x + sqrt_density_r * right->velocity.x);
    velocity.y = weight * (sqrt_density_l * left->velocity.y + sqrt_density_r * right->velocity.y);
    velocity.z = weight * (sqrt_density_l * left->velocity.z + sqrt_density_r * right->velocity.z);
    return velocity;
}

// Roe flux with entropy fix on eigenvalues.
static void roe(void *flux_, const void *left_, const void *right_, const scalar *property,
                const Basis *basis)
{
    Conserved *flux = flux_;
    Euler left = global_to_local(left_, property, basis);
    Euler right = global_to_local(right_, property, basis);
    scalar gamma = property[EULER_HEAT_CAPACITY_RATIO];
    scalar gamma_m1 = gamma - 1;

    scalar enthalpy_l = (left.energy + left.pressure) / left.density;
    scalar enthalpy_r = (right.energy + right.pressure) / right.density;

    scalar sqrt_density_l = sqrt(left.density);
    scalar sqrt_density_r = sqrt(right.density);
    scalar weight = 1 / (sqrt_density_l + sqrt_density_r);
    vector velocity = roe_velocity(weight, sqrt_density_l, sqrt_density_r, &left, &right);
    scalar enthalpy = weight * (sqrt_density_l * enthalpy_l + sqrt_density_r * enthalpy_r);

    scalar velocity2 = vector_norm2(velocity);
    scalar speed_of_sound2 = gamma_m1 * (enthalpy - velocity2 / 2);
    scalar speed_of_sound = sqrt(speed_of_sound2);

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
        gamma_m1 / speed_of_sound2 *
        (jump.density * (enthalpy - sq(velocity.x)) + velocity.x * jump.momentum.x - jump.energy +
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

    for (long i = 0; i < 5; i++) {
        dissipation[i] = 0;
        for (long j = 0; j < 5; j++) {
            dissipation[i] += eigenvector[i][j] * characteristic[j];
        }
    }

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

// HLL two-wave solver using Roe-averaged signal speeds.
static void hll(void *flux_, const void *left_, const void *right_, const scalar *property,
                const Basis *basis)
{
    Conserved *flux = flux_;
    Euler left = global_to_local(left_, property, basis);
    Euler right = global_to_local(right_, property, basis);
    scalar gamma = property[EULER_HEAT_CAPACITY_RATIO];

    scalar enthalpy_l = (left.energy + left.pressure) / left.density;
    scalar enthalpy_r = (right.energy + right.pressure) / right.density;

    scalar sqrt_density_l = sqrt(left.density);
    scalar sqrt_density_r = sqrt(right.density);
    scalar weight = 1 / (sqrt_density_l + sqrt_density_r);
    vector velocity = roe_velocity(weight, sqrt_density_l, sqrt_density_r, &left, &right);
    scalar enthalpy = weight * (sqrt_density_l * enthalpy_l + sqrt_density_r * enthalpy_r);

    scalar velocity2 = vector_norm2(velocity);
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

// HLLC with restored contact/tangential waves.
static void hllc(void *flux_, const void *left_, const void *right_, const scalar *property,
                 const Basis *basis)
{
    Conserved *flux = flux_;
    Euler left = global_to_local(left_, property, basis);
    Euler right = global_to_local(right_, property, basis);
    scalar gamma = property[EULER_HEAT_CAPACITY_RATIO];

    scalar enthalpy_l = (left.energy + left.pressure) / left.density;
    scalar enthalpy_r = (right.energy + right.pressure) / right.density;

    scalar sqrt_density_l = sqrt(left.density);
    scalar sqrt_density_r = sqrt(right.density);
    scalar weight = 1 / (sqrt_density_l + sqrt_density_r);
    vector velocity = roe_velocity(weight, sqrt_density_l, sqrt_density_r, &left, &right);
    scalar enthalpy = weight * (sqrt_density_l * enthalpy_l + sqrt_density_r * enthalpy_r);

    scalar velocity2 = vector_norm2(velocity);
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

// HLLE with Einfeldt speeds to avoid carbuncles.
static void hlle(void *flux_, const void *left_, const void *right_, const scalar *property,
                 const Basis *basis)
{
    Conserved *flux = flux_;
    Euler left = global_to_local(left_, property, basis);
    Euler right = global_to_local(right_, property, basis);
    scalar gamma = property[EULER_HEAT_CAPACITY_RATIO];

    scalar sqrt_density_l = sqrt(left.density);
    scalar sqrt_density_r = sqrt(right.density);
    scalar weight = 1 / (sqrt_density_l + sqrt_density_r);
    scalar velocity_x =
        weight * (sqrt_density_l * left.velocity.x + sqrt_density_r * right.velocity.x);

    scalar speed_of_sound2_l = gamma * left.pressure / left.density;
    scalar speed_of_sound2_r = gamma * right.pressure / right.density;

    scalar speed_of_sound_l = sqrt(speed_of_sound2_l);
    scalar speed_of_sound_r = sqrt(speed_of_sound2_r);

    scalar speed_of_sound2 =
        weight * (sqrt_density_l * speed_of_sound2_l + sqrt_density_r * speed_of_sound2_r);

    scalar velocity_jump =
        sq(weight) * sqrt_density_l * sqrt_density_r * sq(right.velocity.x - left.velocity.x) / 2;

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

static void ausmd(void *flux_, const void *left_, const void *right_, const scalar *property,
                  const Basis *basis)
{
    Conserved *flux = flux_;
    Euler left = global_to_local(left_, property, basis);
    Euler right = global_to_local(right_, property, basis);
    scalar gamma = property[EULER_HEAT_CAPACITY_RATIO];

    scalar speed_of_sound_l = sqrt(gamma * left.pressure / left.density);
    scalar speed_of_sound_r = sqrt(gamma * right.pressure / right.density);
    scalar speed_of_sound = fmax(speed_of_sound_l, speed_of_sound_r);

    scalar pressure_over_density_l = left.pressure / left.density;
    scalar pressure_over_density_r = right.pressure / right.density;
    scalar pressure_over_density_sum = pressure_over_density_l + pressure_over_density_r;

    scalar mach_l = left.velocity.x / speed_of_sound;
    scalar mach_r = right.velocity.x / speed_of_sound;

    scalar velocity_plus_l;
    scalar pressure_plus_l;
    if (fabs(left.velocity.x) <= speed_of_sound) {
        scalar alpha_l = 2 * pressure_over_density_l / pressure_over_density_sum;
        scalar factor = left.velocity.x + speed_of_sound;
        velocity_plus_l = (alpha_l * sq(factor) / (4 * speed_of_sound)) +
                          ((1 - alpha_l) * (left.velocity.x + fabs(left.velocity.x)) / 2);
        pressure_plus_l = (left.pressure * sq(factor) / (4 * sq(speed_of_sound))) * (2 - mach_l);
    }
    else {
        velocity_plus_l = (left.velocity.x + fabs(left.velocity.x)) / 2;
        pressure_plus_l =
            (left.pressure * (left.velocity.x + fabs(left.velocity.x)) / 2) / left.velocity.x;
    }

    scalar velocity_minus_r;
    scalar pressure_minus_r;
    if (fabs(right.velocity.x) <= speed_of_sound) {
        scalar alpha_r = 2 * pressure_over_density_r / pressure_over_density_sum;
        scalar factor = right.velocity.x - speed_of_sound;
        velocity_minus_r = (-alpha_r * sq(factor) / (4 * speed_of_sound)) +
                           ((1 - alpha_r) * (right.velocity.x - fabs(right.velocity.x)) / 2);
        pressure_minus_r = (right.pressure * sq(factor) / (4 * sq(speed_of_sound))) * (2 + mach_r);
    }
    else {
        velocity_minus_r = (right.velocity.x - fabs(right.velocity.x)) / 2;
        pressure_minus_r =
            (right.pressure * (right.velocity.x - fabs(right.velocity.x)) / 2) / right.velocity.x;
    }

    scalar mass_flux = (velocity_plus_l * left.density) + (velocity_minus_r * right.density);
    scalar abs_mass_flux = fabs(mass_flux);
    scalar pressure_flux = pressure_plus_l + pressure_minus_r;

    scalar enthalpy_l = (left.energy + left.pressure) / left.density;
    scalar enthalpy_r = (right.energy + right.pressure) / right.density;

    flux->density = mass_flux;
    flux->momentum.x = ((mass_flux * (left.velocity.x + right.velocity.x) -
                         abs_mass_flux * (right.velocity.x - left.velocity.x)) /
                        2) +
                       pressure_flux;
    flux->momentum.y = (mass_flux * (left.velocity.y + right.velocity.y) -
                        abs_mass_flux * (right.velocity.y - left.velocity.y)) /
                       2;
    flux->momentum.z = (mass_flux * (left.velocity.z + right.velocity.z) -
                        abs_mass_flux * (right.velocity.z - left.velocity.z)) /
                       2;
    flux->energy =
        (mass_flux * (enthalpy_l + enthalpy_r) - abs_mass_flux * (enthalpy_r - enthalpy_l)) / 2;
    local_to_global(flux, basis);
}

static void ausmdv(void *flux_, const void *left_, const void *right_, const scalar *property,
                   const Basis *basis)
{
    Conserved *flux = flux_;
    Euler left = global_to_local(left_, property, basis);
    Euler right = global_to_local(right_, property, basis);
    scalar gamma = property[EULER_HEAT_CAPACITY_RATIO];

    scalar speed_of_sound_l = sqrt(gamma * left.pressure / left.density);
    scalar speed_of_sound_r = sqrt(gamma * right.pressure / right.density);
    scalar speed_of_sound = fmax(speed_of_sound_l, speed_of_sound_r);

    scalar pressure_over_density_l = left.pressure / left.density;
    scalar pressure_over_density_r = right.pressure / right.density;
    scalar pressure_over_density_sum = pressure_over_density_l + pressure_over_density_r;

    scalar mach_l = left.velocity.x / speed_of_sound;
    scalar mach_r = right.velocity.x / speed_of_sound;

    scalar velocity_plus_l;
    scalar pressure_plus_l;
    if (fabs(left.velocity.x) <= speed_of_sound) {
        scalar alpha_l = 2 * pressure_over_density_l / pressure_over_density_sum;
        scalar factor = left.velocity.x + speed_of_sound;
        velocity_plus_l = (alpha_l * sq(factor) / (4 * speed_of_sound)) +
                          ((1 - alpha_l) * (left.velocity.x + fabs(left.velocity.x)) / 2);
        pressure_plus_l = (left.pressure * sq(factor) / (4 * sq(speed_of_sound))) * (2 - mach_l);
    }
    else {
        velocity_plus_l = (left.velocity.x + fabs(left.velocity.x)) / 2;
        pressure_plus_l =
            left.pressure * (left.velocity.x + fabs(left.velocity.x)) / (2 * left.velocity.x);
    }

    scalar velocity_minus_r;
    scalar pressure_minus_r;
    if (fabs(right.velocity.x) <= speed_of_sound) {
        scalar alpha_r = 2 * pressure_over_density_r / pressure_over_density_sum;
        scalar factor = right.velocity.x - speed_of_sound;
        velocity_minus_r = (-alpha_r * sq(factor) / (4 * speed_of_sound)) +
                           ((1 - alpha_r) * (right.velocity.x - fabs(right.velocity.x)) / 2);
        pressure_minus_r = (right.pressure * sq(factor) / (4 * sq(speed_of_sound))) * (2 + mach_r);
    }
    else {
        velocity_minus_r = (right.velocity.x - fabs(right.velocity.x)) / 2;
        pressure_minus_r =
            right.pressure * (right.velocity.x - fabs(right.velocity.x)) / (2 * right.velocity.x);
    }

    scalar mass_flux = (velocity_plus_l * left.density) + (velocity_minus_r * right.density);
    scalar abs_mass_flux = fabs(mass_flux);
    scalar pressure_flux = pressure_plus_l + pressure_minus_r;

    scalar normal_momentum_flux_ausmd = (mass_flux * (left.velocity.x + right.velocity.x) -
                                         abs_mass_flux * (right.velocity.x - left.velocity.x)) /
                                        2;

    scalar normal_momentum_flux_ausmv =
        (velocity_plus_l * left.momentum.x) + (velocity_minus_r * right.momentum.x);

    static const scalar switch_parameter = 10;
    scalar pressure_min = fmin(left.pressure, right.pressure);

    scalar switching = fmin(
        1.0 / 2, fmin(1, (switch_parameter * fabs(right.pressure - left.pressure)) / pressure_min));

    scalar normal_momentum_flux = ((1 + switching) * normal_momentum_flux_ausmv +
                                   (1 - switching) * normal_momentum_flux_ausmd) /
                                  2;

    scalar enthalpy_l = (left.energy + left.pressure) / left.density;
    scalar enthalpy_r = (right.energy + right.pressure) / right.density;

    flux->density = mass_flux;
    flux->momentum.x = normal_momentum_flux + pressure_flux;
    flux->momentum.y = (mass_flux * (left.velocity.y + right.velocity.y) -
                        abs_mass_flux * (right.velocity.y - left.velocity.y)) /
                       2;
    flux->momentum.z = (mass_flux * (left.velocity.z + right.velocity.z) -
                        abs_mass_flux * (right.velocity.z - left.velocity.z)) /
                       2;
    flux->energy =
        (mass_flux * (enthalpy_l + enthalpy_r) - abs_mass_flux * (enthalpy_r - enthalpy_l)) / 2;
    local_to_global(flux, basis);
}

// Lax-Friedrichs/Rusanov with max signal speed estimate.
static void lxf(void *flux_, const void *left_, const void *right_, const scalar *property,
                const Basis *basis)
{
    Conserved *flux = flux_;
    Euler left = global_to_local(left_, property, basis);
    Euler right = global_to_local(right_, property, basis);
    scalar gamma = property[EULER_HEAT_CAPACITY_RATIO];

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
    if (!strcmp(name, "ausmd")) {
        return ausmd;
    }
    if (!strcmp(name, "ausmdv")) {
        return ausmdv;
    }
    if (!strcmp(name, "lxf")) {
        return lxf;
    }
    error("invalid convective flux (%s)", name);
}
