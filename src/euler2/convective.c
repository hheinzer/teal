#include <math.h>
#include <string.h>

#include "euler2.h"
#include "teal2.h"

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
    local.momentum = vector2_mul(local.density, local.velocity);
    local.energy = (local.pressure / (property[EULER_HEAT_CAPACITY_RATIO] - 1)) +
                   (vector2_dot(local.momentum, local.velocity) / 2);
    return local;
}

static void local_to_global(EulerConserved *flux, Matrix basis)
{
    flux->momentum = matrix_vector(matrix_transpose(basis), flux->momentum);
}

static EulerConserved physical_flux(Euler local)
{
    EulerConserved flux;
    flux.density = local.momentum.x;
    flux.momentum.x = (local.momentum.x * local.velocity.x) + local.pressure;
    flux.momentum.y = (local.momentum.x * local.velocity.y);
    flux.momentum.z = (local.momentum.x * local.velocity.z);
    flux.energy = (local.energy + local.pressure) * local.velocity.x;
    return flux;
}

static void godunov(void *flux_, const void *left_, const void *right_, const double *property,
                    Matrix basis)
{
    EulerConserved *flux = flux_;
    Euler left = global_to_local(left_, property, basis);
    Euler right = global_to_local(right_, property, basis);
    double gamma = property[EULER_HEAT_CAPACITY_RATIO];

    EulerPrimitive left_p = {left.density, left.velocity, left.pressure};
    EulerPrimitive right_p = {right.density, right.velocity, right.pressure};
    EulerPrimitive face_p = euler2_riemann(left_p, right_p, gamma, 0);

    Euler face;
    face.density = face_p.density;
    face.velocity = face_p.velocity;
    face.pressure = face_p.pressure;
    face.momentum = vector2_mul(face.density, face.velocity);
    face.energy = (face.pressure / (gamma - 1)) + (vector2_dot(face.momentum, face.velocity) / 2);

    *flux = physical_flux(face);
    local_to_global(flux, basis);
}

static Vector roe_velocity(double weight, double sqrt_density_l, double sqrt_density_r, Euler left,
                           Euler right)
{
    return vector2_mul(weight, vector2_add(vector2_mul(sqrt_density_l, left.velocity),
                                           vector2_mul(sqrt_density_r, right.velocity)));
}

static void entropy_fix(double eigenvalue[3], Euler left, Euler right, double gamma)
{
    double speed_of_sound_l = sqrt(gamma * left.pressure / left.density);
    double speed_of_sound_r = sqrt(gamma * right.pressure / right.density);

    double eigenvalue_l[3] = {
        left.velocity.x - speed_of_sound_l,
        left.velocity.x,
        left.velocity.x + speed_of_sound_l,
    };
    double eigenvalue_r[3] = {
        right.velocity.x - speed_of_sound_r,
        right.velocity.x,
        right.velocity.x + speed_of_sound_r,
    };

    for (long i = 0; i < 3; i++) {
        double lambda = eigenvalue[i];
        double delta = fmax(0, fmax(lambda - eigenvalue_l[i], eigenvalue_r[i] - lambda));
        if (fabs(lambda) < delta) {
            // eigenvalue[i] = delta; // first entropy fix of Harten and Hyman
            eigenvalue[i] =
                ((sq(lambda) / delta) + delta) / 2;  // second entropy fix of Harten and Hyman
        }
        else {
            eigenvalue[i] = fabs(lambda);
        }
    }
}

static EulerConserved compute_jump(Euler left, Euler right)
{
    EulerConserved jump;
    jump.density = right.density - left.density;
    jump.momentum = vector2_sub(right.momentum, left.momentum);
    jump.energy = right.energy - left.energy;
    return jump;
}

static void roe(void *flux_, const void *left_, const void *right_, const double *property,
                Matrix basis)
{
    EulerConserved *flux = flux_;
    Euler left = global_to_local(left_, property, basis);
    Euler right = global_to_local(right_, property, basis);
    double gamma = property[EULER_HEAT_CAPACITY_RATIO];
    double gamma_m1 = gamma - 1;

    double enthalpy_l = (left.energy + left.pressure) / left.density;
    double enthalpy_r = (right.energy + right.pressure) / right.density;

    double sqrt_density_l = sqrt(left.density);
    double sqrt_density_r = sqrt(right.density);
    double weight = 1 / (sqrt_density_l + sqrt_density_r);
    Vector velocity = roe_velocity(weight, sqrt_density_l, sqrt_density_r, left, right);
    double enthalpy = weight * ((sqrt_density_l * enthalpy_l) + (sqrt_density_r * enthalpy_r));

    double velocity2 = vector2_norm2(velocity);
    double speed_of_sound2 = gamma_m1 * (enthalpy - (velocity2 / 2));
    double speed_of_sound = sqrt(speed_of_sound2);

    double eigenvalue[3] = {
        velocity.x - speed_of_sound,
        velocity.x,
        velocity.x + speed_of_sound,
    };
    entropy_fix(eigenvalue, left, right, gamma);

    double wave_strength[5];
    EulerConserved jump = compute_jump(left, right);
    wave_strength[2] = jump.momentum.y - (velocity.y * jump.density);
    wave_strength[3] = jump.momentum.z - (velocity.z * jump.density);
    wave_strength[1] =
        gamma_m1 / speed_of_sound2 *
        ((jump.density * (enthalpy - sq(velocity.x))) + (velocity.x * jump.momentum.x) -
         jump.energy + (wave_strength[2] * velocity.y) + (wave_strength[3] * velocity.z));
    wave_strength[0] = ((jump.density * (velocity.x + speed_of_sound)) - jump.momentum.x -
                        (speed_of_sound * wave_strength[1])) /
                       (2 * speed_of_sound);
    wave_strength[4] = jump.density - (wave_strength[0] + wave_strength[1]);

    double dissipation[5];
    double eigenvector[5][5] = {
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
    double characteristic[5] = {
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

    EulerConserved flux_l = physical_flux(left);
    EulerConserved flux_r = physical_flux(right);
    flux->density = (flux_l.density + flux_r.density - dissipation[0]) / 2;
    flux->momentum.x = (flux_l.momentum.x + flux_r.momentum.x - dissipation[1]) / 2;
    flux->momentum.y = (flux_l.momentum.y + flux_r.momentum.y - dissipation[2]) / 2;
    flux->momentum.z = (flux_l.momentum.z + flux_r.momentum.z - dissipation[3]) / 2;
    flux->energy = (flux_l.energy + flux_r.energy - dissipation[4]) / 2;
    local_to_global(flux, basis);
}

static void average_flux(EulerConserved *flux, Euler left, Euler right, double signal_speed_l,
                         double signal_speed_r)
{
    EulerConserved flux_l = physical_flux(left);
    flux_l.density *= signal_speed_r;
    vector2_imul(&flux_l.momentum, signal_speed_r);
    flux_l.energy *= signal_speed_r;

    EulerConserved flux_r = physical_flux(right);
    flux_r.density *= signal_speed_l;
    vector2_imul(&flux_r.momentum, signal_speed_l);
    flux_r.energy *= signal_speed_l;

    double signal_speed_j = signal_speed_l * signal_speed_r;
    EulerConserved jump = compute_jump(left, right);
    jump.density *= signal_speed_j;
    vector2_imul(&jump.momentum, signal_speed_j);
    jump.energy *= signal_speed_j;

    double factor = 1 / (signal_speed_r - signal_speed_l);
    flux->density = factor * (flux_l.density - flux_r.density + jump.density);
    flux->momentum = vector2_mul(
        factor, vector2_add(vector2_sub(flux_l.momentum, flux_r.momentum), jump.momentum));
    flux->energy = factor * (flux_l.energy - flux_r.energy + jump.energy);
}

static void hll(void *flux_, const void *left_, const void *right_, const double *property,
                Matrix basis)
{
    EulerConserved *flux = flux_;
    Euler left = global_to_local(left_, property, basis);
    Euler right = global_to_local(right_, property, basis);
    double gamma = property[EULER_HEAT_CAPACITY_RATIO];

    double enthalpy_l = (left.energy + left.pressure) / left.density;
    double enthalpy_r = (right.energy + right.pressure) / right.density;

    double sqrt_density_l = sqrt(left.density);
    double sqrt_density_r = sqrt(right.density);
    double weight = 1 / (sqrt_density_l + sqrt_density_r);
    Vector velocity = roe_velocity(weight, sqrt_density_l, sqrt_density_r, left, right);
    double enthalpy = weight * ((sqrt_density_l * enthalpy_l) + (sqrt_density_r * enthalpy_r));

    double velocity2 = vector2_norm2(velocity);
    double speed_of_sound = sqrt((gamma - 1) * (enthalpy - (velocity2 / 2)));

    double speed_of_sound_l = sqrt(gamma * left.pressure / left.density);
    double speed_of_sound_r = sqrt(gamma * right.pressure / right.density);

    double signal_speed_l = fmin(left.velocity.x - speed_of_sound_l, velocity.x - speed_of_sound);
    double signal_speed_r = fmax(right.velocity.x + speed_of_sound_r, velocity.x + speed_of_sound);

    if (signal_speed_l >= 0) {
        *flux = physical_flux(left);
    }
    else if (signal_speed_r <= 0) {
        *flux = physical_flux(right);
    }
    else {
        average_flux(flux, left, right, signal_speed_l, signal_speed_r);
    }
    local_to_global(flux, basis);
}

static void contact_flux(EulerConserved *flux, Euler state_k, double signal_speed_k,
                         double factor_k, double signal_speed)
{
    double factor = factor_k / (signal_speed_k - signal_speed);
    EulerConserved star;
    star.density = factor;
    star.momentum.x = factor * signal_speed;
    star.momentum.y = factor * state_k.velocity.y;
    star.momentum.z = factor * state_k.velocity.z;
    star.energy =
        factor *
        ((state_k.energy / state_k.density) +
         ((signal_speed - state_k.velocity.x) * (signal_speed + (state_k.pressure / factor_k))));

    *flux = physical_flux(state_k);
    flux->density += signal_speed_k * (star.density - state_k.density);
    vector2_iadd(&flux->momentum,
                 vector2_mul(signal_speed_k, vector2_sub(star.momentum, state_k.momentum)));
    flux->energy += signal_speed_k * (star.energy - state_k.energy);
}

static void hllc(void *flux_, const void *left_, const void *right_, const double *property,
                 Matrix basis)
{
    EulerConserved *flux = flux_;
    Euler left = global_to_local(left_, property, basis);
    Euler right = global_to_local(right_, property, basis);
    double gamma = property[EULER_HEAT_CAPACITY_RATIO];

    double enthalpy_l = (left.energy + left.pressure) / left.density;
    double enthalpy_r = (right.energy + right.pressure) / right.density;

    double sqrt_density_l = sqrt(left.density);
    double sqrt_density_r = sqrt(right.density);
    double weight = 1 / (sqrt_density_l + sqrt_density_r);
    Vector velocity = roe_velocity(weight, sqrt_density_l, sqrt_density_r, left, right);
    double enthalpy = weight * ((sqrt_density_l * enthalpy_l) + (sqrt_density_r * enthalpy_r));

    double velocity2 = vector2_norm2(velocity);
    double speed_of_sound = sqrt((gamma - 1) * (enthalpy - (velocity2 / 2)));

    double speed_of_sound_l = sqrt(gamma * left.pressure / left.density);
    double speed_of_sound_r = sqrt(gamma * right.pressure / right.density);

    double signal_speed_l = fmin(left.velocity.x - speed_of_sound_l, velocity.x - speed_of_sound);
    double signal_speed_r = fmax(right.velocity.x + speed_of_sound_r, velocity.x + speed_of_sound);

    if (signal_speed_l >= 0) {
        *flux = physical_flux(left);
    }
    else if (signal_speed_r <= 0) {
        *flux = physical_flux(right);
    }
    else {
        double factor_l = left.density * (signal_speed_l - left.velocity.x);
        double factor_r = right.density * (signal_speed_r - right.velocity.x);
        double signal_speed = (right.pressure - left.pressure + (factor_l * left.velocity.x) -
                               (factor_r * right.velocity.x)) /
                              (factor_l - factor_r);
        if (signal_speed_l <= 0 && 0 <= signal_speed) {
            contact_flux(flux, left, signal_speed_l, factor_l, signal_speed);
        }
        else {
            contact_flux(flux, right, signal_speed_r, factor_r, signal_speed);
        }
    }
    local_to_global(flux, basis);
}

static void hlle(void *flux_, const void *left_, const void *right_, const double *property,
                 Matrix basis)
{
    EulerConserved *flux = flux_;
    Euler left = global_to_local(left_, property, basis);
    Euler right = global_to_local(right_, property, basis);
    double gamma = property[EULER_HEAT_CAPACITY_RATIO];

    double sqrt_density_l = sqrt(left.density);
    double sqrt_density_r = sqrt(right.density);
    double weight = 1 / (sqrt_density_l + sqrt_density_r);
    double velocity_x =
        weight * ((sqrt_density_l * left.velocity.x) + (sqrt_density_r * right.velocity.x));

    double speed_of_sound2_l = gamma * left.pressure / left.density;
    double speed_of_sound2_r = gamma * right.pressure / right.density;

    double speed_of_sound_l = sqrt(speed_of_sound2_l);
    double speed_of_sound_r = sqrt(speed_of_sound2_r);

    double speed_of_sound2 =
        weight * ((sqrt_density_l * speed_of_sound2_l) + (sqrt_density_r * speed_of_sound2_r));

    double velocity_jump =
        sq(weight) * sqrt_density_l * sqrt_density_r * sq(right.velocity.x - left.velocity.x) / 2;

    double speed_of_sound = sqrt(speed_of_sound2 + velocity_jump);

    double signal_speed_l = fmin(left.velocity.x - speed_of_sound_l, velocity_x - speed_of_sound);
    double signal_speed_r = fmax(right.velocity.x + speed_of_sound_r, velocity_x + speed_of_sound);

    if (signal_speed_l >= 0) {
        *flux = physical_flux(left);
    }
    else if (signal_speed_r <= 0) {
        *flux = physical_flux(right);
    }
    else {
        average_flux(flux, left, right, signal_speed_l, signal_speed_r);
    }
    local_to_global(flux, basis);
}

static void ausmd(void *flux_, const void *left_, const void *right_, const double *property,
                  Matrix basis)
{
    EulerConserved *flux = flux_;
    Euler left = global_to_local(left_, property, basis);
    Euler right = global_to_local(right_, property, basis);
    double gamma = property[EULER_HEAT_CAPACITY_RATIO];

    double speed_of_sound_l = sqrt(gamma * left.pressure / left.density);
    double speed_of_sound_r = sqrt(gamma * right.pressure / right.density);
    double speed_of_sound = fmax(speed_of_sound_l, speed_of_sound_r);

    double pressure_over_density_l = left.pressure / left.density;
    double pressure_over_density_r = right.pressure / right.density;
    double pressure_over_density_sum = pressure_over_density_l + pressure_over_density_r;

    double mach_l = left.velocity.x / speed_of_sound;
    double mach_r = right.velocity.x / speed_of_sound;

    double velocity_plus_l;
    double pressure_plus_l;
    if (fabs(left.velocity.x) <= speed_of_sound) {
        double alpha_l = 2 * pressure_over_density_l / pressure_over_density_sum;
        double factor = left.velocity.x + speed_of_sound;
        velocity_plus_l = (alpha_l * sq(factor) / (4 * speed_of_sound)) +
                          ((1 - alpha_l) * (left.velocity.x + fabs(left.velocity.x)) / 2);
        pressure_plus_l = (left.pressure * sq(factor) / (4 * sq(speed_of_sound))) * (2 - mach_l);
    }
    else {
        velocity_plus_l = (left.velocity.x + fabs(left.velocity.x)) / 2;
        pressure_plus_l =
            (left.pressure * (left.velocity.x + fabs(left.velocity.x)) / 2) / left.velocity.x;
    }

    double velocity_minus_r;
    double pressure_minus_r;
    if (fabs(right.velocity.x) <= speed_of_sound) {
        double alpha_r = 2 * pressure_over_density_r / pressure_over_density_sum;
        double factor = right.velocity.x - speed_of_sound;
        velocity_minus_r = (-alpha_r * sq(factor) / (4 * speed_of_sound)) +
                           ((1 - alpha_r) * (right.velocity.x - fabs(right.velocity.x)) / 2);
        pressure_minus_r = (right.pressure * sq(factor) / (4 * sq(speed_of_sound))) * (2 + mach_r);
    }
    else {
        velocity_minus_r = (right.velocity.x - fabs(right.velocity.x)) / 2;
        pressure_minus_r =
            (right.pressure * (right.velocity.x - fabs(right.velocity.x)) / 2) / right.velocity.x;
    }

    double mass_flux = (velocity_plus_l * left.density) + (velocity_minus_r * right.density);
    double abs_mass_flux = fabs(mass_flux);
    double pressure_flux = pressure_plus_l + pressure_minus_r;

    double enthalpy_l = (left.energy + left.pressure) / left.density;
    double enthalpy_r = (right.energy + right.pressure) / right.density;

    flux->density = mass_flux;
    flux->momentum.x = (((mass_flux * (left.velocity.x + right.velocity.x)) -
                         (abs_mass_flux * (right.velocity.x - left.velocity.x))) /
                        2) +
                       pressure_flux;
    flux->momentum.y = ((mass_flux * (left.velocity.y + right.velocity.y)) -
                        (abs_mass_flux * (right.velocity.y - left.velocity.y))) /
                       2;
    flux->momentum.z = ((mass_flux * (left.velocity.z + right.velocity.z)) -
                        (abs_mass_flux * (right.velocity.z - left.velocity.z))) /
                       2;
    flux->energy =
        ((mass_flux * (enthalpy_l + enthalpy_r)) - (abs_mass_flux * (enthalpy_r - enthalpy_l))) / 2;
    local_to_global(flux, basis);
}

static void ausmdv(void *flux_, const void *left_, const void *right_, const double *property,
                   Matrix basis)
{
    EulerConserved *flux = flux_;
    Euler left = global_to_local(left_, property, basis);
    Euler right = global_to_local(right_, property, basis);
    double gamma = property[EULER_HEAT_CAPACITY_RATIO];

    double speed_of_sound_l = sqrt(gamma * left.pressure / left.density);
    double speed_of_sound_r = sqrt(gamma * right.pressure / right.density);
    double speed_of_sound = fmax(speed_of_sound_l, speed_of_sound_r);

    double pressure_over_density_l = left.pressure / left.density;
    double pressure_over_density_r = right.pressure / right.density;
    double pressure_over_density_sum = pressure_over_density_l + pressure_over_density_r;

    double mach_l = left.velocity.x / speed_of_sound;
    double mach_r = right.velocity.x / speed_of_sound;

    double velocity_plus_l;
    double pressure_plus_l;
    if (fabs(left.velocity.x) <= speed_of_sound) {
        double alpha_l = 2 * pressure_over_density_l / pressure_over_density_sum;
        double factor = left.velocity.x + speed_of_sound;
        velocity_plus_l = (alpha_l * sq(factor) / (4 * speed_of_sound)) +
                          ((1 - alpha_l) * (left.velocity.x + fabs(left.velocity.x)) / 2);
        pressure_plus_l = (left.pressure * sq(factor) / (4 * sq(speed_of_sound))) * (2 - mach_l);
    }
    else {
        velocity_plus_l = (left.velocity.x + fabs(left.velocity.x)) / 2;
        pressure_plus_l =
            left.pressure * (left.velocity.x + fabs(left.velocity.x)) / (2 * left.velocity.x);
    }

    double velocity_minus_r;
    double pressure_minus_r;
    if (fabs(right.velocity.x) <= speed_of_sound) {
        double alpha_r = 2 * pressure_over_density_r / pressure_over_density_sum;
        double factor = right.velocity.x - speed_of_sound;
        velocity_minus_r = (-alpha_r * sq(factor) / (4 * speed_of_sound)) +
                           ((1 - alpha_r) * (right.velocity.x - fabs(right.velocity.x)) / 2);
        pressure_minus_r = (right.pressure * sq(factor) / (4 * sq(speed_of_sound))) * (2 + mach_r);
    }
    else {
        velocity_minus_r = (right.velocity.x - fabs(right.velocity.x)) / 2;
        pressure_minus_r =
            right.pressure * (right.velocity.x - fabs(right.velocity.x)) / (2 * right.velocity.x);
    }

    double mass_flux = (velocity_plus_l * left.density) + (velocity_minus_r * right.density);
    double abs_mass_flux = fabs(mass_flux);
    double pressure_flux = pressure_plus_l + pressure_minus_r;

    double normal_momentum_flux_ausmd = ((mass_flux * (left.velocity.x + right.velocity.x)) -
                                         (abs_mass_flux * (right.velocity.x - left.velocity.x))) /
                                        2;

    double normal_momentum_flux_ausmv =
        (velocity_plus_l * left.momentum.x) + (velocity_minus_r * right.momentum.x);

    static const double switch_parameter = 10;
    double pressure_min = fmin(left.pressure, right.pressure);

    double switching = fmin(
        1.0 / 2, fmin(1, (switch_parameter * fabs(right.pressure - left.pressure)) / pressure_min));

    double normal_momentum_flux = (((1 + switching) * normal_momentum_flux_ausmv) +
                                   ((1 - switching) * normal_momentum_flux_ausmd)) /
                                  2;

    double enthalpy_l = (left.energy + left.pressure) / left.density;
    double enthalpy_r = (right.energy + right.pressure) / right.density;

    flux->density = mass_flux;
    flux->momentum.x = normal_momentum_flux + pressure_flux;
    flux->momentum.y = ((mass_flux * (left.velocity.y + right.velocity.y)) -
                        (abs_mass_flux * (right.velocity.y - left.velocity.y))) /
                       2;
    flux->momentum.z = ((mass_flux * (left.velocity.z + right.velocity.z)) -
                        (abs_mass_flux * (right.velocity.z - left.velocity.z))) /
                       2;
    flux->energy =
        ((mass_flux * (enthalpy_l + enthalpy_r)) - (abs_mass_flux * (enthalpy_r - enthalpy_l))) / 2;
    local_to_global(flux, basis);
}

static void lxf(void *flux_, const void *left_, const void *right_, const double *property,
                Matrix basis)
{
    EulerConserved *flux = flux_;
    Euler left = global_to_local(left_, property, basis);
    Euler right = global_to_local(right_, property, basis);
    double gamma = property[EULER_HEAT_CAPACITY_RATIO];

    double speed_of_sound_l = sqrt(gamma * left.pressure / left.density);
    double speed_of_sound_r = sqrt(gamma * right.pressure / right.density);

    double signal_speed =
        fmax(fabs(left.velocity.x) + speed_of_sound_l, fabs(right.velocity.x) + speed_of_sound_r);

    EulerConserved jump = compute_jump(left, right);
    jump.density *= signal_speed;
    vector2_imul(&jump.momentum, signal_speed);
    jump.energy *= signal_speed;

    EulerConserved flux_l = physical_flux(left);
    EulerConserved flux_r = physical_flux(right);
    flux->density = (flux_l.density + flux_r.density - jump.density) / 2;
    flux->momentum =
        vector2_div(vector2_sub(vector2_add(flux_l.momentum, flux_r.momentum), jump.momentum), 2);
    flux->energy = (flux_l.energy + flux_r.energy - jump.energy) / 2;
    local_to_global(flux, basis);
}

static void central(void *flux_, const void *left_, const void *right_, const double *property,
                    Matrix basis)
{
    EulerConserved *flux = flux_;
    Euler left = global_to_local(left_, property, basis);
    Euler right = global_to_local(right_, property, basis);

    EulerConserved flux_l = physical_flux(left);
    EulerConserved flux_r = physical_flux(right);

    flux->density = (flux_l.density + flux_r.density) / 2;
    flux->momentum = vector2_div(vector2_add(flux_l.momentum, flux_r.momentum), 2);
    flux->energy = (flux_l.energy + flux_r.energy) / 2;
    local_to_global(flux, basis);
}

Convective *euler2_convective(const char *name)
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
    if (!strcmp(name, "central")) {
        return central;
    }
    teal2_error("invalid convective flux (%s)", name);
}
