/*
 * Toro 1999, sec. 4.9
 */
#include <assert.h>
#include <math.h>

#include "euler.h"
#include "teal/utils.h"

static const scalar pvrs_threshold = 2;
static const scalar prel_tolerance = 1e-6;
static const int max_newton = 20;

/* Return a robust starting value for the Newton solve of the star pressure. */
static scalar guess_pressure(const Euler *left, const Euler *right, scalar left_c, scalar right_c,
                             const scalar *gam)
{
    scalar cup = (left->density + right->density) * (left_c + right_c) / 4;
    scalar diff = left->velocity.x - right->velocity.x;
    scalar pvrs = fmax((left->pressure + right->pressure) + (diff * cup), prel_tolerance) / 2;
    scalar pmin = fmin(left->pressure, right->pressure);
    scalar pmax = fmax(left->pressure, right->pressure);
    scalar qmax = pmax / pmin;

    if (qmax <= pvrs_threshold && pmin <= pvrs && pvrs <= pmax) {
        return pvrs;  // PVRS: pressures are comparable and PVRS lies between them
    }

    if (pvrs < pmin) {  // two-rarefaction estimate (both sides rarefy)
        scalar ratio = pow(left->pressure / right->pressure, gam[1]);
        scalar guess = (ratio * left->velocity.x / left_c + right->velocity.x / right_c +
                        gam[4] * (ratio - 1)) /
                       (ratio / left_c + 1 / right_c);
        scalar theta_l = 1 + (gam[7] * (left->velocity.x - guess) / left_c);
        scalar theta_r = 1 + (gam[7] * (guess - right->velocity.x) / right_c);
        return (left->pressure * pow(theta_l, gam[3]) + right->pressure * pow(theta_r, gam[3])) / 2;
    }

    // two-shock estimate (both sides shock)
    scalar scale_l = sqrt((gam[5] / left->density) / (gam[6] * left->pressure + pvrs));
    scalar scale_r = sqrt((gam[5] / right->density) / (gam[6] * right->pressure + pvrs));
    return (scale_l * left->pressure + scale_r * right->pressure + diff) / (scale_l + scale_r);
}

/* Compute wave curve f_k(p) and derivative f'_k(p) for a single state. */
static void prefun(scalar *state_f, scalar *state_df, scalar pold, const Euler *state,
                   scalar state_c, const scalar *gam)
{
    if (pold <= state->pressure) {  // rarefaction
        scalar ratio = pold / state->pressure;
        *state_f = gam[4] * state_c * (pow(ratio, gam[1]) - 1);
        *state_df = (1 / (state->density * state_c)) * pow(ratio, -gam[2]);
    }
    else {  // shock
        scalar coef_a = gam[5] / state->density;
        scalar coef_b = gam[6] * state->pressure;
        scalar root = sqrt(coef_a / (coef_b + pold));
        *state_f = (pold - state->pressure) * root;
        *state_df = (1 - (pold - state->pressure) / (coef_b + pold) / 2) * root;
    }
}

/* Newton iteration to solve for (p*, u*) satisfying f_L(p*) + f_R(p*) + (u_R - u_L) = 0. */
static void solve(scalar *star_u, scalar *star_p, const Euler *left, const Euler *right,
                  scalar left_c, scalar right_c, const scalar *gam)
{
    scalar diff = right->velocity.x - left->velocity.x;
    scalar pold = guess_pressure(left, right, left_c, right_c, gam);
    for (int iter = 0; iter < max_newton; iter++) {
        scalar left_f;
        scalar left_df;
        prefun(&left_f, &left_df, pold, left, left_c, gam);

        scalar right_f;
        scalar right_df;
        prefun(&right_f, &right_df, pold, right, right_c, gam);

        scalar pnew = pold - ((left_f + right_f + diff) / (left_df + right_df));
        if (2 * fabs((pnew - pold) / (pnew + pold)) <= prel_tolerance) {
            *star_u = (left->velocity.x + right->velocity.x + right_f - left_f) / 2;
            *star_p = pnew;
            return;
        }
        pold = fmax(prel_tolerance, pnew);
    }
    error("newton solver did not converge after %d iterations", max_newton);
}

/* Return Euler state at the query location loc = x/t. */
static Euler sample(const Euler *left, const Euler *right, scalar left_c, scalar right_c,
                    scalar star_u, scalar star_p, const scalar *gam, scalar location)
{
    if (location <= star_u) {
        if (star_p <= left->pressure) {  // left rarefaction
            scalar loc_head = left->velocity.x - left_c;
            if (location <= loc_head) {
                return *left;
            }

            scalar star_c = left_c * pow(star_p / left->pressure, gam[1]);
            scalar loc_tail = star_u - star_c;
            if (location > loc_tail) {
                return (Euler){
                    .density = left->density * pow(star_p / left->pressure, 1 / gam[0]),
                    .velocity = {.x = star_u},
                    .pressure = star_p,
                };
            }

            scalar fan_c = gam[5] * (left_c + gam[7] * (left->velocity.x - location));
            return (Euler){
                .density = left->density * pow(fan_c / left_c, gam[4]),
                .velocity = {.x = gam[5] * (left_c + gam[7] * left->velocity.x + location)},
                .pressure = left->pressure * pow(fan_c / left_c, gam[3]),
            };
        }

        // left shock
        scalar ratio = star_p / left->pressure;
        scalar loc_shock = left->velocity.x - (left_c * sqrt((gam[2] * ratio) + gam[1]));
        if (location <= loc_shock) {
            return *left;
        }
        return (Euler){
            .density = left->density * (ratio + gam[6]) / (ratio * gam[6] + 1),
            .velocity = {.x = star_u},
            .pressure = star_p,
        };
    }

    if (star_p > right->pressure) {  // right shock
        scalar ratio = star_p / right->pressure;
        scalar loc_shock = right->velocity.x + (right_c * sqrt((gam[2] * ratio) + gam[1]));
        if (location >= loc_shock) {
            return *right;
        }
        return (Euler){
            .density = right->density * (ratio + gam[6]) / (ratio * gam[6] + 1),
            .velocity = {.x = star_u},
            .pressure = star_p,
        };
    }

    // right rarefaction
    scalar loc_head = right->velocity.x + right_c;
    if (location >= loc_head) {
        return *right;
    }

    scalar star_c = right_c * pow(star_p / right->pressure, gam[1]);
    scalar loc_tail = star_u + star_c;
    if (location <= loc_tail) {
        return (Euler){
            .density = right->density * pow(star_p / right->pressure, 1 / gam[0]),
            .velocity = {.x = star_u},
            .pressure = star_p,
        };
    }

    scalar fan_c = gam[5] * (right_c - gam[7] * (right->velocity.x - location));
    return (Euler){
        .density = right->density * pow(fan_c / right_c, gam[4]),
        .velocity = {.x = gam[5] * (-right_c + gam[7] * right->velocity.x + location)},
        .pressure = right->pressure * pow(fan_c / right_c, gam[3]),
    };
}

Euler euler_riemann(const Euler *left, const Euler *right, scalar gamma, scalar location)
{
    scalar gam[9] = {
        gamma,
        (gamma - 1) / (2 * gamma),
        (gamma + 1) / (2 * gamma),
        2 * gamma / (gamma - 1),
        2 / (gamma - 1),
        2 / (gamma + 1),
        (gamma - 1) / (gamma + 1),
        (gamma - 1) / 2,
        gamma - 1,
    };

    scalar left_c = sqrt(gamma * left->pressure / left->density);
    scalar right_c = sqrt(gamma * right->pressure / right->density);
    assert(gam[4] * (left_c + right_c) > (right->velocity.x - left->velocity.x));

    scalar star_u;
    scalar star_p;
    solve(&star_u, &star_p, left, right, left_c, right_c, gam);

    Euler solution = sample(left, right, left_c, right_c, star_u, star_p, gam, location);
    if (location <= star_u) {
        solution.velocity.y = left->velocity.y;
        solution.velocity.z = left->velocity.z;
    }
    else {
        solution.velocity.y = right->velocity.y;
        solution.velocity.z = right->velocity.z;
    }

    return solution;
}
