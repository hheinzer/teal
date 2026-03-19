#include <math.h>

#include "euler2.h"
#include "teal2.h"

static const double pvrs_threshold = 2;
static const double prel_tolerance = 1e-6;
static const int max_newton = 20;

static double guess_pressure(EulerPrimitive left, EulerPrimitive right, double left_c,
                             double right_c, const double *gam)
{
    double cup = (left.density + right.density) * (left_c + right_c) / 4;
    double diff = left.velocity.x - right.velocity.x;
    double pvrs = fmax((left.pressure + right.pressure) + (diff * cup), prel_tolerance) / 2;
    double pmin = fmin(left.pressure, right.pressure);
    double pmax = fmax(left.pressure, right.pressure);
    double qmax = pmax / pmin;

    if (qmax <= pvrs_threshold && pmin <= pvrs && pvrs <= pmax) {
        return pvrs;  // PVRS: pressures are comparable and PVRS lies between them
    }

    if (pvrs < pmin) {  // two-rarefaction estimate (both sides rarefy)
        double ratio = pow(left.pressure / right.pressure, gam[1]);
        double guess = ((ratio * left.velocity.x / left_c) + (right.velocity.x / right_c) +
                        (gam[4] * (ratio - 1))) /
                       ((ratio / left_c) + (1 / right_c));
        double theta_l = 1 + (gam[7] * (left.velocity.x - guess) / left_c);
        double theta_r = 1 + (gam[7] * (guess - right.velocity.x) / right_c);
        return ((left.pressure * pow(theta_l, gam[3])) + (right.pressure * pow(theta_r, gam[3]))) /
               2;
    }

    // two-shock estimate (both sides shock)
    double scale_l = sqrt((gam[5] / left.density) / ((gam[6] * left.pressure) + pvrs));
    double scale_r = sqrt((gam[5] / right.density) / ((gam[6] * right.pressure) + pvrs));
    return ((scale_l * left.pressure) + (scale_r * right.pressure) + diff) / (scale_l + scale_r);
}

static void prefun(double *state_f, double *state_df, double pold, EulerPrimitive state,
                   double state_c, const double *gam)
{
    if (pold <= state.pressure) {  // rarefaction
        double ratio = pold / state.pressure;
        *state_f = gam[4] * state_c * (pow(ratio, gam[1]) - 1);
        *state_df = (1 / (state.density * state_c)) * pow(ratio, -gam[2]);
    }
    else {  // shock
        double coef_a = gam[5] / state.density;
        double coef_b = gam[6] * state.pressure;
        double root = sqrt(coef_a / (coef_b + pold));
        *state_f = (pold - state.pressure) * root;
        *state_df = (1 - ((pold - state.pressure) / (coef_b + pold) / 2)) * root;
    }
}

static void solve(double *star_u, double *star_p, EulerPrimitive left, EulerPrimitive right,
                  double left_c, double right_c, const double *gam)
{
    double diff = right.velocity.x - left.velocity.x;
    double pold = guess_pressure(left, right, left_c, right_c, gam);
    for (int iter = 0; iter < max_newton; iter++) {
        double left_f;
        double left_df;
        prefun(&left_f, &left_df, pold, left, left_c, gam);

        double right_f;
        double right_df;
        prefun(&right_f, &right_df, pold, right, right_c, gam);

        double pnew = pold - ((left_f + right_f + diff) / (left_df + right_df));
        if (2 * fabs((pnew - pold) / (pnew + pold)) <= prel_tolerance) {
            *star_u = (left.velocity.x + right.velocity.x + right_f - left_f) / 2;
            *star_p = pnew;
            return;
        }
        pold = fmax(prel_tolerance, pnew);
    }
    teal2_error("newton solver did not converge (%d)", max_newton);
}

static EulerPrimitive sample(EulerPrimitive left, EulerPrimitive right, double left_c,
                             double right_c, double star_u, double star_p, const double *gam,
                             double loc)
{
    if (loc <= star_u) {
        if (star_p <= left.pressure) {  // left rarefaction
            double loc_head = left.velocity.x - left_c;
            if (loc <= loc_head) {
                return left;
            }

            double star_c = left_c * pow(star_p / left.pressure, gam[1]);
            double loc_tail = star_u - star_c;
            if (loc > loc_tail) {
                return (EulerPrimitive){
                    .density = left.density * pow(star_p / left.pressure, 1 / gam[0]),
                    .velocity = {.x = star_u},
                    .pressure = star_p,
                };
            }

            double fan_c = gam[5] * (left_c + (gam[7] * (left.velocity.x - loc)));
            return (EulerPrimitive){
                .density = left.density * pow(fan_c / left_c, gam[4]),
                .velocity = {.x = gam[5] * (left_c + (gam[7] * left.velocity.x) + loc)},
                .pressure = left.pressure * pow(fan_c / left_c, gam[3]),
            };
        }

        // left shock
        double ratio = star_p / left.pressure;
        double loc_shock = left.velocity.x - (left_c * sqrt((gam[2] * ratio) + gam[1]));
        if (loc <= loc_shock) {
            return left;
        }
        return (EulerPrimitive){
            .density = left.density * (ratio + gam[6]) / ((ratio * gam[6]) + 1),
            .velocity = {.x = star_u},
            .pressure = star_p,
        };
    }

    if (star_p > right.pressure) {  // right shock
        double ratio = star_p / right.pressure;
        double loc_shock = right.velocity.x + (right_c * sqrt((gam[2] * ratio) + gam[1]));
        if (loc >= loc_shock) {
            return right;
        }
        return (EulerPrimitive){
            .density = right.density * (ratio + gam[6]) / ((ratio * gam[6]) + 1),
            .velocity = {.x = star_u},
            .pressure = star_p,
        };
    }

    // right rarefaction
    double loc_head = right.velocity.x + right_c;
    if (loc >= loc_head) {
        return right;
    }

    double star_c = right_c * pow(star_p / right.pressure, gam[1]);
    double loc_tail = star_u + star_c;
    if (loc <= loc_tail) {
        return (EulerPrimitive){
            .density = right.density * pow(star_p / right.pressure, 1 / gam[0]),
            .velocity = {.x = star_u},
            .pressure = star_p,
        };
    }

    double fan_c = gam[5] * (right_c - (gam[7] * (right.velocity.x - loc)));
    return (EulerPrimitive){
        .density = right.density * pow(fan_c / right_c, gam[4]),
        .velocity = {.x = gam[5] * (-right_c + (gam[7] * right.velocity.x) + loc)},
        .pressure = right.pressure * pow(fan_c / right_c, gam[3]),
    };
}

EulerPrimitive euler2_riemann(EulerPrimitive left, EulerPrimitive right, double gamma, double loc)
{
    double gam[9] = {
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

    double left_c = sqrt(gamma * left.pressure / left.density);
    double right_c = sqrt(gamma * right.pressure / right.density);
    if (gam[4] * (left_c + right_c) <= (right.velocity.x - left.velocity.x)) {
        teal2_error("vacuum state detected in riemann solver");
    }

    double star_u;
    double star_p;
    solve(&star_u, &star_p, left, right, left_c, right_c, gam);

    EulerPrimitive solution = sample(left, right, left_c, right_c, star_u, star_p, gam, loc);
    if (loc <= star_u) {
        solution.velocity.y = left.velocity.y;
        solution.velocity.z = left.velocity.z;
    }
    else {
        solution.velocity.y = right.velocity.y;
        solution.velocity.z = right.velocity.z;
    }

    return solution;
}
