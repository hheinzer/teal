#include "advance.h"

#include <math.h>

#include "core/memory.h"
#include "core/sync.h"
#include "core/utils.h"
#include "equations.h"

double explicit_euler(Simulation *sim, const double max_dt)
{
    const long n_cons = sim->eqns->n_cons;
    const long n_vars = sim->eqns->n_vars;
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    const double(*dudt)[n_cons] = (void *)sim->eqns->vars.dudt;
    ALIAS(update, sim->eqns->advance);
    double(*u)[n_vars] = (void *)sim->eqns->vars.u;

    const double dt0 = sim->cfl * equations_timestep(sim->eqns);
    const double dt = min(max_dt, dt0);

    equations_derivative(sim->eqns, sim->time);
    for (long i = 0; i < n_inner_cells; ++i) {
        for (long v = 0; v < n_cons; ++v) u[i][v] += dudt[i][v] * dt;
        update(sim->eqns, u[i]);
    }

    sim->time += dt;

    return dt0;
}

double lserk(Simulation *sim, const double max_dt)
{
    // van Leer, 1989
    static const double ALPHA[3][6][7] = {
        {
            {0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.3333, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1481, 0.4000, 1.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.0833, 0.2069, 0.4265, 1.0000, 0.0000, 0.0000},
            {0.0000, 0.0533, 0.1263, 0.2375, 0.4414, 1.0000, 0.0000},
            {0.0000, 0.0370, 0.0851, 0.1521, 0.2562, 0.4512, 1.0000},
        },
        {
            {0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},  // not second order
            {0.0000, 0.4242, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1918, 0.4929, 1.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1084, 0.2602, 0.5052, 1.0000, 0.0000, 0.0000},
            {0.0000, 0.0695, 0.1602, 0.2898, 0.5060, 1.0000, 0.0000},
            {0.0000, 0.0482, 0.1085, 0.1885, 0.3050, 0.5063, 1.0000},
        },
        {
            {0.0000, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000},  // not third order
            {0.0000, 0.6612, 1.0000, 0.0000, 0.0000, 0.0000, 0.0000},  // maybe third order
            {0.0000, 0.2884, 0.5010, 1.0000, 0.0000, 0.0000, 0.0000},
            {0.0000, 0.1666, 0.3027, 0.5275, 1.0000, 0.0000, 0.0000},
            {0.0000, 0.1067, 0.1979, 0.3232, 0.5201, 1.0000, 0.0000},
            {0.0000, 0.0742, 0.1393, 0.2198, 0.3302, 0.5181, 1.0000},
        },
    };
    const long time_order = sim->time_order;
    const long n_stages = sim->n_stages;
    const double *alpha = ALPHA[time_order - 1][n_stages - 1];

    const long n_cons = sim->eqns->n_cons;
    const long n_vars = sim->eqns->n_vars;
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    const double(*dudt)[n_cons] = (void *)sim->eqns->vars.dudt;
    ALIAS(update, sim->eqns->advance);
    double(*u0)[n_cons] = (void *)sim->advance.u0;
    double(*u)[n_vars] = (void *)sim->eqns->vars.u;

    const double dt0 = sim->cfl * equations_timestep(sim->eqns);
    const double dt = min(max_dt, dt0);

    for (long i = 0; i < n_inner_cells; ++i)
        for (long v = 0; v < n_cons; ++v) u0[i][v] = u[i][v];

    for (long k = 0; k < n_stages; ++k) {
        equations_derivative(sim->eqns, sim->time + alpha[k] * dt);
        for (long i = 0; i < n_inner_cells; ++i) {
            for (long v = 0; v < n_cons; ++v) u[i][v] = u0[i][v] + alpha[k + 1] * dudt[i][v] * dt;
            update(sim->eqns, u[i]);
        }
    }

    sim->time += dt;

    return dt0;
}

double implicit_euler(Simulation *sim, const double max_dt)
{
    const long n_newton = sim->n_newton;
    const long n_krylov = sim->n_krylov;
    const double tol_newton = sim->tol_newton;
    const double tol_krylov = sim->tol_krylov;

    const long n_cons = sim->eqns->n_cons;
    const long n_vars = sim->eqns->n_vars;
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    const double(*dudt)[n_cons] = (void *)sim->eqns->vars.dudt;
    ALIAS(update, sim->eqns->advance);
    double(*u0)[n_cons] = (void *)sim->advance.u0;
    double(*xk)[n_cons] = (void *)sim->advance.xk;
    double(*f0)[n_cons] = (void *)sim->advance.f0;
    double(*fk)[n_cons] = (void *)sim->advance.fk;
    double(*rk)[n_cons] = (void *)sim->advance.rk;
    double(*dx)[n_cons] = (void *)sim->advance.dx;
    double(*V)[n_inner_cells][n_cons] = (void *)sim->advance.V;
    double *g = (void *)sim->advance.g;
    double(*w)[n_cons] = (void *)sim->advance.w;
    double(*H)[n_krylov] = (void *)sim->advance.H;
    double *s = (void *)sim->advance.s;
    double *c = (void *)sim->advance.c;
    double *y = (void *)sim->advance.y;
    double(*u)[n_vars] = (void *)sim->eqns->vars.u;

    const double dt0 = sim->cfl * equations_timestep(sim->eqns);
    const double dt = min(max_dt, dt0);

    sim->time += dt;

    // initialize Newton solver
    equations_derivative(sim->eqns, sim->time);
    for (long i = 0; i < n_inner_cells; ++i) {
        for (long v = 0; v < n_cons; ++v) {
            xk[i][v] = u0[i][v] = u[i][v];
            fk[i][v] = f0[i][v] = -dt * dudt[i][v];
            rk[i][v] = dudt[i][v];
        }
    }
    const double norm_f0 = sync_norm(*f0, n_inner_cells * n_cons);
    double norm_fk = norm_f0;

    // Newton solver
    for (long n = 0; n < n_newton && norm_fk > tol_newton * norm_f0; ++n) {
        // initialize GMRES(m)
        for (long i = 0; i < n_inner_cells; ++i)
            for (long v = 0; v < n_cons; ++v) V[0][i][v] = -fk[i][v] / norm_fk;
        g[0] = norm_fk;

        // GMRES(m) solver
        long m;
        for (m = 0; m < n_krylov; ++m) {
            // matrix-vector
            const double eps = EPS / sync_norm(*V[m], n_inner_cells * n_cons);
            for (long i = 0; i < n_inner_cells; ++i) {
                for (long v = 0; v < n_cons; ++v) u[i][v] = xk[i][v] + eps * V[m][i][v];
                update(sim->eqns, u[i]);
            }
            equations_derivative(sim->eqns, sim->time);
            for (long i = 0; i < n_inner_cells; ++i)
                for (long v = 0; v < n_cons; ++v)
                    w[i][v] = V[m][i][v] - dt * (dudt[i][v] - rk[i][v]) / eps;

            // Arnoldi's method
            for (long j = 0; j < m + 1; ++j) {
                H[j][m] = sync_dot(*V[j], *w, n_inner_cells * n_cons);
                for (long i = 0; i < n_inner_cells; ++i)
                    for (long v = 0; v < n_cons; ++v) w[i][v] -= H[j][m] * V[j][i][v];
            }
            H[m + 1][m] = sync_norm(*w, n_inner_cells * n_cons);

            // apply past Givens rotations
            for (long j = 0; j < m; ++j) {
                const double tmp = c[j] * H[j][m] + s[j] * H[j + 1][m];
                H[j + 1][m] = -s[j] * H[j][m] + c[j] * H[j + 1][m];
                H[j][m] = tmp;
            }

            // compute Givens rotation
            const double r = sqrt(sq(H[m][m]) + sq(H[m + 1][m]));
            c[m] = H[m][m] / r;
            s[m] = H[m + 1][m] / r;
            H[m][m] = r;

            // update residual vector
            g[m + 1] = -s[m] * g[m];
            g[m] = c[m] * g[m];

            // check for convergence
            if (fabs(g[m + 1]) < tol_krylov * norm_fk) {
                m += 1;
                break;
            }

            // compute next Krylov basis vector
            for (long i = 0; i < n_inner_cells; ++i)
                for (long v = 0; v < n_cons; ++v) V[m + 1][i][v] = w[i][v] / H[m + 1][m];
        }

        // compute GMRES(m) result
        memory_setzero(dx, n_inner_cells, sizeof(*dx));
        for (long j = m - 1; j >= 0; --j) {
            y[j] = g[j];
            for (long k = j + 1; k < m; ++k) y[j] -= H[j][k] * y[k];
            y[j] /= H[j][j];
            for (long i = 0; i < n_inner_cells; ++i)
                for (long v = 0; v < n_cons; ++v) dx[i][v] += y[j] * V[j][i][v];
        }

        // advance Newton solver
        for (long i = 0; i < n_inner_cells; ++i) {
            for (long v = 0; v < n_cons; ++v) u[i][v] = xk[i][v] += dx[i][v];
            update(sim->eqns, u[i]);
        }
        equations_derivative(sim->eqns, sim->time);
        for (long i = 0; i < n_inner_cells; ++i) {
            for (long v = 0; v < n_cons; ++v) {
                fk[i][v] = u[i][v] - u0[i][v] - dt * dudt[i][v];
                rk[i][v] = dudt[i][v];
            }
        }
        norm_fk = sync_norm(*fk, n_inner_cells * n_cons);
    }
    ensure(norm_fk <= tol_newton * norm_f0);

    return dt0;
}
