#include "gmres.h"

#include <math.h>

#include "teal/memory.h"
#include "teal/sync.h"
#include "teal/utils.h"

static void matvec(Simulation *sim, double omega, double dt, long m);

void gmres_m(Simulation *sim, double norm_fk, double omega, double dt)
{
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    const long n_cons = sim->eqns->n_cons;
    const long nn = n_inner_cells * n_cons;
    const long dim_krylov = sim->dim_krylov;
    const double tol_krylov = sim->tol_krylov * norm_fk;
    const double *fk = (void *)sim->advance.buf[2];
    double *dx = (void *)sim->advance.buf[4];
    double *w = (void *)sim->advance.buf[5];
    double(*V)[nn] = (void *)sim->advance.buf[6];
    double(*H)[dim_krylov] = (void *)sim->advance.buf[7];
    double *c = (void *)sim->advance.buf[8];
    double *s = (void *)sim->advance.buf[9];
    double *g = (void *)sim->advance.buf[10];
    double *y = (void *)sim->advance.buf[11];

    for (long i = 0; i < nn; ++i) V[0][i] = -fk[i] / norm_fk;
    g[0] = norm_fk;

    long m;
    for (m = 0; m < dim_krylov; ++m) {
        matvec(sim, omega, dt, m);

        // Arnoldi's method
        for (long j = 0; j < m + 1; ++j) {
            H[j][m] = sync_dot(V[j], w, nn);
            for (long i = 0; i < nn; ++i) w[i] -= H[j][m] * V[j][i];
        }
        H[m + 1][m] = sync_norm(w, nn);

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
        if (fabs(g[m + 1]) < tol_krylov) {
            m += 1;
            break;
        }

        // compute next Krylov basis vector
        for (long i = 0; i < nn; ++i) V[m + 1][i] = w[i] / H[m + 1][m];
    }
    sim->iter_krylov += m;

    memory_setzero(dx, nn, sizeof(*dx));
    for (long j = m - 1; j >= 0; --j) {
        y[j] = g[j];
        for (long k = j + 1; k < m; ++k) y[j] -= H[j][k] * y[k];
        y[j] /= H[j][j];
        for (long i = 0; i < nn; ++i) dx[i] += y[j] * V[j][i];
    }
}

static void matvec(Simulation *sim, double omega, double dt, long m)
{
    const long n_inner_cells = sim->eqns->mesh->n_inner_cells;
    const long n_cons = sim->eqns->n_cons;
    const long n_vars = sim->eqns->n_vars;
    const double(*dudt)[n_cons] = (void *)sim->eqns->vars.dudt;
    const double(*xk)[n_cons] = (void *)sim->advance.buf[1];
    const double(*rk)[n_cons] = (void *)sim->advance.buf[3];
    const double(*V)[n_inner_cells][n_cons] = (void *)sim->advance.buf[6];
    alias(update, sim->eqns->update.advance);
    double(*u)[n_vars] = (void *)sim->eqns->vars.u;
    double(*w)[n_cons] = (void *)sim->advance.buf[5];

    const double eps = omega / sync_norm(*V[m], n_inner_cells * n_cons);
    for (long i = 0; i < n_inner_cells; ++i) {
        for (long v = 0; v < n_cons; ++v) u[i][v] = xk[i][v] + eps * V[m][i][v];
        update(sim->eqns, u[i]);
    }
    equations_derivative(sim->eqns, sim->time);
    for (long i = 0; i < n_inner_cells; ++i)
        for (long v = 0; v < n_cons; ++v) w[i][v] = V[m][i][v] - dt * (dudt[i][v] - rk[i][v]) / eps;
}
