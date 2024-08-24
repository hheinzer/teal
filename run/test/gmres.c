#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "core/memory.h"
#include "core/utils.h"

#define N 1000

void init(double *A, double *b);

double residual(const double *A, const double *x, const double *b);

void gmres(double *x, const double *A, const double *b, long m, double tol);

int main(void)
{
    smart double *A = memory_calloc(N * N, sizeof(*A));
    smart double *x = memory_calloc(N, sizeof(*x));
    smart double *b = memory_calloc(N, sizeof(*b));
    init(A, b);
    printf("%g\n", residual(A, x, b));
    for (long i = 0; i < 100; ++i) gmres(x, A, b, 5, 1e-5);
    printf("%g (%g)\n", residual(A, x, b), 9.21086e-06);
}

void init(double *A, double *b)
{
    double(*a)[N] = (void *)A;
    for (long i = 0; i < N; ++i) {
        for (long j = 0; j < N; ++j) a[i][j] = 0.1 * rand() / (double)RAND_MAX;
        a[i][i] += 1;
        b[i] = rand() / (double)RAND_MAX;
    }
}

double residual(const double *A, const double *x, const double *b)
{
    const double(*a)[N] = (void *)A;
    double norm_r = 0, norm_b = 0;
    for (long i = 0; i < N; ++i) {
        double r = b[i];
        for (long j = 0; j < N; ++j) r -= a[i][j] * x[j];
        norm_r += sq(r);
        norm_b += sq(b[i]);
    }
    return sqrt(norm_r) / sqrt(norm_b);
}

void gmres(double *x, const double *A, const double *b, long m, double tol)
{
    const double(*a)[N] = (void *)A;
    smart double(*v)[N] = memory_calloc(m + 1, sizeof(*v));
    smart double *g = memory_calloc(m + 1, sizeof(*g));
    smart double *w = memory_calloc(N, sizeof(*w));
    smart double(*h)[m] = memory_calloc(m + 1, sizeof(*h));
    smart double *s = memory_calloc(m, sizeof(*s));
    smart double *c = memory_calloc(m, sizeof(*c));
    smart double *y = memory_calloc(m, sizeof(*y));

    double norm_b = 0, norm_v = 0;
    for (long i = 0; i < N; ++i) {
        v[0][i] = b[i];
        for (long j = 0; j < N; ++j) v[0][i] -= a[i][j] * x[j];
        norm_b += sq(b[i]);
        norm_v += sq(v[0][i]);
    }
    norm_b = sqrt(norm_b);
    norm_v = sqrt(norm_v);
    if (norm_v < tol * norm_b) return;

    for (long i = 0; i < N; ++i) v[0][i] /= norm_v;
    g[0] = norm_v;

    long k;
    for (k = 0; k < m; ++k) {
        // Arnoldi's method
        for (long i = 0; i < N; ++i) {
            w[i] = 0;
            for (long j = 0; j < N; ++j) w[i] += a[i][j] * v[k][j];
        }
        for (long i = 0; i < k + 1; ++i) {
            h[i][k] = 0;
            for (long j = 0; j < N; ++j) h[i][k] += v[i][j] * w[j];
            for (long j = 0; j < N; ++j) w[j] -= h[i][k] * v[i][j];
        }
        double norm_w = 0;
        for (long i = 0; i < N; ++i) norm_w += sq(w[i]);
        norm_w = sqrt(norm_w);
        h[k + 1][k] = norm_w;

        // Givens rotation
        for (long i = 0; i < k; ++i) {
            const double tmp = c[i] * h[i][k] + s[i] * h[i + 1][k];
            h[i + 1][k] = -s[i] * h[i][k] + c[i] * h[i + 1][k];
            h[i][k] = tmp;
        }
        const double t = sqrt(sq(h[k][k]) + sq(h[k + 1][k]));
        c[k] = h[k][k] / t;
        s[k] = h[k + 1][k] / t;
        h[k][k] = t;

        // update residual
        g[k + 1] = -s[k] * g[k];
        g[k] = c[k] * g[k];

        // check for convergence
        if (fabs(g[k + 1]) < tol * norm_b) {
            k += 1;
            break;
        };

        // compute next basis vector of Krylov subspace
        for (long i = 0; i < N; ++i) v[k + 1][i] = w[i] / h[k + 1][k];
    }

    // compute result
    for (long i = k - 1; i >= 0; --i) {
        y[i] = g[i];
        for (long j = i + 1; j < k; ++j) y[i] -= h[i][j] * y[j];
        y[i] /= h[i][i];
        for (long j = 0; j < N; ++j) x[j] += y[i] * v[i][j];
    }
}
