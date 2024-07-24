#include "core/memory.h"
#include "core/sync.h"
#include "core/utils.h"
#include "equations.h"

static void initialize_derivative(Equations *eqns);

static void apply_boundary_conditions(Equations *eqns, double time);

static void integrate_fluxes(Equations *eqns);

static void compute_derivative(Equations *eqns, double time);

static void integrate_reconstructed_fluxes(Equations *eqns);

void equations_derivative(Equations *eqns, double time)
{
    // dUdt = (-int_dV (F^C - F^V) dS + int_V Q dV ) / V
    switch (eqns->space_order) {
        case 1:
            initialize_derivative(eqns);
            apply_boundary_conditions(eqns, time);
            integrate_fluxes(eqns);
            compute_derivative(eqns, time);
            break;
        case 2:
            initialize_derivative(eqns);
            apply_boundary_conditions(eqns, time);
            equations_gradient(eqns);
            if (eqns->limiter.func) equations_limiter(eqns);
            integrate_reconstructed_fluxes(eqns);
            compute_derivative(eqns, time);
            break;
        default: error("unsupported space order '%ld'", eqns->space_order);
    }
}

static void initialize_derivative(Equations *eqns)
{
    const long n_vars = eqns->n_vars;
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    double(*dudt)[n_vars] = (void *)eqns->vars.dudt;
    memory_setzero(dudt, n_inner_cells, sizeof(*dudt));
}

static void apply_boundary_conditions(Equations *eqns, double time)
{
    const long n_vars = eqns->n_vars;
    const long n_entities = eqns->mesh->n_entities;
    const ALIAS(j_face, eqns->mesh->entity.j_face);
    const ALIAS(n, eqns->mesh->face.normal);
    const ALIAS(cell, eqns->mesh->face.cell);
    const ALIAS(x, eqns->mesh->cell.center);
    const ALIAS(apply, eqns->bc.apply);
    const ALIAS(state, eqns->bc.state);
    const ALIAS(custom, eqns->bc.custom);
    ALIAS(update, eqns->boundary);
    double(*u)[n_vars] = (void *)eqns->vars.u;

    for (long e = 0; e < n_entities; ++e) {
        if (apply[e]) {
            for (long j = j_face[e]; j < j_face[e + 1]; ++j) {
                apply[e](eqns, n[j], state[e], u[cell[j][L]], u[cell[j][R]]);
                update(eqns, u[cell[j][R]]);
            }
        }
        else if (custom[e]) {
            for (long j = j_face[e]; j < j_face[e + 1]; ++j) {
                custom[e](u[cell[j][R]], u[cell[j][L]], x[cell[j][R]], time);
                update(eqns, u[cell[j][R]]);
            }
        }
    }
}

static void integrate_fluxes(Equations *eqns)
{
    const long n_vars = eqns->n_vars;
    const long n_inner_faces = eqns->mesh->n_inner_faces;
    const long n_bound_faces = eqns->mesh->n_bound_faces;
    const long n_faces = eqns->mesh->n_faces;
    const ALIAS(n, eqns->mesh->face.normal);
    const ALIAS(cell, eqns->mesh->face.cell);
    const ALIAS(area, eqns->mesh->face.area);
    ALIAS(flux, eqns->flux.func);
    double(*u)[n_vars] = (void *)eqns->vars.u;
    double(*dudt)[n_vars] = (void *)eqns->vars.dudt;
    double f[n_vars];
    memory_setzero(f, n_vars, sizeof(*f));

    sync_begin(eqns, *u, n_vars);

    for (long i = 0; i < n_inner_faces; ++i) {
        flux(eqns, n[i], u[cell[i][L]], u[cell[i][R]], f);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cell[i][L]][v] += f[v] * area[i];
            dudt[cell[i][R]][v] -= f[v] * area[i];
        }
    }
    for (long i = n_inner_faces; i < n_inner_faces + n_bound_faces; ++i) {
        flux(eqns, n[i], u[cell[i][L]], u[cell[i][R]], f);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cell[i][L]][v] += f[v] * area[i];
        }
    }

    sync_wait(eqns);

    for (long i = n_inner_faces + n_bound_faces; i < n_faces; ++i) {
        flux(eqns, n[i], u[cell[i][L]], u[cell[i][R]], f);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cell[i][L]][v] += f[v] * area[i];
        }
    }

    sync_end(eqns);
}

static void compute_derivative(Equations *eqns, double time)
{
    const long n_vars = eqns->n_vars;
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const ALIAS(x, eqns->mesh->cell.center);
    const ALIAS(volume, eqns->mesh->cell.volume);
    ALIAS(source, eqns->source);
    double(*dudt)[n_vars] = (void *)eqns->vars.dudt;

    if (!source) {
        for (long i = 0; i < n_inner_cells; ++i)
            for (long v = 0; v < n_vars; ++v) dudt[i][v] = -dudt[i][v] / volume[i];
    }
    else {
        const double(*u)[n_vars] = (void *)eqns->vars.u;
        double q[n_vars];
        memory_setzero(q, n_vars, sizeof(*q));

        for (long i = 0; i < n_inner_cells; ++i) {
            source(q, u[i], x[i], time);
            for (long v = 0; v < n_vars; ++v) dudt[i][v] = -dudt[i][v] / volume[i] + q[v];
        }
    }
}

static void integrate_reconstructed_fluxes(Equations *eqns)
{
    const long n_vars = eqns->n_vars;
    const long n_inner_faces = eqns->mesh->n_inner_faces;
    const long n_bound_faces = eqns->mesh->n_bound_faces;
    const long n_faces = eqns->mesh->n_faces;
    const ALIAS(n, eqns->mesh->face.normal);
    const ALIAS(cell, eqns->mesh->face.cell);
    const ALIAS(dx, eqns->mesh->face.reconstruction);
    const ALIAS(area, eqns->mesh->face.area);
    const double(*u)[n_vars] = (void *)eqns->vars.u;
    ALIAS(flux, eqns->flux.func);
    double(*dudx)[n_vars][N_DIMS] = (void *)eqns->vars.dudx;
    double(*dudt)[n_vars] = (void *)eqns->vars.dudt;
    double ul[n_vars], ur[n_vars], f[n_vars];
    memory_setzero(f, n_vars, sizeof(*f));

    sync_begin(eqns, **dudx, n_vars * N_DIMS);

    for (long i = 0; i < n_inner_faces; ++i) {
        for (long v = 0; v < n_vars; ++v) {
            ul[v] = u[cell[i][L]][v];
            ur[v] = u[cell[i][R]][v];
            for (long d = 0; d < N_DIMS; ++d) {
                ul[v] += dudx[cell[i][L]][v][d] * dx[i][L][d];
                ur[v] += dudx[cell[i][R]][v][d] * dx[i][R][d];
            }
        }
        flux(eqns, n[i], ul, ur, f);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cell[i][L]][v] += f[v] * area[i];
            dudt[cell[i][R]][v] -= f[v] * area[i];
        }
    }

    for (long i = n_inner_faces; i < n_inner_faces + n_bound_faces; ++i) {
        for (long v = 0; v < n_vars; ++v) {
            ul[v] = u[cell[i][L]][v];
            for (long d = 0; d < N_DIMS; ++d) {
                ul[v] += dudx[cell[i][L]][v][d] * dx[i][L][d];
            }
        }
        flux(eqns, n[i], ul, u[cell[i][R]], f);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cell[i][L]][v] += f[v] * area[i];
        }
    }

    sync_wait(eqns);

    for (long i = n_inner_faces + n_bound_faces; i < n_faces; ++i) {
        for (long v = 0; v < n_vars; ++v) {
            ul[v] = u[cell[i][L]][v];
            ur[v] = u[cell[i][R]][v];
            for (long d = 0; d < N_DIMS; ++d) {
                ul[v] += dudx[cell[i][L]][v][d] * dx[i][L][d];
                ur[v] += dudx[cell[i][R]][v][d] * dx[i][R][d];
            }
        }
        flux(eqns, n[i], ul, ur, f);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cell[i][L]][v] += f[v] * area[i];
        }
    }

    sync_end(eqns);
}
