#include "core/memory.h"
#include "core/sync.h"
#include "core/utils.h"
#include "equations.h"

static void initialize_derivative(Equations *eqns);

static void apply_boundary_conditions(Equations *eqns, double time);

static void integrate_conv_fluxes(Equations *eqns);

static void compute_derivative(Equations *eqns, double time);

static void integrate_reconstructed_conv_visc_fluxes(Equations *eqns);

static void integrate_reconstructed_conv_fluxes(Equations *eqns);

static void integrate_reconstructed_visc_fluxes(Equations *eqns);

void equations_derivative(Equations *eqns, double time)
{
    // dUdt = (-int_dV (F^C - F^V) dS + int_V Q dV ) / V
    switch (eqns->space_order) {
        case 1:
            initialize_derivative(eqns);
            apply_boundary_conditions(eqns, time);
            integrate_conv_fluxes(eqns);
            compute_derivative(eqns, time);
            break;
        case 2:
            initialize_derivative(eqns);
            apply_boundary_conditions(eqns, time);
            equations_gradient(eqns);
            if (eqns->limiter.func) equations_limiter(eqns);
            if (eqns->flux.conv && eqns->flux.visc)
                integrate_reconstructed_conv_visc_fluxes(eqns);
            else if (eqns->flux.conv)
                integrate_reconstructed_conv_fluxes(eqns);
            else
                integrate_reconstructed_visc_fluxes(eqns);
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
                const long cL = cell[j][L], cR = cell[j][R];
                apply[e](eqns, n[j], state[e], u[cL], u[cR]);
                update(eqns, u[cR]);
            }
        }
        else if (custom[e]) {
            for (long j = j_face[e]; j < j_face[e + 1]; ++j) {
                const long cL = cell[j][L], cR = cell[j][R];
                custom[e](u[cR], u[cL], x[cR], time);
                update(eqns, u[cR]);
            }
        }
    }
}

static void integrate_conv_fluxes(Equations *eqns)
{
    const long n_vars = eqns->n_vars;
    const long n_inner_faces = eqns->mesh->n_inner_faces;
    const long n_bound_faces = eqns->mesh->n_bound_faces;
    const long n_faces = eqns->mesh->n_faces;
    const ALIAS(n, eqns->mesh->face.normal);
    const ALIAS(cell, eqns->mesh->face.cell);
    const ALIAS(area, eqns->mesh->face.area);
    ALIAS(conv, eqns->flux.conv);
    double(*u)[n_vars] = (void *)eqns->vars.u;
    double(*dudt)[n_vars] = (void *)eqns->vars.dudt;
    double f[n_vars];
    memory_setzero(f, n_vars, sizeof(*f));

    sync_begin(eqns, *u, n_vars);

    for (long i = 0; i < n_inner_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        conv(eqns, n[i], u[cL], u[cR], f);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cL][v] += f[v] * area[i];
            dudt[cR][v] -= f[v] * area[i];
        }
    }
    for (long i = n_inner_faces; i < n_inner_faces + n_bound_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        conv(eqns, n[i], u[cL], u[cR], f);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cL][v] += f[v] * area[i];
        }
    }

    sync_wait(eqns);

    for (long i = n_inner_faces + n_bound_faces; i < n_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        conv(eqns, n[i], u[cL], u[cR], f);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cL][v] += f[v] * area[i];
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

static void reconstruct(const Equations *eqns, const double *dx, const double *u,
                        const double (*dudx)[N_DIMS], double *us)
{
    const long n_vars = eqns->n_vars;
    for (long v = 0; v < n_vars; ++v) {
        us[v] = u[v];
        for (long d = 0; d < N_DIMS; ++d) us[v] += dudx[v][d] * dx[d];
    }
}

static void average(const Equations *eqns, const double *tl, const double *ul, const double *ur,
                    const double (*dudxl)[N_DIMS], const double (*dudxr)[N_DIMS], double *um,
                    double (*dudxm)[N_DIMS])
{
    const long n_vars = eqns->n_vars;
    for (long v = 0; v < n_vars; ++v) {
        um[v] = 0.5 * (ul[v] + ur[v]);
        for (long d = 0; d < N_DIMS; ++d) dudxm[v][d] = 0.5 * (dudxl[v][d] + dudxr[v][d]);
        double corr = -(ur[v] - ul[v]) / tl[N_DIMS];
        for (long d = 0; d < N_DIMS; ++d) corr += dudxm[v][d] * tl[d];
        for (long d = 0; d < N_DIMS; ++d) dudxm[v][d] -= corr * tl[d];
    }
}

static void integrate_reconstructed_conv_visc_fluxes(Equations *eqns)
{
    const long n_vars = eqns->n_vars;
    const long n_inner_faces = eqns->mesh->n_inner_faces;
    const long n_bound_faces = eqns->mesh->n_bound_faces;
    const long n_faces = eqns->mesh->n_faces;
    const ALIAS(n, eqns->mesh->face.normal);
    const ALIAS(cell, eqns->mesh->face.cell);
    const ALIAS(dx, eqns->mesh->face.reconstruction);
    const ALIAS(tl, eqns->mesh->face.gradient_correction);
    const ALIAS(area, eqns->mesh->face.area);
    const double(*u)[n_vars] = (void *)eqns->vars.u;
    ALIAS(conv, eqns->flux.conv);
    ALIAS(visc, eqns->flux.visc);
    double(*dudx)[n_vars][N_DIMS] = (void *)eqns->vars.dudx;
    double(*dudt)[n_vars] = (void *)eqns->vars.dudt;
    double ul[n_vars], ur[n_vars], fc[n_vars], um[n_vars], dudxm[n_vars][N_DIMS], fv[n_vars];
    memory_setzero(fc, n_vars, sizeof(*fc));
    memory_setzero(fv, n_vars, sizeof(*fv));

    sync_begin(eqns, **dudx, n_vars * N_DIMS);

    for (long i = 0; i < n_inner_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        reconstruct(eqns, dx[i][L], u[cL], dudx[cL], ul);
        reconstruct(eqns, dx[i][R], u[cR], dudx[cR], ur);
        average(eqns, tl[i], u[cL], u[cR], dudx[cL], dudx[cR], um, dudxm);
        conv(eqns, n[i], ul, ur, fc);
        visc(eqns, n[i], um, dudxm, fv);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cL][v] += (fc[v] - fv[v]) * area[i];
            dudt[cR][v] -= (fc[v] - fv[v]) * area[i];
        }
    }
    for (long i = n_inner_faces; i < n_inner_faces + n_bound_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        reconstruct(eqns, dx[i][L], u[cL], dudx[cL], ul);
        average(eqns, tl[i], u[cL], u[cR], dudx[cL], dudx[cR], um, dudxm);
        conv(eqns, n[i], ul, u[cR], fc);
        visc(eqns, n[i], um, dudxm, fv);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cL][v] += (fc[v] - fv[v]) * area[i];
        }
    }

    sync_wait(eqns);

    for (long i = n_inner_faces + n_bound_faces; i < n_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        reconstruct(eqns, dx[i][L], u[cL], dudx[cL], ul);
        reconstruct(eqns, dx[i][R], u[cR], dudx[cR], ur);
        average(eqns, tl[i], u[cL], u[cR], dudx[cL], dudx[cR], um, dudxm);
        conv(eqns, n[i], ul, ur, fc);
        visc(eqns, n[i], um, dudxm, fv);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cL][v] += (fc[v] - fv[v]) * area[i];
        }
    }

    sync_end(eqns);
}

static void integrate_reconstructed_conv_fluxes(Equations *eqns)
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
    ALIAS(conv, eqns->flux.conv);
    double(*dudx)[n_vars][N_DIMS] = (void *)eqns->vars.dudx;
    double(*dudt)[n_vars] = (void *)eqns->vars.dudt;
    double ul[n_vars], ur[n_vars], fc[n_vars];
    memory_setzero(fc, n_vars, sizeof(*fc));

    sync_begin(eqns, **dudx, n_vars * N_DIMS);

    for (long i = 0; i < n_inner_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        reconstruct(eqns, dx[i][L], u[cL], dudx[cL], ul);
        reconstruct(eqns, dx[i][R], u[cR], dudx[cR], ur);
        conv(eqns, n[i], ul, ur, fc);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cL][v] += fc[v] * area[i];
            dudt[cR][v] -= fc[v] * area[i];
        }
    }
    for (long i = n_inner_faces; i < n_inner_faces + n_bound_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        reconstruct(eqns, dx[i][L], u[cL], dudx[cL], ul);
        conv(eqns, n[i], ul, u[cR], fc);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cL][v] += fc[v] * area[i];
        }
    }

    sync_wait(eqns);

    for (long i = n_inner_faces + n_bound_faces; i < n_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        reconstruct(eqns, dx[i][L], u[cL], dudx[cL], ul);
        reconstruct(eqns, dx[i][R], u[cR], dudx[cR], ur);
        conv(eqns, n[i], ul, ur, fc);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cL][v] += fc[v] * area[i];
        }
    }

    sync_end(eqns);
}

static void integrate_reconstructed_visc_fluxes(Equations *eqns)
{
    const long n_vars = eqns->n_vars;
    const long n_inner_faces = eqns->mesh->n_inner_faces;
    const long n_bound_faces = eqns->mesh->n_bound_faces;
    const long n_faces = eqns->mesh->n_faces;
    const ALIAS(n, eqns->mesh->face.normal);
    const ALIAS(cell, eqns->mesh->face.cell);
    const ALIAS(tl, eqns->mesh->face.gradient_correction);
    const ALIAS(area, eqns->mesh->face.area);
    const double(*u)[n_vars] = (void *)eqns->vars.u;
    ALIAS(visc, eqns->flux.visc);
    double(*dudx)[n_vars][N_DIMS] = (void *)eqns->vars.dudx;
    double(*dudt)[n_vars] = (void *)eqns->vars.dudt;
    double um[n_vars], dudxm[n_vars][N_DIMS], fv[n_vars];
    memory_setzero(fv, n_vars, sizeof(*fv));

    sync_begin(eqns, **dudx, n_vars * N_DIMS);

    for (long i = 0; i < n_inner_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        average(eqns, tl[i], u[cL], u[cR], dudx[cL], dudx[cR], um, dudxm);
        visc(eqns, n[i], um, dudxm, fv);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cL][v] += -fv[v] * area[i];
            dudt[cR][v] -= -fv[v] * area[i];
        }
    }
    for (long i = n_inner_faces; i < n_inner_faces + n_bound_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        average(eqns, tl[i], u[cL], u[cR], dudx[cL], dudx[cR], um, dudxm);
        visc(eqns, n[i], um, dudxm, fv);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cL][v] += -fv[v] * area[i];
        }
    }

    sync_wait(eqns);

    for (long i = n_inner_faces + n_bound_faces; i < n_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        average(eqns, tl[i], u[cL], u[cR], dudx[cL], dudx[cR], um, dudxm);
        visc(eqns, n[i], um, dudxm, fv);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cL][v] += -fv[v] * area[i];
        }
    }

    sync_end(eqns);
}
