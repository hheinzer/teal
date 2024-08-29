#include <stdlib.h>

#include "equations.h"
#include "sync.h"
#include "teal/memory.h"
#include "teal/utils.h"

static void initialize_derivative(Equations *eqns);

static void apply_boundary_conditions(Equations *eqns, double time);

static void integrate_conv_fluxes(Equations *eqns);

static void integrate_reconstructed_conv_visc_fluxes(Equations *eqns);

static void integrate_reconstructed_conv_fluxes(Equations *eqns);

static void integrate_reconstructed_visc_fluxes(Equations *eqns);

static void compute_derivative(Equations *eqns, double time);

static void reconstruct(const Equations *eqns, const Vector3d dx, const double *u,
                        const Vector3d *dudx, double *us);

static void average(const Equations *eqns, const Vector4d c, const double *ul, const double *ur,
                    const Vector3d *dudxl, const Vector3d *dudxr, double *um, Vector3d *dudxm);

void equations_derivative(Equations *eqns, double time)
{
    initialize_derivative(eqns);
    apply_boundary_conditions(eqns, time);
    switch (eqns->space_order) {
        case 1: integrate_conv_fluxes(eqns); break;
        case 2:
            equations_gradient(eqns);
            if (eqns->limiter.limiter) equations_limiter(eqns);
            if (eqns->conv.flux && eqns->visc.flux)
                integrate_reconstructed_conv_visc_fluxes(eqns);
            else if (eqns->conv.flux)
                integrate_reconstructed_conv_fluxes(eqns);
            else if (eqns->visc.flux)
                integrate_reconstructed_visc_fluxes(eqns);
            else
                abort();
            break;
        default: abort();
    }
    compute_derivative(eqns, time);
}

static void initialize_derivative(Equations *eqns)
{
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const long n_cons = eqns->n_cons;
    double(*dudt)[n_cons] = (void *)eqns->vars.dudt;
    memory_setzero(dudt, n_inner_cells, sizeof(*dudt));
}

static void apply_boundary_conditions(Equations *eqns, double time)
{
    const long n_entities = eqns->mesh->n_entities;
    const long n_vars = eqns->n_vars;
    const alias(x, eqns->mesh->cell.center);
    const alias(cell, eqns->mesh->face.cell);
    const alias(r, eqns->mesh->face.basis);
    const alias(j_face, eqns->mesh->entity.j_face);
    const alias(state, eqns->bc.state);
    const alias(apply, eqns->bc.apply);
    const alias(compute, eqns->bc.compute);
    alias(update, eqns->update.boundary);
    double(*u)[n_vars] = (void *)eqns->vars.u;

    for (long e = 0; e < n_entities; ++e) {
        if (apply[e]) {
            for (long j = j_face[e]; j < j_face[e + 1]; ++j) {
                const long cL = cell[j][L], cR = cell[j][R];
                apply[e](eqns, r[j], state[e], u[cL], u[cR]);
                update(eqns, u[cR]);
            }
        }
        else if (compute[e]) {
            for (long j = j_face[e]; j < j_face[e + 1]; ++j) {
                const long cL = cell[j][L], cR = cell[j][R];
                compute[e](u[cR], u[cL], x[cR], time);
                update(eqns, u[cR]);
            }
        }
    }
}

static void integrate_conv_fluxes(Equations *eqns)
{
    const long n_inner_faces = eqns->mesh->n_inner_faces;
    const long n_bound_faces = eqns->mesh->n_ghost_faces;
    const long n_faces = eqns->mesh->n_faces;
    const long n_cons = eqns->n_cons;
    const long n_vars = eqns->n_vars;
    const alias(cell, eqns->mesh->face.cell);
    const alias(area, eqns->mesh->face.area);
    const alias(r, eqns->mesh->face.basis);
    alias(conv, eqns->conv.flux);
    double(*u)[n_vars] = (void *)eqns->vars.u;
    double(*dudt)[n_cons] = (void *)eqns->vars.dudt;
    double f[n_cons];

    sync_begin(eqns, *u, n_vars);

    for (long i = 0; i < n_inner_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        conv(eqns, r[i], u[cL], u[cR], f);
        for (long v = 0; v < n_cons; ++v) {
            dudt[cL][v] += f[v] * area[i];
            dudt[cR][v] -= f[v] * area[i];
        }
    }
    for (long i = n_inner_faces; i < n_inner_faces + n_bound_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        conv(eqns, r[i], u[cL], u[cR], f);
        for (long v = 0; v < n_cons; ++v) {
            dudt[cL][v] += f[v] * area[i];
        }
    }

    sync_wait(eqns);

    for (long i = n_inner_faces + n_bound_faces; i < n_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        conv(eqns, r[i], u[cL], u[cR], f);
        for (long v = 0; v < n_cons; ++v) {
            dudt[cL][v] += f[v] * area[i];
        }
    }

    sync_end(eqns);
}

static void integrate_reconstructed_conv_visc_fluxes(Equations *eqns)
{
    const long n_inner_faces = eqns->mesh->n_inner_faces;
    const long n_bound_faces = eqns->mesh->n_ghost_faces;
    const long n_faces = eqns->mesh->n_faces;
    const long n_cons = eqns->n_cons;
    const long n_vars = eqns->n_vars;
    const alias(cell, eqns->mesh->face.cell);
    const alias(area, eqns->mesh->face.area);
    const alias(r, eqns->mesh->face.basis);
    const alias(dx, eqns->mesh->face.to_cell);
    const alias(tl, eqns->mesh->face.correction);
    const double(*u)[n_vars] = (void *)eqns->vars.u;
    alias(conv, eqns->conv.flux);
    alias(visc, eqns->visc.flux);
    Vector3d(*dudx)[n_vars] = (void *)eqns->vars.dudx;
    double(*dudt)[n_cons] = (void *)eqns->vars.dudt;
    double ul[n_vars], ur[n_vars], um[n_vars], fc[n_cons], fv[n_cons];
    Vector3d dudxm[n_vars];

    sync_begin(eqns, **dudx, n_vars * N_DIMS);

    for (long i = 0; i < n_inner_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        reconstruct(eqns, dx[i][L], u[cL], dudx[cL], ul);
        reconstruct(eqns, dx[i][R], u[cR], dudx[cR], ur);
        average(eqns, tl[i], u[cL], u[cR], dudx[cL], dudx[cR], um, dudxm);
        conv(eqns, r[i], ul, ur, fc);
        visc(eqns, r[i], um, dudxm, fv);
        for (long v = 0; v < n_cons; ++v) {
            dudt[cL][v] += (fc[v] - fv[v]) * area[i];
            dudt[cR][v] -= (fc[v] - fv[v]) * area[i];
        }
    }
    for (long i = n_inner_faces; i < n_inner_faces + n_bound_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        reconstruct(eqns, dx[i][L], u[cL], dudx[cL], ul);
        average(eqns, tl[i], u[cL], u[cR], dudx[cL], dudx[cR], um, dudxm);
        conv(eqns, r[i], ul, u[cR], fc);
        visc(eqns, r[i], um, dudxm, fv);
        for (long v = 0; v < n_cons; ++v) {
            dudt[cL][v] += (fc[v] - fv[v]) * area[i];
        }
    }

    sync_wait(eqns);

    for (long i = n_inner_faces + n_bound_faces; i < n_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        reconstruct(eqns, dx[i][L], u[cL], dudx[cL], ul);
        reconstruct(eqns, dx[i][R], u[cR], dudx[cR], ur);
        average(eqns, tl[i], u[cL], u[cR], dudx[cL], dudx[cR], um, dudxm);
        conv(eqns, r[i], ul, ur, fc);
        visc(eqns, r[i], um, dudxm, fv);
        for (long v = 0; v < n_cons; ++v) {
            dudt[cL][v] += (fc[v] - fv[v]) * area[i];
        }
    }

    sync_end(eqns);
}

static void integrate_reconstructed_conv_fluxes(Equations *eqns)
{
    const long n_inner_faces = eqns->mesh->n_inner_faces;
    const long n_bound_faces = eqns->mesh->n_ghost_faces;
    const long n_faces = eqns->mesh->n_faces;
    const long n_cons = eqns->n_cons;
    const long n_vars = eqns->n_vars;
    const alias(cell, eqns->mesh->face.cell);
    const alias(area, eqns->mesh->face.area);
    const alias(r, eqns->mesh->face.basis);
    const alias(dx, eqns->mesh->face.to_cell);
    const double(*u)[n_vars] = (void *)eqns->vars.u;
    alias(conv, eqns->conv.flux);
    Vector3d(*dudx)[n_vars] = (void *)eqns->vars.dudx;
    double(*dudt)[n_cons] = (void *)eqns->vars.dudt;
    double ul[n_vars], ur[n_vars], fc[n_cons];

    sync_begin(eqns, **dudx, n_vars * N_DIMS);

    for (long i = 0; i < n_inner_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        reconstruct(eqns, dx[i][L], u[cL], dudx[cL], ul);
        reconstruct(eqns, dx[i][R], u[cR], dudx[cR], ur);
        conv(eqns, r[i], ul, ur, fc);
        for (long v = 0; v < n_cons; ++v) {
            dudt[cL][v] += fc[v] * area[i];
            dudt[cR][v] -= fc[v] * area[i];
        }
    }
    for (long i = n_inner_faces; i < n_inner_faces + n_bound_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        reconstruct(eqns, dx[i][L], u[cL], dudx[cL], ul);
        conv(eqns, r[i], ul, u[cR], fc);
        for (long v = 0; v < n_cons; ++v) {
            dudt[cL][v] += fc[v] * area[i];
        }
    }

    sync_wait(eqns);

    for (long i = n_inner_faces + n_bound_faces; i < n_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        reconstruct(eqns, dx[i][L], u[cL], dudx[cL], ul);
        reconstruct(eqns, dx[i][R], u[cR], dudx[cR], ur);
        conv(eqns, r[i], ul, ur, fc);
        for (long v = 0; v < n_cons; ++v) {
            dudt[cL][v] += fc[v] * area[i];
        }
    }

    sync_end(eqns);
}

static void integrate_reconstructed_visc_fluxes(Equations *eqns)
{
    const long n_inner_faces = eqns->mesh->n_inner_faces;
    const long n_bound_faces = eqns->mesh->n_ghost_faces;
    const long n_faces = eqns->mesh->n_faces;
    const long n_cons = eqns->n_cons;
    const long n_vars = eqns->n_vars;
    const alias(cell, eqns->mesh->face.cell);
    const alias(area, eqns->mesh->face.area);
    const alias(r, eqns->mesh->face.basis);
    const alias(tl, eqns->mesh->face.correction);
    const double(*u)[n_vars] = (void *)eqns->vars.u;
    alias(visc, eqns->visc.flux);
    Vector3d(*dudx)[n_vars] = (void *)eqns->vars.dudx;
    double(*dudt)[n_cons] = (void *)eqns->vars.dudt;
    double um[n_vars], fv[n_cons];
    Vector3d dudxm[n_vars];

    sync_begin(eqns, **dudx, n_vars * N_DIMS);

    for (long i = 0; i < n_inner_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        average(eqns, tl[i], u[cL], u[cR], dudx[cL], dudx[cR], um, dudxm);
        visc(eqns, r[i], um, dudxm, fv);
        for (long v = 0; v < n_cons; ++v) {
            dudt[cL][v] += -fv[v] * area[i];
            dudt[cR][v] -= -fv[v] * area[i];
        }
    }
    for (long i = n_inner_faces; i < n_inner_faces + n_bound_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        average(eqns, tl[i], u[cL], u[cR], dudx[cL], dudx[cR], um, dudxm);
        visc(eqns, r[i], um, dudxm, fv);
        for (long v = 0; v < n_cons; ++v) {
            dudt[cL][v] += -fv[v] * area[i];
        }
    }

    sync_wait(eqns);

    for (long i = n_inner_faces + n_bound_faces; i < n_faces; ++i) {
        const long cL = cell[i][L], cR = cell[i][R];
        average(eqns, tl[i], u[cL], u[cR], dudx[cL], dudx[cR], um, dudxm);
        visc(eqns, r[i], um, dudxm, fv);
        for (long v = 0; v < n_cons; ++v) {
            dudt[cL][v] += -fv[v] * area[i];
        }
    }

    sync_end(eqns);
}

static void compute_derivative(Equations *eqns, double time)
{
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const long n_cons = eqns->n_cons;
    const long n_vars = eqns->n_vars;
    const alias(x, eqns->mesh->cell.center);
    const alias(cv, eqns->mesh->cell.volume);
    alias(compute, eqns->source.compute);
    alias(prepare, eqns->source.prepare);
    double(*dudt)[n_cons] = (void *)eqns->vars.dudt;

    if (!compute) {
        for (long i = 0; i < n_inner_cells; ++i)
            for (long v = 0; v < n_cons; ++v) dudt[i][v] = -dudt[i][v] / cv[i];
    }
    else {
        const double(*u)[n_vars] = (void *)eqns->vars.u;
        double q[n_cons];
        memory_setzero(q, n_cons, sizeof(*q));

        if (prepare) prepare(eqns);
        for (long i = 0; i < n_inner_cells; ++i) {
            compute(q, u[i], x[i], time);
            for (long v = 0; v < n_cons; ++v) dudt[i][v] = -dudt[i][v] / cv[i] + q[v];
        }
    }
}

static void reconstruct(const Equations *eqns, const Vector3d dx, const double *u,
                        const Vector3d *dudx, double *us)
{
    const long n_vars = eqns->n_vars;
    for (long v = 0; v < n_vars; ++v) {
        us[v] = u[v];
        for (long d = 0; d < N_DIMS; ++d) us[v] += dudx[v][d] * dx[d];
    }
}

static void average(const Equations *eqns, const Vector4d c, const double *ul, const double *ur,
                    const Vector3d *dudxl, const Vector3d *dudxr, double *um, Vector3d *dudxm)
{
    const long n_vars = eqns->n_vars;
    for (long v = 0; v < n_vars; ++v) {
        um[v] = 0.5 * (ul[v] + ur[v]);
        for (long d = 0; d < N_DIMS; ++d) dudxm[v][d] = 0.5 * (dudxl[v][d] + dudxr[v][d]);
        double corr = -(ur[v] - ul[v]) / c[N_DIMS];
        for (long d = 0; d < N_DIMS; ++d) corr += dudxm[v][d] * c[d];
        for (long d = 0; d < N_DIMS; ++d) dudxm[v][d] -= corr * c[d];
    }
}
