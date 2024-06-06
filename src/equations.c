#include "equations.h"

#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "array.h"
#include "global.h"
#include "hdf5_io.h"
#include "limiter.h"
#include "memory.h"
#include "mesh.h"
#include "sync.h"
#include "utils.h"

static char m_limiter[128];  // limiter function

static void output_create(Fields *fields);
static void initialize_derivative(Equations *eqns);
static void apply_boundary_conditions(Equations *eqns, const double time);
static void integrate_sources(Equations *eqns, const double time);
static void integrate_local_fluxes(Equations *eqns);
static void integrate_non_local_fluxes(Equations *eqns);
static void compute_local_gradients(Equations *eqns);
static void compute_non_local_gradients(Equations *eqns);
static void limit_gradients(Equations *eqns);
static void integrate_local_reconstructed_fluxes(Equations *eqns);
static void integrate_non_local_reconstructed_fluxes(Equations *eqns);
static void scale_derivative(Equations *eqns);

Equations equations_create(const Mesh *mesh, const Fields *vars, const Fields *user)
{
    Equations eqns = {};

    const long n_vars = vars->n_fields;
    eqns.vars = *vars;
    eqns.vars.u = memory_calloc(mesh->n_cells * n_vars, sizeof(*eqns.vars.u));
    eqns.vars.dudt = memory_calloc(mesh->n_cells * n_vars, sizeof(*eqns.vars.dudt));
    output_create(&eqns.vars);

    if (user) {
        eqns.user = *user;
        output_create(&eqns.user);
    }

    eqns.mesh = mesh;

    const ALIAS(bc, eqns.mesh->entity.bc.name);
    for (long e = 0; e < eqns.mesh->n_entities; ++e) {
        if (bc[e] && !strcmp(bc[e], "periodic")) {
            free(bc[e]);
            bc[e] = 0;
        }
    }

    equations_set_space_order(&eqns, 2, "minmod", 0);

    const long n_send = eqns.mesh->sync.i_send[eqns.mesh->size];
    eqns.sync.buf_u = memory_calloc(n_send * n_vars, sizeof(*eqns.sync.buf_u));
    eqns.sync.recv_u = memory_calloc(mesh->size, sizeof(*eqns.sync.recv_u));
    eqns.sync.send_u = memory_calloc(mesh->size, sizeof(*eqns.sync.send_u));
    eqns.sync.buf_dudx = memory_calloc(n_send * n_vars * N_DIMS, sizeof(*eqns.sync.buf_dudx));
    eqns.sync.recv_dudx = memory_calloc(mesh->size, sizeof(*eqns.sync.recv_dudx));
    eqns.sync.send_dudx = memory_calloc(mesh->size, sizeof(*eqns.sync.send_dudx));

    return eqns;
}

void equations_free(Equations *eqns)
{
    free(eqns->vars.u);
    free(eqns->vars.dudx);
    free(eqns->vars.dudt);
    for (long i = 0; i < eqns->vars.n_fields; ++i) free(eqns->vars.name[i]);
    free(eqns->vars.name);
    free(eqns->vars.output.dim);
    if (eqns->vars.output.name)
        for (long i = 0; i < eqns->vars.output.n_dims; ++i) free(eqns->vars.output.name[i]);
    free(eqns->vars.output.name);
    eqns->vars = (Fields){};

    free(eqns->user.output.dim);
    if (eqns->user.output.name)
        for (long i = 0; i < eqns->user.output.n_dims; ++i) free(eqns->user.output.name[i]);
    free(eqns->user.output.name);
    eqns->vars = (Fields){};

    free(eqns->buf);

    free(eqns->sync.buf_u);
    free(eqns->sync.recv_u);
    free(eqns->sync.send_u);
    free(eqns->sync.buf_dudx);
    free(eqns->sync.recv_dudx);
    free(eqns->sync.send_dudx);

    *eqns = (Equations){};
}

void equations_set_initial_condition(Equations *eqns, Function *initial)
{
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const ALIAS(x, eqns->mesh->cell.center);
    FIELDS(u, eqns->vars);

    for (long i = 0; i < n_inner_cells; ++i) initial(x[i], 0, 0, u[i]);
}

void equations_set_initial_state(Equations *eqns, const long nu, const long *u, const double *state)
{
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    FIELDS(vars, eqns->vars);

    for (long i = 0; i < n_inner_cells; ++i)
        for (long j = 0; j < nu; ++j) vars[i][u[j]] = state[j];
}

void equations_set_space_order(Equations *eqns, const long space_order, const char *limiter,
                               const double k)
{
    assert(1 <= space_order && space_order <= 2 && "invalid space order");
    eqns->space_order = space_order;
    if (space_order == 2) {
        const long n_cells = eqns->mesh->n_cells;
        const long n_vars = eqns->vars.n_fields;
        eqns->vars.dudx = memory_calloc(n_cells * n_vars * N_DIMS, sizeof(*eqns->vars.dudx));
    }
    else {
        free(eqns->vars.dudx);
        eqns->vars.dudx = 0;
    }

    if (!limiter) {
        eqns->limiter = 0;
        strncpy(m_limiter, "none", sizeof(m_limiter) - 1);
    }
    else if (!strcmp(limiter, "barth jespersen") || !strcmp(limiter, "minmod")) {
        eqns->limiter = limiter_barth_jespersen;
        strncpy(m_limiter, limiter, sizeof(m_limiter) - 1);

        const long n_inner_cells = eqns->mesh->n_inner_cells;
        eqns->buf = memory_calloc(n_inner_cells, sizeof(*eqns->buf));
    }
    else if (!strcmp(limiter, "venkatakrishnan") || !strcmp(limiter, "venk")) {
        eqns->limiter = limiter_venkatakrishnan;
        strncpy(m_limiter, limiter, sizeof(m_limiter) - 1);

        eqns->k = k;
        const long n_inner_cells = eqns->mesh->n_inner_cells;
        const ALIAS(v, eqns->mesh->cell.volume);
        eqns->buf = memory_calloc(n_inner_cells, sizeof(*eqns->buf));
        for (long i = 0; i < n_inner_cells; ++i) eqns->buf[i] = pow(k * pow(v[i], 1.0 / N_DIMS), 3);
    }
    else {
        assert("unsupported limiter function");
    }
}

void equations_print(const Equations *eqns, const char *name)
{
    const FIELDS(u, eqns->vars);
    const long n_vars = eqns->vars.n_fields;
    double vars_min[n_vars], vars_max[n_vars];
    sync_min(array_min_s(*u, eqns->mesh->n_inner_cells, n_vars, vars_min), n_vars);
    sync_max(array_max_s(*u, eqns->mesh->n_inner_cells, n_vars, vars_max), n_vars);

    if (eqns->mesh->rank == 0) {
        printf("%s equations summary:\n", name);

        printf(" | " FMT_KEY ": %ld\n", "space order", eqns->space_order);
        printf(" | " FMT_KEY ": %s (k = %g)\n", "limiter", m_limiter, eqns->k);

        for (long n = 0, i = 0; i < eqns->vars.output.n_dims; ++i) {
            char key[128];
            snprintf(key, sizeof(key), "min/max %s", eqns->vars.output.name[i]);
            printf(" | " FMT_KEY ": ", key);
            array_print(&vars_min[n], eqns->vars.output.dim[i], " / ");
            array_print(&vars_max[n], eqns->vars.output.dim[i], "\n");
            n += eqns->vars.output.dim[i];
        }
    }
}

void equations_write(const Equations *eqns, const char *prefix, const long count, const double time)
{
    char fname[128], lname[128];
    snprintf(fname, sizeof(fname), "%s_%05ld.vtkhdf", prefix, count);
    snprintf(lname, sizeof(lname), "%s_mesh.vtkhdf", utils_basename(prefix));

    // https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html#unstructured-grid
    hid_t file = hdf5_file_create(fname);
    hid_t vtkhdf = hdf5_group_create(file, "VTKHDF");

    const long version[] = {2, 0};
    hdf5_write_attribute(vtkhdf, "Version", version, 1, HDF5_DIMS(2));
    hdf5_write_attribute(vtkhdf, "Type", "UnstructuredGrid");

    hdf5_link_create(vtkhdf, "NumberOfPoints", lname, "/VTKHDF/NumberOfPoints");
    hdf5_link_create(vtkhdf, "NumberOfConnectivityIds", lname, "/VTKHDF/NumberOfConnectivityIds");
    hdf5_link_create(vtkhdf, "NumberOfCells", lname, "/VTKHDF/NumberOfCells");
    hdf5_link_create(vtkhdf, "Points", lname, "/VTKHDF/Points");
    hdf5_link_create(vtkhdf, "Connectivity", lname, "/VTKHDF/Connectivity");
    hdf5_link_create(vtkhdf, "Offsets", lname, "/VTKHDF/Offsets");
    hdf5_link_create(vtkhdf, "Types", lname, "/VTKHDF/Types");

    hid_t field_data = hdf5_group_create(vtkhdf, "FieldData");

    const long rank = eqns->mesh->rank;
    hdf5_write_dataset(field_data, "TimeValue", &time, 1, HDF5_DIMS(rank == 0));

    hdf5_group_close(field_data);

    hid_t cell_data = hdf5_group_create(vtkhdf, "CellData");

    const long n_user = eqns->user.n_fields;
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const long n_vdims = eqns->vars.output.n_dims;
    const long n_udims = eqns->user.output.n_dims;
    const ALIAS(vdim, eqns->vars.output.dim);
    const ALIAS(udim, eqns->user.output.dim);
    long n_buf = array_max(vdim, n_vdims);
    if (n_user) n_buf = MAX(n_buf, array_max(udim, n_udims));
    cleanup double *buf = memory_calloc(n_inner_cells * n_buf, sizeof(*buf));

    const ALIAS(vname, eqns->vars.output.name);
    const FIELDS(vars, eqns->vars);
    for (long n = 0, j = 0; j < n_vdims; ++j) {
        for (long i = 0; i < n_inner_cells; ++i)
            for (long d = 0; d < vdim[j]; ++d) buf[i * vdim[j] + d] = vars[i][n + d];
        hdf5_write_dataset(cell_data, vname[j], buf, 2, HDF5_DIMS(n_inner_cells, vdim[j]));
        n += vdim[j];
    }

    if (eqns->user.n_fields) {
        const ALIAS(x, eqns->mesh->cell.center);
        const ALIAS(uname, eqns->user.output.name);
        double user[n_user] = {};
        for (long n = 0, j = 0; j < n_udims; ++j) {
            for (long i = 0; i < n_inner_cells; ++i) {
                eqns->user.compute(x[i], time, vars[i], user);
                for (long d = 0; d < udim[j]; ++d) buf[i * udim[j] + d] = user[n + d];
            }
            hdf5_write_dataset(cell_data, uname[j], buf, 2, HDF5_DIMS(n_inner_cells, udim[j]));
            n += udim[j];
        }
    }

    hdf5_group_close(cell_data);
    hdf5_group_close(vtkhdf);
    hdf5_file_close(file);
}

void equations_time_derivative(Equations *eqns, const double time)
{
    // dU_dt = (int_V Q dV - int_A F(U) * n dA) / V
    sync_u_begin(eqns);
    initialize_derivative(eqns);
    apply_boundary_conditions(eqns, time);
    if (eqns->source) integrate_sources(eqns, time);
    if (eqns->space_order == 1) {
        integrate_local_fluxes(eqns);
        sync_u_wait(eqns);
        integrate_non_local_fluxes(eqns);
    }
    else {
        compute_local_gradients(eqns);
        sync_u_wait(eqns);
        compute_non_local_gradients(eqns);
        if (eqns->limiter) limit_gradients(eqns);
        sync_dudx_begin(eqns);
        integrate_local_reconstructed_fluxes(eqns);
        sync_dudx_wait(eqns);
        integrate_non_local_reconstructed_fluxes(eqns);
        sync_dudx_end(eqns);
    }
    scale_derivative(eqns);
    sync_u_end(eqns);
}

void equations_residual(const Equations *eqns, double *residual)
{
    const long n_vars = eqns->vars.n_fields;
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const double total_volume = eqns->mesh->volume;
    const ALIAS(volume, eqns->mesh->cell.volume);
    const DERIVS(dudt, eqns->vars);

    for (long v = 0; v < n_vars; ++v) residual[v] = 0;
    for (long i = 0; i < n_inner_cells; ++i)
        for (long v = 0; v < n_vars; ++v) residual[v] += volume[i] * dudt[i][v] * dudt[i][v];
    sync_sum(residual, n_vars);
    for (long v = 0; v < n_vars; ++v) residual[v] = sqrt(residual[v] / total_volume);
}

static void output_create(Fields *fields)
{
    long ndim = 0;
    long *dim = memory_calloc(fields->n_fields, sizeof(*dim));
    char **name = memory_calloc(fields->n_fields, sizeof(*name));
    for (long i = 0; i < fields->n_fields; ++i) {
        dim[ndim] += 1;
        if (!name[ndim]) {
            name[ndim] = utils_strdup(fields->name[i]);
            const char *dash = strrchr(fields->name[i], '-');
            const long pos = dash - fields->name[i] - strlen(fields->name[i]);
            if (dash && pos == -2) name[ndim][strlen(name[ndim]) - 2] = 0;
        }
        if (i == fields->n_fields - 1 ||
            strncmp(fields->name[i + 1], name[ndim], strlen(name[ndim])))
            ndim += 1;
    }
    fields->output.n_dims = ndim;
    fields->output.dim = memory_realloc(dim, ndim, sizeof(*dim));
    fields->output.name = memory_realloc(name, ndim, sizeof(*name));
}

static void initialize_derivative(Equations *eqns)
{
    const long n_vars = eqns->vars.n_fields;
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    DERIVS(dudt, eqns->vars);

    for (long i = 0; i < n_inner_cells; ++i)
        for (long v = 0; v < n_vars; ++v) dudt[i][v] = 0;
}

static void apply_boundary_conditions(Equations *eqns, const double time)
{
    const long n_entities = eqns->mesh->n_entities;
    const ALIAS(apply, eqns->mesh->entity.bc.apply);
    const ALIAS(j_face, eqns->mesh->entity.j_face);
    const ALIAS(context, eqns->mesh->entity.bc.context);
    const ALIAS(cell, eqns->mesh->face.cell);
    const ALIAS(n, eqns->mesh->face.normal);
    const ALIAS(x, eqns->mesh->cell.center);
    FIELDS(u, eqns->vars);

    for (long e = 0; e < n_entities; ++e) {
        if (!apply[e]) continue;  // no boundary condition
        if (!context[e].custom) {
            // regular boundary condition
            for (long j = j_face[e]; j < j_face[e + 1]; ++j)
                apply[e](context[e], n[j], u[cell[j][0]], u[cell[j][1]]);
        }
        else {
            // custom boundary condition
            context[e].state = &time;
            for (long j = j_face[e]; j < j_face[e + 1]; ++j)
                apply[e](context[e], x[cell[j][1]], u[cell[j][0]], u[cell[j][1]]);
        }
    }
}

static void integrate_sources(Equations *eqns, const double time)
{
    const long n_vars = eqns->vars.n_fields;
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const ALIAS(xc, eqns->mesh->cell.center);
    const ALIAS(wc, eqns->mesh->cell.gauss_weight);
    DERIVS(dudt, eqns->vars);
    double q[n_vars] = {};

    for (long i = 0; i < n_inner_cells; ++i) {
        eqns->source(xc[i], time, 0, q);
        for (long v = 0; v < n_vars; ++v) {
            dudt[i][v] += wc[i] * q[v];
        }
    }

    const long n_inner_faces = eqns->mesh->n_inner_faces;
    const long n_faces = eqns->mesh->n_faces;
    const ALIAS(cell, eqns->mesh->face.cell);
    const ALIAS(xf, eqns->mesh->face.center);
    const ALIAS(wf, eqns->mesh->face.gauss_weight);

    for (long i = 0; i < n_inner_faces; ++i) {
        eqns->source(xf[i], time, 0, q);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cell[i][0]][v] += wf[i][0] * q[v];
            dudt[cell[i][1]][v] += wf[i][1] * q[v];
        }
    }

    for (long i = n_inner_faces; i < n_faces; ++i) {
        eqns->source(xf[i], time, 0, q);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cell[i][0]][v] += wf[i][0] * q[v];
        }
    }
}

static void integrate_local_fluxes(Equations *eqns)
{
    const long n_vars = eqns->vars.n_fields;
    const long n_inner_faces = eqns->mesh->n_inner_faces;
    const long n_bound_faces = eqns->mesh->n_bound_faces;
    const ALIAS(cell, eqns->mesh->face.cell);
    const ALIAS(normal, eqns->mesh->face.normal);
    const ALIAS(area, eqns->mesh->face.area);
    const FIELDS(u, eqns->vars);
    DERIVS(dudt, eqns->vars);
    double f[n_vars] = {};

    for (long i = 0; i < n_inner_faces; ++i) {
        eqns->flux(normal[i], u[cell[i][0]], u[cell[i][1]], f);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cell[i][0]][v] -= f[v] * area[i];
            dudt[cell[i][1]][v] += f[v] * area[i];
        }
    }

    for (long i = n_inner_faces; i < n_inner_faces + n_bound_faces; ++i) {
        eqns->flux(normal[i], u[cell[i][0]], u[cell[i][1]], f);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cell[i][0]][v] -= f[v] * area[i];
        }
    }
}

static void integrate_non_local_fluxes(Equations *eqns)
{
    const long n_vars = eqns->vars.n_fields;
    const long n_inner_faces = eqns->mesh->n_inner_faces;
    const long n_bound_faces = eqns->mesh->n_bound_faces;
    const long n_faces = eqns->mesh->n_faces;
    const ALIAS(cell, eqns->mesh->face.cell);
    const ALIAS(normal, eqns->mesh->face.normal);
    const ALIAS(area, eqns->mesh->face.area);
    const FIELDS(u, eqns->vars);
    DERIVS(dudt, eqns->vars);
    double f[n_vars] = {};

    for (long i = n_inner_faces + n_bound_faces; i < n_faces; ++i) {
        eqns->flux(normal[i], u[cell[i][0]], u[cell[i][1]], f);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cell[i][0]][v] -= f[v] * area[i];
        }
    }
}

static void compute_local_gradients(Equations *eqns)
{
    const long n_vars = eqns->vars.n_fields;
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const long n_inner_faces = eqns->mesh->n_inner_faces;
    const long n_bound_faces = eqns->mesh->n_bound_faces;
    const ALIAS(cell, eqns->mesh->face.cell);
    const ALIAS(weight, eqns->mesh->face.gradient_weight);
    const FIELDS(u, eqns->vars);
    GRADS(dudx, eqns->vars);

    for (long i = 0; i < n_inner_cells; ++i)
        for (long v = 0; v < n_vars; ++v)
            for (long d = 0; d < N_DIMS; ++d) dudx[i][v][d] = 0;

    for (long i = 0; i < n_inner_faces; ++i) {
        for (long v = 0; v < n_vars; ++v) {
            const double du = u[cell[i][1]][v] - u[cell[i][0]][v];
            for (long d = 0; d < N_DIMS; ++d) {
                dudx[cell[i][0]][v][d] += weight[i][0][d] * du;
                dudx[cell[i][1]][v][d] -= weight[i][1][d] * du;
            }
        }
    }

    for (long i = n_inner_faces; i < n_inner_faces + n_bound_faces; ++i) {
        for (long v = 0; v < n_vars; ++v) {
            const double du = u[cell[i][1]][v] - u[cell[i][0]][v];
            for (long d = 0; d < N_DIMS; ++d) {
                dudx[cell[i][0]][v][d] += weight[i][0][d] * du;
            }
        }
    }
}

static void compute_non_local_gradients(Equations *eqns)
{
    const long n_vars = eqns->vars.n_fields;
    const long n_inner_faces = eqns->mesh->n_inner_faces;
    const long n_bound_faces = eqns->mesh->n_bound_faces;
    const long n_faces = eqns->mesh->n_faces;
    const ALIAS(cell, eqns->mesh->face.cell);
    const ALIAS(weight, eqns->mesh->face.gradient_weight);
    const FIELDS(u, eqns->vars);
    GRADS(dudx, eqns->vars);

    for (long i = n_inner_faces + n_bound_faces; i < n_faces; ++i) {
        for (long v = 0; v < n_vars; ++v) {
            const double du = u[cell[i][1]][v] - u[cell[i][0]][v];
            for (long d = 0; d < N_DIMS; ++d) {
                dudx[cell[i][0]][v][d] += weight[i][0][d] * du;
                dudx[cell[i][1]][v][d] -= weight[i][1][d] * du;
            }
        }
    }
}

static void limit_gradients(Equations *eqns)
{
    const long n_vars = eqns->vars.n_fields;
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const ALIAS(i_cell, eqns->mesh->cell.i_cell);
    const ALIAS(cell, eqns->mesh->cell.cell);
    const ALIAS(dx, eqns->mesh->cell.reconstruction);
    const ALIAS(eps2, eqns->buf);
    const FIELDS(u, eqns->vars);
    GRADS(dudx, eqns->vars);

    for (long i = 0; i < n_inner_cells; ++i) {
        for (long v = 0; v < n_vars; ++v) {
            double umin = u[i][v];
            double umax = u[i][v];
            for (long j = i_cell[i]; j < i_cell[i + 1]; ++j) {
                umin = MIN(umin, u[cell[j]][v]);
                umax = MAX(umax, u[cell[j]][v]);
            }
            const long nj = i_cell[i + 1] - i_cell[i];
            const double psi =
                eqns->limiter(eps2[i], u[i][v], umin, umax, dudx[i][v], &dx[i_cell[i]], nj);
            for (long d = 0; d < N_DIMS; ++d) dudx[i][v][d] *= psi;
        }
    }
}

static void integrate_local_reconstructed_fluxes(Equations *eqns)
{
    const long n_vars = eqns->vars.n_fields;
    const long n_inner_faces = eqns->mesh->n_inner_faces;
    const long n_bound_faces = eqns->mesh->n_bound_faces;
    const ALIAS(cell, eqns->mesh->face.cell);
    const ALIAS(normal, eqns->mesh->face.normal);
    const ALIAS(area, eqns->mesh->face.area);
    const ALIAS(dx, eqns->mesh->face.reconstruction);
    const FIELDS(u, eqns->vars);
    const GRADS(dudx, eqns->vars);
    DERIVS(dudt, eqns->vars);
    double ul[n_vars], ur[n_vars], f[n_vars] = {};

    for (long i = 0; i < n_inner_faces; ++i) {
        for (long v = 0; v < n_vars; ++v) {
            ul[v] = u[cell[i][0]][v];
            ur[v] = u[cell[i][1]][v];
            for (long d = 0; d < N_DIMS; ++d) {
                ul[v] += dudx[cell[i][0]][v][d] * dx[i][0][d];
                ur[v] += dudx[cell[i][1]][v][d] * dx[i][1][d];
            }
        }
        eqns->flux(normal[i], ul, ur, f);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cell[i][0]][v] -= f[v] * area[i];
            dudt[cell[i][1]][v] += f[v] * area[i];
        }
    }

    for (long i = n_inner_faces; i < n_inner_faces + n_bound_faces; ++i) {
        for (long v = 0; v < n_vars; ++v) {
            ul[v] = u[cell[i][0]][v];
            for (long d = 0; d < N_DIMS; ++d) {
                ul[v] += dudx[cell[i][0]][v][d] * dx[i][0][d];
            }
        }
        eqns->flux(normal[i], ul, u[cell[i][1]], f);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cell[i][0]][v] -= f[v] * area[i];
        }
    }
}

static void integrate_non_local_reconstructed_fluxes(Equations *eqns)
{
    const long n_vars = eqns->vars.n_fields;
    const long n_inner_faces = eqns->mesh->n_inner_faces;
    const long n_bound_faces = eqns->mesh->n_bound_faces;
    const long n_faces = eqns->mesh->n_faces;
    const ALIAS(cell, eqns->mesh->face.cell);
    const ALIAS(normal, eqns->mesh->face.normal);
    const ALIAS(area, eqns->mesh->face.area);
    const ALIAS(dx, eqns->mesh->face.reconstruction);
    const FIELDS(u, eqns->vars);
    const GRADS(dudx, eqns->vars);
    DERIVS(dudt, eqns->vars);
    double ul[n_vars], ur[n_vars], f[n_vars] = {};

    for (long i = n_inner_faces + n_bound_faces; i < n_faces; ++i) {
        for (long v = 0; v < n_vars; ++v) {
            ul[v] = u[cell[i][0]][v];
            ur[v] = u[cell[i][1]][v];
            for (long d = 0; d < N_DIMS; ++d) {
                ul[v] += dudx[cell[i][0]][v][d] * dx[i][0][d];
                ur[v] += dudx[cell[i][1]][v][d] * dx[i][1][d];
            }
        }
        eqns->flux(normal[i], ul, ur, f);
        for (long v = 0; v < n_vars; ++v) {
            dudt[cell[i][0]][v] -= f[v] * area[i];
            dudt[cell[i][1]][v] += f[v] * area[i];
        }
    }
}

static void scale_derivative(Equations *eqns)
{
    const long n_vars = eqns->vars.n_fields;
    const long n_inner_cells = eqns->mesh->n_inner_cells;
    const ALIAS(volume, eqns->mesh->cell.volume);
    DERIVS(dudt, eqns->vars);

    for (long i = 0; i < n_inner_cells; ++i)
        for (long v = 0; v < n_vars; ++v) dudt[i][v] /= volume[i];
}
