#pragma once

#include <mpi.h>

#include "mesh.h"

typedef void Function(double *, const double *u, const double *x, double time);

typedef struct Equations Equations;

typedef double TimeStep(const Equations *eqns, const double *u, const double *projection,
                        double volume);

typedef void ConvFlux(const Equations *eqns, const double *n, const double *g_ul,
                      const double *g_ur, double *g_f);

typedef void ViscFlux(const Equations *eqns, const double *n, const double *um,
                      const double (*dudxm)[N_DIMS], double *f);

typedef ConvFlux *SelectConvFlux(const char *name);

typedef ViscFlux *SelectViscFlux(const char *name);

typedef void ApplyBC(const Equations *eqns, const double *n, const double *state, const double *ui,
                     double *ug);

typedef ApplyBC *SelectBC(const char *name);

typedef void Update(const Equations *eqns, double *u);

typedef double Limiter(const double (*dx)[N_DIMS], const double *dudx, double u, double umin,
                       double umax, double eps2, long n);

struct Equations {
    long n_vars, n_scalars, n_user;
    long space_order;
    char name[NAMELEN];
    const Mesh *mesh;

    struct {
        char (*name)[NAMELEN];
        double *u, *dudx, *dudt;
        long *dim;
    } vars;

    struct {
        char (*name)[NAMELEN];
        double *value;
    } scalar;

    struct {
        SelectConvFlux *select_conv;
        SelectViscFlux *select_visc;
        char name_conv[NAMELEN], name_visc[NAMELEN];
        ConvFlux *conv;
        ViscFlux *visc;
    } flux;

    TimeStep *timestep;
    Update *boundary, *advance;

    struct {
        char name[NAMELEN];
        Limiter *func;
        double k, *eps2;
    } limiter;

    struct {
        char (*name)[NAMELEN];
        Function *compute;
        long *dim;
    } user;

    Function *source;

    struct {
        SelectBC *select;
        char (*name)[NAMELEN];
        const double **state;
        ApplyBC **apply;
        Function **custom;
    } bc;

    struct {
        double *buf;
        MPI_Request *recv, *send;
    } sync;
};

void equations_create(Equations *eqns, long space_order);

void equations_free(Equations *eqns);

void equations_create_user(Equations *eqns, Function *compute, const char **name, long n_user);

void equations_create_exact(Equations *eqns, Function *compute, long n_user);

void equations_set_scalar(Equations *eqns, long scalar, double value);

void equations_set_convective_flux(Equations *eqns, const char *name);

void equations_set_viscous_flux(Equations *eqns, const char *name);

void equations_set_source(Equations *eqns, Function *source);

void equations_set_limiter(Equations *eqns, const char *limiter, double k);

void equations_set_initial_condition(Equations *eqns, Function *initial);

void equations_set_initial_state(Equations *eqns, const double *state);

void equations_set_boundary_condition(Equations *eqns, const char *entity, const char *bc,
                                      const double *state, Function *custom);

void equations_print(const Equations *eqns);

void equations_write(const Equations *eqns, const char *prefix, long count, double time, long iter);

double equations_timestep(const Equations *eqns);

void equations_derivative(Equations *eqns, double time);

void equations_gradient(Equations *eqns);

void equations_limiter(Equations *eqns);

void equations_residual(const Equations *eqns, double *residual);
