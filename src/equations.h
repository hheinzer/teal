#pragma once

#include <mpi.h>

#include "mesh.h"

typedef struct Equations Equations;

typedef double TimeStep(const Equations *eqns, const double *u, const Vector3d projection,
                        double volume);

typedef void Update(const Equations *eqns, double *u);

typedef void ConvFlux(const Equations *eqns, const Matrix3d b, const double *ul, const double *ur,
                      double *f);

typedef void ViscFlux(const Equations *eqns, const Matrix3d b, const double *um,
                      const Vector3d *dudxm, double *f);

typedef double Limiter(const Vector3d *dx, const Vector3d dudx, double u, double umin, double umax,
                       double eps2, long n);

typedef void ApplyBC(const Equations *eqns, const Matrix3d b, const double *state, const double *ui,
                     double *ug);

typedef ConvFlux *SelectConvFlux(const char *name);

typedef ViscFlux *SelectViscFlux(const char *name);

typedef ApplyBC *SelectApplyBC(const char *name);

typedef void Compute(double *, const double *, const Vector3d x, double time);

typedef void Prepare(const Equations *eqns);

struct Equations {
    const Mesh *mesh;
    String name;

    long n_cons, n_vars, n_scalars, n_user;
    long space_order;

    struct {
        TimeStep *compute;
        double *value;
    } timestep;

    struct {
        Update *boundary;
        Update *advance;
    } update;

    struct {
        String *name;
        double *u, *dudt;
        Vector3d *dudx;
        long *dim;
    } vars;

    struct {
        String *name;
        double *value;
    } scalar;

    struct {
        String *name;
        Compute *compute;
        long *dim;
    } user;

    struct {
        SelectConvFlux *select;
        String name;
        ConvFlux *flux;
    } conv;

    struct {
        SelectViscFlux *select;
        String name;
        ViscFlux *flux;
    } visc;

    struct {
        String name;
        Limiter *limiter;
        double k, *eps2;
    } limiter;

    struct {
        Prepare *prepare;
        Compute *compute;
    } source;

    struct {
        SelectApplyBC *select;
        String *name;
        const double **state;
        ApplyBC **apply;
        Compute **compute;
    } bc;

    struct {
        double *buf;
        MPI_Request *recv, *send;
    } sync;
};

Equations *equations_create(const Mesh *mesh, const char *name);

void equations_create_vars(Equations *eqns, const String *name, long n_cons, long n_vars);

void equations_create_scalar(Equations *eqns, const String *name, long n_scalars);

void equations_create_user(Equations *eqns, const String *name, long n_user, Compute *compute);

void equations_create_exact(Equations *eqns, Compute *compute);

void equations_set_timestep(Equations *eqns, TimeStep *compute);

void equations_set_update(Equations *eqns, Update *boundary, Update *advance);

void equations_set_select_convective_flux(Equations *eqns, SelectConvFlux *select);

void equations_set_select_viscous_flux(Equations *eqns, SelectViscFlux *select);

void equations_set_select_boundary_condition(Equations *eqns, SelectApplyBC *select);

void equations_set_scalar(Equations *eqns, long name, double value);

void equations_set_convective_flux(Equations *eqns, const char *name);

void equations_set_viscous_flux(Equations *eqns, const char *name);

void equations_set_source(Equations *eqns, Compute *compute, Prepare *prepare);

void equations_set_space_order(Equations *eqns, long space_order);

void equations_set_limiter(Equations *eqns, const char *name, double k);

void equations_set_boundary_condition(Equations *eqns, const char *entity, const char *name,
                                      const double *state, Compute *compute);

void equations_set_initial_condition(Equations *eqns, Compute *compute);

void equations_set_initial_state(Equations *eqns, const double *state);

void equations_print(const Equations *eqns);

void equations_write(const Equations *eqns, const char *prefix, long count, double time);

double equations_timestep(Equations *eqns);

void equations_derivative(Equations *eqns, double time);

void equations_gradient(Equations *eqns);

void equations_limiter(Equations *eqns);

void equations_residual(const Equations *eqns, double *residual);

void equations_average(const Equations *eqns, double *average);

void equations_free(Equations **eqns);
