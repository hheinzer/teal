#pragma once

#include <mpi.h>

#include "mesh.h"

typedef struct Equations Equations;

/* Compute a general function of space 'x' and 'time'. The first argument is the result and the
 * second argument is a context variable. */
typedef void Function(double *, const double *, const double *x, double time);

/* Compute the largest possible time step of a cell in the mesh of 'eqns', depending on the
 * solution variables 'u', the cell 'projection', and the cell 'volume'. */
typedef double TimeStep(const Equations *eqns, const double *u, const double *projection,
                        double volume);

/* Evaluate the convective flux across a face in the mesh of 'eqns', 'n' is the normal vector of
 * the face, 'ul' and 'ur' are the left and right state, and 'f' is the flux. */
typedef void ConvFlux(const Equations *eqns, const double *n, const double *ul, const double *ur,
                      double *f);

/* Evaluate the viscous flux across a face in the mesh of 'eqns', 'n' is the normal vector of the
 * face, 'um' is the average state at the face, 'dudxm' is the average state gradients at the face,
 * and 'f' is the flux. */
typedef void ViscFlux(const Equations *eqns, const double *n, const double *um,
                      const double (*dudxm)[N_DIMS], double *f);

/* Select a convective flux by its 'name'. */
typedef ConvFlux *SelectConvFlux(const char *name);

/* Select a viscous flux by its 'name'. */
typedef ViscFlux *SelectViscFlux(const char *name);

/* Apply a boundary condition at a ghost cell in the mesh of 'eqns', 'n' is the normal vector of
 * the face between the ghost cell and the corresponding inner cell, 'state' is a reference state,
 * 'ui' is the inner state and 'ug' is the ghost state. */
typedef void ApplyBC(const Equations *eqns, const double *n, const double *state, const double *ui,
                     double *ug);

/* Select a boundary condition by its 'name'. */
typedef ApplyBC *SelectBC(const char *name);

/* General purpose update function for the state 'u' of a cell in the mesh of 'eqns'. */
typedef void Update(const Equations *eqns, double *u);

/* Compute the limiter value for the gradient 'dudx', given the cell center to face center offsets
 * 'dx', the state variable 'u', the minimum 'umin' and maximum 'umax' of the state variable over
 * the cell and its 'n' neighbors, and the 'eps2' parameter of the cell. */
typedef double Limiter(const double (*dx)[N_DIMS], const double *dudx, double u, double umin,
                       double umax, double eps2, long n);

/* General purpose function that is called before the source term of 'eqns' is evaluated. */
typedef void Prepare(const Equations *eqns);

struct Equations {
    long n_vars, n_scalars, n_user;
    long space_order;
    char name[NAMELEN];
    const Mesh *mesh;

    TimeStep *timestep;
    Update *boundary, *advance;

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

    struct {
        char name[NAMELEN];
        Limiter *func;
        double k, *eps2;
    } limiter;

    struct {
        char (*name)[NAMELEN];
        Function *func;
        long *dim;
    } user;

    struct {
        Prepare *prepare;
        Function *func;
    } source;

    struct {
        SelectBC *select;
        char (*name)[NAMELEN];
        const double **state;
        ApplyBC **apply;
        Function **func;
    } bc;

    struct {
        double *buf;
        MPI_Request *recv, *send;
    } sync;
};

/* Create the equation system 'eqns' with the specified 'space_order'. This function should not be
 * called directly by the user, but rather from the creation function of a specific equation system
 * implementation. */
void equations_create(Equations *eqns, long space_order);

/* Deallocate memory of 'eqns' and set all fields to 0. */
void equations_free(Equations *eqns);

/* Create 'n_user' user variables for 'eqns' with the identifiers 'name'. The variables are
 * computed according to 'func'. */
void equations_create_user(Equations *eqns, const char **name, Function *func, long n_user);

/* Create 'n_user' exact user variables for 'eqns'. The names are taken from the solution variables
 * prefixed with "exact". The variables are computed according to 'func'. */
void equations_create_exact(Equations *eqns, Function *func, long n_user);

/* Set the 'scalar' of 'eqns' to 'value'. */
void equations_set_scalar(Equations *eqns, long scalar, double value);

/* Set the convective flux function of 'eqns' to 'name'. */
void equations_set_convective_flux(Equations *eqns, const char *name);

/* Set the viscous flux function of 'eqns' to 'name'. */
void equations_set_viscous_flux(Equations *eqns, const char *name);

/* Set the source 'prepare' function of 'eqns'. This function will be called before the source
 * function is evaluated. */
void equations_set_prepare(Equations *eqns, Prepare *prepare);

/* Set the source 'func' of 'eqns'. */
void equations_set_source(Equations *eqns, Function *func);

/* Set the limiter function of 'eqns' to 'name' and set the limiter constant 'k'. */
void equations_set_limiter(Equations *eqns, const char *name, double k);

/* Set the initial condition of 'eqns' according to the initial 'func'. */
void equations_set_initial_condition(Equations *eqns, Function *func);

/* Set the initial condition of 'eqns' uniformly to 'state'. */
void equations_set_initial_state(Equations *eqns, const double *state);

/* Set the boundary condition for 'entity' of 'eqns' to 'bc'. Depending on the selected boundary
 * condition 'state' and/or 'func' may be specified. */
void equations_set_boundary_condition(Equations *eqns, const char *entity, const char *bc,
                                      const double *state, Function *func);

/* Print a summary of 'eqns'. */
void equations_print(const Equations *eqns);

/* Write 'eqns' at 'time' and 'iter' to an HDF5 file in the VTKHDF format. The file name is based
 * on 'prefix' and 'count' determines the postfix of the file. */
void equations_write(const Equations *eqns, const char *prefix, long count, double time, long iter);

/* Compute the largest possible time step of 'eqns'. */
double equations_timestep(const Equations *eqns);

/* Compute the time derivative of 'eqns' at 'time'. */
void equations_derivative(Equations *eqns, double time);

/* Compute the spatial derivatives of 'eqns'. */
void equations_gradient(Equations *eqns);

/* Apply the limiter function to the spatial derivatives of 'eqns'. */
void equations_limiter(Equations *eqns);

/* Compute the volume weighted time derivative 'residual' of 'eqns'. */
void equations_residual(const Equations *eqns, double *residual);

/* Compute the volume 'avg' of the solution variables in 'eqns'. */
void equations_average(const Equations *eqns, double *avg);
