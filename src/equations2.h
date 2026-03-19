#pragma once

#include "matrix2.h"
#include "mesh2.h"

// Compute the maximum stable time step for a cell given its volume and face projection.
typedef double Timestep(const void *primitive, const double *property, double volume,
                        Vector projection);

// Compute the numerical flux at a face given left/right primitive states and face geometry.
typedef void Convective(void *flux, const void *left, const void *right, const double *property,
                        Matrix basis);

// Select a convective flux scheme by name.
typedef Convective *ConvectiveSelect(const char *name);

// Compute the viscous flux at a face given the primitive state, gradient, and face geometry.
typedef void Viscous(void *flux, const void *primitive, const void *gradient,
                     const double *property, Matrix basis);

// Select a viscous flux scheme by name.
typedef Viscous *ViscousSelect(const char *name);

// Compute a quantity given location, time, and optional context.
typedef void Compute(void *, const double *property, Vector center, double time,
                     const void *context);

// Compute the outer (ghost) primitive state from the inner state and face geometry.
typedef void Boundary(void *outer, const void *inner, const double *property, Matrix basis,
                      const void *context);

// Select a boundary condition by name.
typedef Boundary *BoundarySelect(const char *name);

typedef struct {
    Timestep *compute;
} EquationsTimestep;

typedef struct {
    String name;
    Convective *flux;
    ConvectiveSelect *select;
} EquationsConvective;

typedef struct {
    String name;
    Viscous *flux;
    ViscousSelect *select;
} EquationsViscous;

typedef struct {
    int kind;
    String name;
    double *parameter;
} EquationsLimiter;

typedef struct {
    int num;
    String *name;
    Compute **compute;
    Boundary **condition;
    const void **context;
    BoundarySelect *select;
} EquationsBoundary;

typedef struct {
    int num;
    int stride;
    String *name;
    int *dimension;
    void *data;
    Compute *compute;
} EquationsVariables;

typedef struct {
    int num;
    String *name;
    double *data;
} EquationsProperties;

typedef struct {
    Compute *compute;
} EquationsSource;

typedef struct {
    const Mesh *mesh;
    String name;
    int space_order;
    EquationsTimestep timestep;
    EquationsConvective convective;
    EquationsViscous viscous;
    EquationsLimiter limiter;
    EquationsBoundary boundary;
    EquationsVariables primitive;
    EquationsVariables conserved;
    EquationsVariables reference;
    EquationsProperties properties;
    EquationsSource source;
} Equations;

// Create an empty equation system.
Equations *equations2_create(const Mesh *mesh, const char *name, Timestep *timestep,
                             ConvectiveSelect *convective, ViscousSelect *viscous,
                             BoundarySelect *boundary);

// Register `num` primitive variables with names, dimensions, and an initial condition callback.
void equations2_create_primitive(Equations *eqns, const char **name, const int *dimension,
                                 Compute *compute, int num);

// Register `num` conserved variables with names, dimensions, and a conversion callback.
void equations2_create_conserved(Equations *eqns, const char **name, const int *dimension,
                                 Compute *compute, int num);

// Register a reference solution callback for error norms.
void equations2_create_reference(Equations *eqns, Compute *compute);

// Register `num` physical properties with names and values.
void equations2_create_properties(Equations *eqns, const char **name, const double *property,
                                  int num);

// Set the spatial order of the reconstruction (1 or 2).
void equations2_set_space_order(Equations *eqns, int space_order);

// Select the convective flux scheme by name.
void equations2_set_convective_flux(Equations *eqns, const char *name);

// Select the viscous flux scheme by name.
void equations2_set_viscous_flux(Equations *eqns, const char *name);

// Select the gradient limiter by name with a scaling coefficient.
void equations2_set_limiter(Equations *eqns, const char *name, double coefficient);

// Set the boundary condition for an entity by name; pass either `name` or `compute`, not both.
void equations2_set_boundary_condition(Equations *eqns, const char *entity, const char *name,
                                       Compute *compute, const void *context);

// Set the initial condition for an entity; pass either `compute` or `initial`, not both.
void equations2_set_initial_condition(Equations *eqns, const char *entity, Compute *compute,
                                      const void *initial);

// Override a named physical property value.
void equations2_set_property(Equations *eqns, const char *name, double property);

// Set the source term callback.
void equations2_set_source(Equations *eqns, Compute *compute);

// Print a summary of the equation system.
void equations2_summary(const Equations *eqns);

// Compute the global stable time step over all inner cells.
double equations2_timestep(const Equations *eqns, const void *primitive);

// Apply boundary conditions to primitive variables.
void equations2_boundary(const Equations *eqns, void *primitive, double time);

// Compute cell-centered gradients of primitive variables.
void equations2_gradient(const Equations *eqns, void *primitive, void *gradient);

// Apply the selected gradient limiter.
void equations2_limiter(const Equations *eqns, const void *primitive, void *gradient);

// Compute the time derivative of the conserved variables.
void equations2_derivative(const Equations *eqns, void *primitive, void *derivative, double time);

// Compute the volume-weighted L2 norm of the derivative.
void equations2_residual(const Equations *eqns, const void *derivative, void *residual);

// Compute the L2 error norm against the reference solution.
void equations2_norm(const Equations *eqns, void *norm, double time);

// Write primitive variables and optional reference solution to disk.
void equations2_write(const Equations *eqns, const char *name, double time, int index);

// Release all equation system resources.
void equations2_destroy(Equations *eqns);
