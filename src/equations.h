#pragma once

#include "matrix.h"
#include "mesh.h"

typedef struct Equations Equations;

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

// Convert from one variable representation to another.
typedef void Convert(void *, const void *, const double *property);

// Pre-pass called once per derivative evaluation before per-cell source terms.
typedef void Prepare(const Equations *, const void *primitive);

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
    union {
        Convert *convert;
        Compute *compute;
    } func;
} EquationsVariables;

typedef struct {
    int num;
    String *name;
    double *data;
} EquationsProperties;

typedef struct {
    Compute *compute;
    Prepare *prepare;
} EquationsSource;

struct Equations {
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
};

// Create an empty equation system.
Equations *equations_create(const Mesh *mesh, const char *name, Timestep *timestep,
                            ConvectiveSelect *convective, ViscousSelect *viscous,
                            BoundarySelect *boundary);

// Register `num` primitive variables with names, dimensions, and a conserved-to-primitive callback.
void equations_create_primitive(Equations *eqns, const char **name, const int *dimension,
                                Convert *convert, int num);

// Register `num` conserved variables with names, dimensions, and a primitive-to-conserved callback.
void equations_create_conserved(Equations *eqns, const char **name, const int *dimension,
                                Convert *convert, int num);

// Register a reference solution callback for error norms.
void equations_create_reference(Equations *eqns, Compute *compute);

// Register `num` physical properties with names and values.
void equations_create_properties(Equations *eqns, const char **name, const double *property,
                                 int num);

// Set the spatial order of the reconstruction (1 or 2).
void equations_set_space_order(Equations *eqns, int space_order);

// Select the convective flux scheme by name.
void equations_set_convective_flux(Equations *eqns, const char *name);

// Select the viscous flux scheme by name.
void equations_set_viscous_flux(Equations *eqns, const char *name);

// Select the gradient limiter by name with a scaling coefficient.
void equations_set_limiter(Equations *eqns, const char *name, double coefficient);

// Set the boundary condition for an entity by name; pass either `name` or `compute`, not both.
void equations_set_boundary_condition(Equations *eqns, const char *entity, const char *name,
                                      Compute *compute, const void *context);

// Set the initial condition for an entity; pass either `compute` or `initial`, not both.
void equations_set_initial_condition(Equations *eqns, const char *entity, Compute *compute,
                                     const void *initial);

// Override a named physical property value.
void equations_set_property(Equations *eqns, const char *name, double property);

// Set the source term callback and optionally a prepare function.
void equations_set_source(Equations *eqns, Compute *compute, Prepare *prepare);

// Print a summary of the equation system.
void equations_summary(const Equations *eqns);

// Compute the global stable time step over all inner cells.
double equations_timestep(const Equations *eqns, const void *primitive);

// Apply boundary conditions to primitive variables.
void equations_boundary(const Equations *eqns, void *primitive, double time);

// Compute cell-centered gradients of primitive variables.
void equations_gradient(const Equations *eqns, void *primitive, void *gradient);

// Apply the selected gradient limiter.
void equations_limiter(const Equations *eqns, const void *primitive, void *gradient);

// Compute the time derivative of the conserved variables.
void equations_derivative(const Equations *eqns, void *primitive, void *derivative, double time);

// Compute the volume-weighted L2 norm of the derivative.
void equations_residual(const Equations *eqns, const void *derivative, void *residual);

// Compute the volume-weighted average of primitive variables over a named entity.
void equations_average(const Equations *eqns, const char *entity, const void *primitive,
                       void *average);

// Compute the L2 error norm against the reference solution.
void equations_norm(const Equations *eqns, void *norm, double time);

// Write primitive variables and optional reference solution to disk.
void equations_write(const Equations *eqns, const char *name, double time, int index);

// Release all equation system resources.
void equations_destroy(Equations *eqns);
