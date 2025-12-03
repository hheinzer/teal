// Equations couples discretized PDEs to a mesh and wires callbacks used during time integration.
//
// Components:
// - Variables: cell-wise state vectors (conserved first, then primitive)
// - Properties: global scalar parameters passed into all callbacks
// - Time stepper: per-cell CFL-like restriction
// - Boundary conditions: ghost-cell update rules per mesh ghost entity
// - Fluxes: convective and viscous face flux implementations
// - Limiter: gradient limiting on reconstructed slopes
// - User variables: derived fields (exact solution, diagnostics, etc.)
// - Source term: optional cell-wise source contribution
//
// Conserved variables must precede primitive ones; if a variable is both, place it in the conserved
// set. All callbacks operate on plain scalar/vector arrays; Equations only provides the common
// plumbing to the mesh and to the time integration driver.
#pragma once

#include "mesh.h"

// Update the variables of a single cell in-place; may switch between conserved/primitive forms or
// enforce algebraic constraints.
typedef void Update(void *variable_, const scalar *property);

// Compute the local stability time step for a single cell; return a non-negative finite value. The
// global time step is the minimum over all cells.
typedef scalar TimeStep(const void *variable_, const scalar *property, scalar volume,
                        vector projection);

// Compute cell variables at (center, time) with optional context: ignored for initial conditions,
// current state for user variables, inner-cell state for custom boundary conditions.
typedef void Compute(void *variable_, const scalar *property, vector center, scalar time,
                     const void *ctx_);

// Compute the ghost-cell state from the inner-cell state and an optional reference.
typedef void Boundary(void *ghost_, const void *inner_, const void *reference_,
                      const scalar *property, const Basis *basis);

// Select a boundary condition implementation by name.
typedef Boundary *BoundarySelect(const char *name);

// Compute the convective flux across a face given the left and right states.
typedef void Convective(void *flux_, const void *left_, const void *right_, const scalar *property,
                        const Basis *basis);

// Select a convective flux implementation by name.
typedef Convective *ConvectiveSelect(const char *name);

// Compute the viscous flux across a face given the average cell state and gradient.
typedef void Viscous(void *flux_, const void *variable_, const void *gradient_,
                     const scalar *property, const Basis *basis);

// Select a viscous flux implementation by name.
typedef Viscous *ViscousSelect(const char *name);

// Limit a cell gradient in-place so reconstructed face values stay within [minimum, maximum], using
// an optional limiter parameter and cell-to-face offsets.
typedef void Limiter(vector *gradient, scalar variable, scalar minimum, scalar maximum,
                     scalar parameter, const vector *offset, long num);

// Compute the source term in a single cell at (center, time) and optionally the cell state.
typedef void Source(void *source_, const void *variable_, const scalar *property, vector center,
                    scalar time);

typedef struct {
    long num;
    long len;     // sum of components for conserved variables
    long stride;  // sum of components for all variables
    long *dim;    // components per variable
    Name *name;
    void *data;
    Update *conserved;  // called after initial and boundary conditions
    Update *primitive;  // called after advancing conserved variables in time
} EquationsVariables;

typedef struct {
    long num;
    Name *name;
    scalar *data;
} EquationsProperties;

typedef struct {
    long num;
    Name *name;
    const void **reference;
    Compute **custom;
    Boundary **condition;
    BoundarySelect *select;
} EquationsBoundary;

typedef struct {
    Name name;
    Convective *flux;
    ConvectiveSelect *select;
} EquationsConvective;

typedef struct {
    Name name;
    Viscous *flux;
    ViscousSelect *select;
} EquationsViscous;

typedef struct {
    Name name;
    scalar *parameter;
    Limiter *compute;
} EquationsLimiter;

typedef struct {
    long num;
    long stride;  // sum of components for all variables
    long *dim;    // components per variable
    Name *name;
    Compute *compute;
    Update *conserved;  // called after initial and boundary conditions (if set)
} EquationsUserVariables;

typedef struct {
    const Mesh *mesh;
    Name name;
    long space_order;
    EquationsVariables variables;
    EquationsProperties properties;
    TimeStep *time_step;
    EquationsBoundary boundary;
    EquationsConvective convective;
    EquationsViscous viscous;
    EquationsLimiter limiter;
    EquationsUserVariables user;
    Source *source;
} Equations;

// Create an empty equation system based on a mesh.
Equations *equations_create(const Mesh *mesh, const char *name);

// Create the primary variables for all mesh cells, initialized to zero.
void equations_create_variables(Equations *eqns, const long *dim, const char **name,
                                Update *conserved, Update *primitive, long num_conserved, long num);

// Create global scalar material properties.
void equations_create_properties(Equations *eqns, const char **name, const scalar *property,
                                 long num);

// Create the user variables for all mesh cells, only computed during `equations_write()`.
void equations_create_user_variables(Equations *eqns, const long *dim, const char **name,
                                     Compute *compute, long num);

// Create the user variables as an exact solution for the primary variables.
void equations_create_exact_solution(Equations *eqns, Compute *compute);

// Set the spatial order of the reconstruction (1 or 2).
void equations_set_space_order(Equations *eqns, long space_order);

// Set the local time step function used by `equations_time_step()`.
void equations_set_time_step(Equations *eqns, TimeStep *time_step);

// Set the boundary condition factory used by `equations_set_boundary_condition()`.
void equations_set_boundary_select(Equations *eqns, BoundarySelect *select);

// Set the convective flux factory used by `equations_set_convective_flux()`.
void equations_set_convective_select(Equations *eqns, ConvectiveSelect *select);

// Set the viscous flux factory used by `equations_set_viscous_flux()`.
void equations_set_viscous_select(Equations *eqns, ViscousSelect *select);

// Set the gradient limiter function used by `equations_limiter()`, optionally with a parameter.
void equations_set_limiter(Equations *eqns, const char *name, scalar parameter);

// Assign a boundary condition to a ghost entity, optionally with a reference or custom function.
void equations_set_boundary_condition(Equations *eqns, const char *entity, const char *name,
                                      const void *reference, Compute *custom);

// Select the convective flux implementation.
void equations_set_convective_flux(Equations *eqns, const char *name);

// Select the viscous flux implementation.
void equations_set_viscous_flux(Equations *eqns, const char *name);

// Set the initial condition of all cells of an inner entity by evaluating `compute`.
void equations_set_initial_condition(Equations *eqns, const char *entity, Compute *compute,
                                     scalar time);

// Set the initial condition of all cells of an inner entity to a constant state.
void equations_set_initial_state(Equations *eqns, const char *entity, const void *state);

// Set the value of a material property.
void equations_set_property(Equations *eqns, const char *name, scalar property);

// Set a user-defined source term callback.
void equations_set_user_source(Equations *eqns, Source *source);

// Print a summary of the equation system.
void equations_summary(const Equations *eqns);

// Restart the equation system and update time and index from the restart state.
void equations_restart(const Equations *eqns, scalar *time, long *index);

// Compute the global stable time step.
scalar equations_time_step(const Equations *eqns, const void *variable_);

// Apply boundary conditions to all ghost cells at the specified time.
void equations_boundary(const Equations *eqns, void *variable_, scalar time);

// Compute the time derivative of the conserved variables.
void *equations_derivative(const Equations *eqns, void *variable_, void *derivative_, scalar time);

// Compute the gradients of all variables in all cells.
void *equations_gradient(const Equations *eqns, void *variable_);

// Apply the selected limiter function to the gradients.
void equations_limiter(const Equations *eqns, const void *variable_, void *gradient_);

// Compute the residual from the time derivative.
void equations_residual(const Equations *eqns, const void *derivative_, void *residual_);

// Compute the average of the solution over a mesh entity.
void equations_average(const Equations *eqns, const char *entity, void *average_);

// Compute the norm of the solution at the specified time.
void equations_norm(const Equations *eqns, scalar time, void *norm_);

// Write the current solution to a VTKHDF file.
void equations_write(const Equations *eqns, const char *prefix, scalar time, long index);
