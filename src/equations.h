#pragma once

#include "mesh.h"

typedef void Update(void *variable_, const scalar *property);

typedef scalar Timestep(const void *variable_, const scalar *property, scalar volume,
                        vector projection);

typedef void Compute(void *variable_, const scalar *property, vector center, scalar time,
                     const void *ctx_);

typedef void Boundary(void *ghost_, const void *inner_, const void *reference_,
                      const scalar *property, const matrix *basis);

typedef Boundary *BoundarySelect(const char *name);

typedef void Convective(void *flux_, const void *left_, const void *right_, const scalar *property,
                        const matrix *basis);

typedef Convective *ConvectiveSelect(const char *name);

typedef void Viscous(void *flux_, const void *variable_, const void *gradient_,
                     const scalar *property, const matrix *basis);

typedef Viscous *ViscousSelect(const char *name);

typedef void Limiter(vector *gradient, scalar variable, scalar minimum, scalar maximum,
                     scalar parameter, const vector *offset, long num);

typedef void Source(void *source_, const void *variable_, const scalar *property, vector center,
                    scalar time);

typedef struct {
    long num;
    long len;
    long stride;
    long *size;
    Name *name;
    void *data;
    Update *conserved;
    Update *primitive;
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
    long stride;
    long *size;
    Name *name;
    Compute *compute;
    Update *conserved;
} EquationsUserVariables;

typedef struct {
    const Mesh *mesh;
    Name name;
    long space_order;
    EquationsVariables variables;
    EquationsProperties properties;
    Timestep *timestep;
    EquationsBoundary boundary;
    EquationsConvective convective;
    EquationsViscous viscous;
    EquationsLimiter limiter;
    EquationsUserVariables user;
    Source *source;
} Equations;

Equations *equations_create(const Mesh *mesh, const char *name);

void equations_create_variables(Equations *eqns, const long *size, const char **name,
                                Update *conserved, Update *primitive, long num_conserved, long num);

void equations_create_properties(Equations *eqns, const char **name, const scalar *property,
                                 long num);

void equations_create_user_variables(Equations *eqns, const long *size, const char **name,
                                     Compute *compute, long num);

void equations_create_exact_solution(Equations *eqns, Compute *compute);

void equations_set_space_order(Equations *eqns, long space_order);

void equations_set_timestep(Equations *eqns, Timestep *timestep);

void equations_set_boundary_select(Equations *eqns, BoundarySelect *select);

void equations_set_convective_select(Equations *eqns, ConvectiveSelect *select);

void equations_set_viscous_select(Equations *eqns, ViscousSelect *select);

void equations_set_limiter(Equations *eqns, const char *name, scalar parameter);

void equations_set_boundary_condition(Equations *eqns, const char *entity, const char *name,
                                      const void *reference, Compute *custom);

void equations_set_convective_flux(Equations *eqns, const char *name);

void equations_set_viscous_flux(Equations *eqns, const char *name);

void equations_set_initial_condition(Equations *eqns, const char *entity, Compute *compute,
                                     scalar time);

void equations_set_initial_state(Equations *eqns, const char *entity, const void *state);

void equations_set_property(Equations *eqns, const char *name, scalar property);

void equations_set_user_source(Equations *eqns, Source *source);

void equations_summary(const Equations *eqns);

void equations_restart(const Equations *eqns, scalar *time, long *index);

scalar equations_timestep(const Equations *eqns, const void *variable_, scalar *step);

void equations_boundary(const Equations *eqns, void *variable_, scalar time);

void *equations_derivative(const Equations *eqns, void *variable_, void *derivative_, scalar time);

void *equations_gradient(const Equations *eqns, void *variable_);

void equations_limiter(const Equations *eqns, const void *variable_, void *gradient_);

void equations_residual(const Equations *eqns, const void *derivative_, void *residual_);

void equations_average(const Equations *eqns, const char *entity, void *average_);

void equations_norm(const Equations *eqns, scalar time, void *norm_);

void equations_write(const Equations *eqns, const char *prefix, scalar time, long index);
