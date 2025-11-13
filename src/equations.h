#pragma once

#include "mesh.h"

typedef enum {
    SCALAR = 1,
    VECTOR = 3,
    MATRIX = 9,
} Type;

typedef void Update(void *variable_, const scalar *property);

typedef struct {
    number num;
    number len;
    number stride;
    Name *name;
    Type *type;
    void *data;
    Update *conserved;
    Update *primitive;
} EquationsVariables;

typedef struct {
    number num;
    Name *name;
    scalar *data;
} EquationsProperties;

typedef scalar Timestep(const void *variable_, const scalar *property, scalar volume,
                        vector projection);

typedef void Compute(void *variable_, const void *context_, const scalar *property, vector center,
                     scalar time);

typedef void Boundary(void *ghost_, const void *inner_, const void *reference_,
                      const scalar *property, matrix basis);

typedef Boundary *BoundarySelect(const char *name);

typedef struct {
    number num;
    const Name *entity;
    const number *cell_off;
    const number *face_off;
    Name *name;
    const void **reference;
    Compute **custom;
    Boundary **condition;
    BoundarySelect *select;
} EquationsBoundary;

typedef void Convective(void *flux_, const void *left_, const void *right_, const scalar *property,
                        matrix basis);

typedef Convective *ConvectiveSelect(const char *name);

typedef struct {
    Name name;
    Convective *flux;
    ConvectiveSelect *select;
} EquationsConvective;

typedef void Viscous(void *flux_, const void *gradient_, const void *variable_,
                     const scalar *property, matrix basis);

typedef Viscous *ViscousSelect(const char *name);

typedef struct {
    Name name;
    Viscous *flux;
    ViscousSelect *select;
} EquationsViscous;

typedef scalar Limiter(vector gradient, scalar variable, scalar minimum, scalar maximum,
                       scalar parameter, const vector *offset, number beg, number end);

typedef struct {
    Name name;
    scalar *parameter;
    Limiter *compute;
} EquationsLimiter;

typedef struct {
    number num;
    number stride;
    Name *name;
    Type *type;
    Compute *compute;
    Update *conserved;
} EquationsUserVariables;

typedef void Source(void *source_, const void *variable_, const scalar *property, vector center,
                    scalar time);

typedef struct {
    const Mesh *mesh;
    Name name;
    number space_order;
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

void equations_create_variables(Equations *eqns, const char **name, const Type *type,
                                Update *conserved, Update *primitive, number num_conserved,
                                number num);

void equations_create_properties(Equations *eqns, const char **name, const scalar *property,
                                 number num);

void equations_create_user_variables(Equations *eqns, const char **name, const Type *type,
                                     Compute *compute, number num);

void equations_create_exact_solution(Equations *eqns, Compute *compute);

void equations_set_space_order(Equations *eqns, number space_order);

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

void equations_restart(const Equations *eqns, scalar *time, number *index);

scalar equations_timestep(const Equations *eqns, const void *variable_, scalar *value);

void equations_boundary(const Equations *eqns, void *variable_, scalar time);

void *equations_derivative(const Equations *eqns, void *variable_, scalar time);

void *equations_gradient(const Equations *eqns, void *variable_);

void equations_limiter(const Equations *eqns, const void *variable_, void *gradient_);

void equations_residual(const Equations *eqns, const void *derivative_, void *residual_);

void *equations_average(const Equations *eqns, const char *entity, const void *variable_);

void equations_write(const Equations *eqns, scalar time, const char *prefix, number index);
