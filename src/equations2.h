#pragma once

#include "matrix2.h"
#include "mesh2.h"

typedef double Timestep(const void *primitive, const double *property, double volume,
                        Vector projection);

typedef void Convective(void *flux, const void *left, const void *right, const double *property,
                        Matrix basis);

typedef Convective *ConvectiveSelect(const char *name);

typedef void Viscous(void *flux, const void *primitive, const void *gradient,
                     const double *property, Matrix basis);

typedef Viscous *ViscousSelect(const char *name);

typedef void Compute(void *, const double *property, Vector center, double time,
                     const void *context);

typedef void Boundary(void *outer, const void *inner, const double *property, Matrix basis,
                      const void *context);

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

Equations *equations2_create(const Mesh *mesh, const char *name, Timestep *timestep,
                             BoundarySelect *boundary, ConvectiveSelect *convective,
                             ViscousSelect *viscous, int space_order);

void equations2_create_primitive(Equations *eqns, const char **name, const int *dimension,
                                 Compute *compute, int num);

void equations2_create_conserved(Equations *eqns, const char **name, const int *dimension,
                                 Compute *compute, int num);

void equations2_create_reference(Equations *eqns, Compute *compute);

void equations2_create_properties(Equations *eqns, const char **name, const double *property,
                                  int num);

void equations2_set_convective_flux(Equations *eqns, const char *name);

void equations2_set_viscous_flux(Equations *eqns, const char *name);

void equations2_set_limiter(Equations *eqns, const char *name, double coefficient);

void equations2_set_boundary_condition(Equations *eqns, const char *entity, const char *name,
                                       Compute *compute, const void *context);

void equations2_set_initial_condition(Equations *eqns, const char *entity, Compute *compute,
                                      const void *initial);

void equations2_set_property(Equations *eqns, const char *name, double property);

void equations2_set_source(Equations *eqns, Compute *compute);

void equations2_summary(const Equations *eqns);

double equations2_timestep(const Equations *eqns, const void *primitive);

void equations2_boundary(const Equations *eqns, void *primitive, double time);

void equations2_gradient(const Equations *eqns, void *primitive, void *gradient);

void equations2_limiter(const Equations *eqns, const void *primitive, void *gradient);

void equations2_derivative(const Equations *eqns, void *primitive, void *derivative, double time);

void equations2_residual(const Equations *eqns, const void *derivative, void *residual);

void equations2_norm(const Equations *eqns, void *norm, double time);

void equations2_write(const Equations *eqns, const char *name, double time, int index);

void equations2_destroy(Equations *eqns);
