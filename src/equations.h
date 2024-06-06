#ifndef EQUATIONS_H
#define EQUATIONS_H

#include <mpi.h>

#include "mesh.h"
#include "utils.h"

#define FIELDS(a, f) double(*a)[f.n_fields] = TCAST(a, f.u)
#define GRADS(a, f) double(*a)[f.n_fields][N_DIMS] = TCAST(a, f.dudx)
#define DERIVS(a, f) double(*a)[f.n_fields] = TCAST(a, f.dudt)

typedef struct Fields {
    long n_fields;
    double *u, *dudx, *dudt;
    char **name;
    struct {
        long n_dims, *dim;
        char **name;
    } output;
    Function *compute;
} Fields;

typedef struct Equations Equations;
typedef double TimeStep(const Equations *eqns);
typedef void Update(double *u);
typedef void Flux(const double *n, const double *ul_g, const double *ur_g, double *f_g);
typedef double Limiter(const double eps2, const double u, const double umin, const double umax,
                       const double *dudx, const double (*dx)[N_DIMS], const long nj);
struct Equations {
    Fields vars, user;
    const Mesh *mesh;

    TimeStep *time_step;
    Update *update;
    Flux *flux;
    Function *source;

    long space_order;
    Limiter *limiter;
    double k, *buf;

    struct {
        double *buf_u, *buf_dudx;
        MPI_Request *recv_u, *send_u, *recv_dudx, *send_dudx;
    } sync;
};

Equations equations_create(const Mesh *mesh, const Fields *vars, const Fields *user);

void equations_free(Equations *eqns);

void equations_set_initial_condition(Equations *eqns, Function *initial);

void equations_set_initial_state(Equations *eqns, const long nu, const long *u,
                                 const double *state);

void equations_set_space_order(Equations *eqns, const long space_order, const char *limiter,
                               const double k);

void equations_print(const Equations *eqns, const char *name);

void equations_write(const Equations *eqns, const char *prefix, const long count,
                     const double time);

void equations_time_derivative(Equations *eqns, const double time);

void equations_residual(const Equations *eqns, double *residual);

#endif
