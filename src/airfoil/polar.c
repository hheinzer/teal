#include <math.h>
#include <string.h>

#include "airfoil.h"
#include "core/h5io.h"
#include "core/memory.h"
#include "core/sync.h"
#include "core/utils.h"
#include "teal.h"

static long find_entities(const Mesh *mesh, const char *airfoil);

static void extract(const Equations *eqns, double (**x)[N_DIMS], double **p, double **area,
                    double (**n)[N_DIMS], long E, long n_faces, long P);

static void compute(double **cp, double *cl, double *cd, long n_faces, const double *p,
                    const double *a, const double (*n)[N_DIMS], const double *uref, double lref,
                    long D, long U, long V, long W, long P);

void airfoil_polar(const Simulation *sim, const char *airfoil, const double *uref, double lref,
                   long D, long U, long V, long W, long P)
{
    const ALIAS(j_face, sim->eqns->mesh->entity.j_face);
    const long E = find_entities(sim->eqns->mesh, airfoil);
    const long n_faces = j_face[E + 1] - j_face[E];

    smart double(*x)[N_DIMS], *p, *a, (*n)[N_DIMS];
    extract(sim->eqns, &x, &p, &a, &n, E, n_faces, P);

    smart double *cp;
    double cl = 0, cd = 0;
    compute(&cp, &cl, &cd, n_faces, p, a, n, uref, lref, D, U, V, W, P);

    if (teal.rank == 0) {
        printf(" | " KEYFMT ": %g\n", "lift coefficient", cl);
        printf(" | " KEYFMT ": %g\n", "drag coefficient", cd);
    }

    char fname[128];
    sprintf(fname, "%s_polar.hdf", sim->prefix);
    hid_t file = h5_file_create(fname);
    h5_dataset_write("x", *x, H5DIMS(n_faces, N_DIMS), file);
    h5_dataset_write("cp", cp, H5DIMS(n_faces), file);
    h5_dataset_write("cl", &cl, H5DIMS(teal.rank == 0), file);
    h5_dataset_write("cd", &cd, H5DIMS(teal.rank == 0), file);
    h5_file_close(file);
}

static long find_entities(const Mesh *mesh, const char *airfoil)
{
    long E = -1;
    for (long e = 0; e < mesh->n_entities; ++e)
        if (!strcmp(mesh->entity.name[e], airfoil)) E = e;
    if (E == -1) error("could not find entity '%s'", airfoil);
    return E;
}

static void extract(const Equations *eqns, double (**x)[N_DIMS], double **p, double **a,
                    double (**n)[N_DIMS], long E, long n_faces, long P)
{
    const long n_vars = eqns->n_vars;
    const ALIAS(j_face, eqns->mesh->entity.j_face);
    const ALIAS(center, eqns->mesh->face.center);
    const ALIAS(area, eqns->mesh->face.area);
    const ALIAS(normal, eqns->mesh->face.normal);
    const ALIAS(cell, eqns->mesh->face.cell);
    double(*u)[n_vars] = (void *)eqns->vars.u;

    *x = memory_calloc(n_faces, sizeof(**x));
    *p = memory_calloc(n_faces, sizeof(**p));
    *a = memory_calloc(n_faces, sizeof(**a));
    *n = memory_calloc(n_faces, sizeof(**n));
    for (long m = 0, j = j_face[E]; j < j_face[E + 1]; ++j) {
        for (long d = 0; d < N_DIMS; ++d) (*x)[m][d] = center[j][d];
        (*p)[m] = 0.5 * (u[cell[j][L]][P] + u[cell[j][R]][P]);
        (*a)[m] = area[j];
        for (long d = 0; d < N_DIMS; ++d) (*n)[m][d] = normal[j][d];
        m += 1;
    }
}

static void compute(double **cp, double *cl, double *cd, long n_faces, const double *p,
                    const double *a, const double (*n)[N_DIMS], const double *uref, double lref,
                    long D, long U, long V, long W, long P)
{
    const double Vref = sqrt(sq(uref[U]) + sq(uref[V]) + sq(uref[W]));
    const double qref = 0.5 * uref[D] * sq(Vref);
    double c[2] = {0};

    *cp = memory_calloc(n_faces, sizeof(**cp));
    for (long i = 0; i < n_faces; ++i) {
        (*cp)[i] = (p[i] - uref[P]) / qref;
        c[0] += (p[i] * a[i] * n[i][X]) / (qref * lref);
        c[1] += (p[i] * a[i] * n[i][Y]) / (qref * lref);
    }
    for (long i = 0; i < 2; ++i) c[i] = sync_sum(c[i]);

    const double alpha = atan2(uref[V], uref[U]);
    *cl += cos(alpha) * c[1] - sin(alpha) * c[0];
    *cd += cos(alpha) * c[0] + sin(alpha) * c[1];
}
