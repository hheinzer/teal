#include <math.h>

#include "euler.h"
#include "mesh/find.h"
#include "teal/h5io.h"
#include "teal/memory.h"
#include "teal/option.h"
#include "teal/print.h"
#include "teal/sync.h"
#include "teal/utils.h"

static void extract(const Equations *eqns, Vector3d **x, Vector3d **n, double **p, double **a,
                    long E, long n_faces);

static void compute(double **cp, double *cl, double *cd, const Vector3d *n, const double *p,
                    const double *a, const double *uref, double lref, long n_faces);

void euler_polar(const Simulation *sim, const char *airfoil, const double *uref, double lref)
{
    const alias(j_face, sim->eqns->mesh->entity.j_face);
    const long E = mesh_find_entity(sim->eqns->mesh, airfoil);
    const long n_faces = j_face[E + 1] - j_face[E];

    smart Vector3d *x, *n;
    smart double *p, *a;
    extract(sim->eqns, &x, &n, &p, &a, E, n_faces);

    smart double *cp;
    double cl = 0, cd = 0;
    compute(&cp, &cl, &cd, n, p, a, uref, lref, n_faces);

    if (sync.rank == 0 && !option.quiet) {
        print_key("lift coefficient", "%g", cl);
        print_key("drag coefficient", "%g", cd);
    }

    String fname;
    sprintf(fname, "%s_polar.hdf", sim->prefix);
    hid_t file = h5_file_create(fname);
    h5_dataset_write("x", *x, h5_dims(n_faces, N_DIMS), file);
    h5_dataset_write("cp", cp, h5_dims(n_faces), file);
    h5_dataset_write("cl", &cl, h5_dims(sync.rank == 0), file);
    h5_dataset_write("cd", &cd, h5_dims(sync.rank == 0), file);
    h5_file_close(file);
}

static void extract(const Equations *eqns, Vector3d **x, Vector3d **n, double **p, double **a,
                    long E, long n_faces)
{
    const long n_vars = eqns->n_vars;
    const alias(j_face, eqns->mesh->entity.j_face);
    const alias(center, eqns->mesh->face.center);
    const alias(area, eqns->mesh->face.area);
    const alias(basis, eqns->mesh->face.basis);
    const alias(cell, eqns->mesh->face.cell);
    const double(*u)[n_vars] = (void *)eqns->vars.u;

    *x = memory_calloc(n_faces, sizeof(**x));
    *n = memory_calloc(n_faces, sizeof(**n));
    *p = memory_calloc(n_faces, sizeof(**p));
    *a = memory_calloc(n_faces, sizeof(**a));
    for (long m = 0, j = j_face[E]; j < j_face[E + 1]; ++j) {
        const long cL = cell[j][L], cR = cell[j][R];
        for (long d = 0; d < N_DIMS; ++d) (*x)[m][d] = center[j][d];
        for (long d = 0; d < N_DIMS; ++d) (*n)[m][d] = basis[j][X][d];
        (*p)[m] = 0.5 * (u[cL][P] + u[cR][P]);
        (*a)[m] = area[j];
        m += 1;
    }
}

static void compute(double **cp, double *cl, double *cd, const Vector3d *n, const double *p,
                    const double *a, const double *uref, double lref, long n_faces)
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
