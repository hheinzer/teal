#include "airfoil.h"

#include <assert.h>
#include <gmshc.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "equations.h"
#include "global.h"
#include "gmsh_io.h"
#include "hdf5_io.h"
#include "memory.h"
#include "mesh.h"
#include "simulation.h"
#include "sync.h"

Mesh airfoil_mesh(const char *fname, const long n_upper, const long n_lower, const long n_outer,
                  const double radius) {
    int ierr;

    gmshInitialize(0, 0, 0, 0, &ierr);
    gmshOptionSetNumber("General.Verbosity", 3, &ierr);
    gmshModelAdd("airfoil", &ierr);

    long n_points = 0;
    double point[1000][2];
    FILE *file = fopen(fname, "r");
    char line[256];
    while (fgets(line, sizeof(line), file)) {
        if (sscanf(line, " %lg %lg ", &point[n_points][X], &point[n_points][Y]) != 2) continue;
        n_points += 1;
    }
    fclose(file);
    assert(fabs(point[0][X] - point[n_points - 1][X]) < EPS &&
           fabs(point[0][Y] - point[n_points - 1][Y]) < EPS && "airfoild shape must be closed");
    n_points -= 1;

    long inose = 0;
    double dnose = sqrt(point[0][X] * point[0][X] + point[0][Y] * point[0][Y]);
    for (long i = 0; i < n_points; ++i) {
        gmshModelOccAddPoint(point[i][X], point[i][Y], 0, 1, i + 1, &ierr);
        const double d = sqrt(point[i][X] * point[i][X] + point[i][Y] * point[i][Y]);
        if (d < dnose) {
            dnose = d;
            inose = i;
        }
    }

    cleanup int *pointTags = memory_calloc(n_points, sizeof(*pointTags));
    for (long n = 0, i = 0; i <= inose; ++i) pointTags[n++] = i + 1;
    gmshModelOccAddSpline(pointTags, inose + 1, 1, 0, 0, &ierr);

    for (long n = 0, i = inose; i < n_points; ++i) pointTags[n++] = i + 1;
    pointTags[n_points - inose] = 1;
    gmshModelOccAddSpline(pointTags, n_points - inose + 1, 2, 0, 0, &ierr);

    const double pi = acos(-1);
    gmshModelOccAddCircle(0.5, 0, 0, radius, 3, 0, 2 * pi, 0, 0, 0, 0, &ierr);

    gmshModelOccAddCurveLoop((int[]){1, 2}, 2, 1, &ierr);
    gmshModelOccAddCurveLoop((int[]){3}, 1, 2, &ierr);
    gmshModelOccAddPlaneSurface((int[]){2, 1}, 2, 1, &ierr);

    gmshModelOccSynchronize(&ierr);

    gmshModelAddPhysicalGroup(1, (int[]){1, 2}, 2, 1, "airfoil", &ierr);
    gmshModelAddPhysicalGroup(1, (int[]){3}, 1, 2, "farfield", &ierr);
    gmshModelAddPhysicalGroup(2, (int[]){1}, 1, 1, "domain", &ierr);

    gmshModelMeshSetTransfiniteCurve(1, n_upper, "Progression", 1, &ierr);
    gmshModelMeshSetTransfiniteCurve(2, n_lower, "Progression", 1, &ierr);
    gmshModelMeshSetTransfiniteCurve(3, n_outer, "Progression", 1, &ierr);

    gmshModelMeshGenerate(2, &ierr);
    gmshModelMeshOptimize("Laplace2D", 0, 10, 0, 0, &ierr);

    Mesh mesh = {};
    gmsh_export(&mesh);
    mesh_finalize(&mesh);

    gmshFinalize(&ierr);
    return mesh;
}

void airfoil_polar(const Simulation *sim, const char *name, const long P, const double *uref,
                   const double lref) {
    // find entity
    const long n_entities = sim->eqns->mesh->n_entities;
    const ALIAS(ename, sim->eqns->mesh->entity.name);
    long E = -1;
    for (long e = 0; e < n_entities; ++e)
        if (!strcmp(ename[e], (name ? name : "airfoil"))) E = e;
    assert(E != -1 && "entity not found");

    // extract coordinates, area * normals, and pressure from faces
    const ALIAS(j_face, sim->eqns->mesh->entity.j_face);
    const ALIAS(cell, sim->eqns->mesh->face.cell);
    const ALIAS(center, sim->eqns->mesh->face.center);
    const ALIAS(area, sim->eqns->mesh->face.area);
    const ALIAS(normal, sim->eqns->mesh->face.normal);
    const FIELDS(u, sim->eqns->vars);
    const long n_faces = j_face[E + 1] - j_face[E];
    cleanup double(*x)[N_DIMS] = memory_calloc(n_faces * (1 + MAX_FACE_NODES), sizeof(*x));
    cleanup double *p = memory_calloc(n_faces * (1 + MAX_FACE_NODES), sizeof(*p));
    cleanup double(*an)[N_DIMS] = memory_calloc(n_faces * (1 + MAX_FACE_NODES), sizeof(*an));
    long n = 0;
    for (long j = j_face[E]; j < j_face[E + 1]; ++j) {
        for (long d = 0; d < N_DIMS; ++d) x[n][d] = center[j][d];
        p[n] = u[cell[j][0]][P];
        for (long d = 0; d < N_DIMS; ++d) an[n][d] = area[j] * normal[j][d];
        n += 1;
    }
    const long n_points = n;

    // compute pressure, lift, and drag coefficients
    const double vref = sqrt(uref[1] * uref[1] + uref[2] * uref[2]);
    const double qref = 0.5 * uref[0] * vref * vref;
    cleanup double *cp = memory_calloc(n_points, sizeof(*cp));
    double c[2] = {};
    for (long i = 0; i < n_points; ++i) {
        cp[i] = (p[i] - uref[3]) / qref;
        c[0] += (p[i] * an[i][1]) / (qref * lref);
        c[1] += (p[i] * an[i][0]) / (qref * lref);
    }
    sync_sum(c, 2);
    const double alpha = atan2(uref[2], uref[1]);
    const double cl = cos(alpha) * c[0] - sin(alpha) * c[1];
    const double cd = cos(alpha) * c[1] + sin(alpha) * c[0];

    const int rank = sim->eqns->mesh->rank;
    if (rank == 0) {
        printf(" | " FMT_KEY ": %g\n", "lift coefficient", cl);
        printf(" | " FMT_KEY ": %g\n", "drag coefficient", cd);
    }

    char fname[256];
    snprintf(fname, sizeof(fname), "%s_polar.h5", sim->prefix);
    hid_t file = hdf5_file_create(fname);
    hdf5_write_dataset(file, "x", *x, 2, HDF5_DIMS(n_points, N_DIMS));
    hdf5_write_dataset(file, "cp", cp, 1, HDF5_DIMS(n_points));
    hdf5_write_dataset(file, "cl", &cl, 1, HDF5_DIMS(rank == 0));
    hdf5_write_dataset(file, "cd", &cd, 1, HDF5_DIMS(rank == 0));
    hdf5_file_close(file);
}
