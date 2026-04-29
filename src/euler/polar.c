#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "euler.h"
#include "h5io.h"
#include "mesh.h"
#include "sync.h"
#include "teal.h"

static int entity_index(const Mesh *mesh, const char *entity)
{
    for (int i = 0; i < mesh->entities.off_boundary; i++) {
        if (!strcmp(mesh->entities.name[i], entity)) {
            return i;
        }
    }
    teal_error("invalid entity (%s)", entity);
}

void euler_polar(const Simulation *sim, const char *entity, EulerPrimitive reference, double length,
                 double time)
{
    assert(sim && entity && length > 0);

    const Mesh *mesh = sim->eqns->mesh;
    int *face_off = mesh->entities.face_off;

    const Equations *eqns = sim->eqns;
    EulerPrimitive *primitive = eqns->primitive.data;

    int idx = entity_index(mesh, entity);
    int num = face_off[idx + 1] - face_off[idx];

    equations_boundary(eqns, primitive, time);

    Vector *center = teal_calloc(num, sizeof(*center));
    double *pressure_c = teal_calloc(num, sizeof(*pressure_c));
    Vector coef = {0};
    double dynamic_pressure = reference.density * vector_norm2(reference.velocity) / 2;
    for (int i = face_off[idx], j = 0; i < face_off[idx + 1]; i++, j++) {
        int left = mesh->faces.cell_idx[i].left;
        int right = mesh->faces.cell_idx[i].right;
        double area = mesh->faces.area[i];
        Vector normal = mesh->faces.basis[i].x;
        double pressure = (primitive[left].pressure + primitive[right].pressure) / 2;
        center[j] = mesh->faces.center[i];
        pressure_c[j] = (pressure - reference.pressure) / dynamic_pressure;
        coef.x += (pressure * area * normal.x) / (dynamic_pressure * length);
        coef.y += (pressure * area * normal.y) / (dynamic_pressure * length);
    }
    sync_sum(&coef, 3, MPI_DOUBLE);

    double alpha = atan2(reference.velocity.y, reference.velocity.x);
    double lift_c = (cos(alpha) * coef.y) - (sin(alpha) * coef.x);
    double drag_c = (cos(alpha) * coef.x) + (sin(alpha) * coef.y);

    teal_print("Airfoil polar");
    teal_print("\t lift coefficient : %g", lift_c);
    teal_print("\t drag coefficient : %g", drag_c);

    if (sim->prefix[0]) {
        double alpha_deg = alpha * 180 / acos(-1.0);

        char fname[256];
        sprintf(fname, "%s_polar_%g.h5", sim->prefix, alpha_deg);

        hid_t file = h5io_file_create(fname);

        int root = (sync.rank == 0);
        h5io_dataset_write("alpha", &alpha_deg, root, 1, H5T_NATIVE_DOUBLE, file);
        h5io_dataset_write("center", center, num, 3, H5T_NATIVE_DOUBLE, file);
        h5io_dataset_write("pressure_c", pressure_c, num, 1, H5T_NATIVE_DOUBLE, file);
        h5io_dataset_write("lift_c", &lift_c, root, 1, H5T_NATIVE_DOUBLE, file);
        h5io_dataset_write("drag_c", &drag_c, root, 1, H5T_NATIVE_DOUBLE, file);

        h5io_file_close(file);
    }

    teal_free(pressure_c);
    teal_free(center);
}
