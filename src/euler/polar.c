#include <assert.h>
#include <math.h>
#include <string.h>

#include "euler.h"
#include "teal/arena.h"
#include "teal/h5io.h"
#include "teal/sync.h"
#include "teal/utils.h"
#include "teal/vector.h"

static long find_entity(const MeshEntities *entities, const char *entity)
{
    for (long i = 0; i < entities->num; i++) {
        if (!strcmp(entities->name[i], entity)) {
            return i;
        }
    }
    error("invalid entity (%s)", entity);
}

void euler_polar(const Simulation *sim, const char *entity, const Euler *reference, scalar length,
                 scalar time)
{
    assert(sim && entity && reference && length > 0);
    Arena save = arena_save();

    const Mesh *mesh = sim->eqns->mesh;
    long *face_off = mesh->entities.face_off;

    const Equations *eqns = sim->eqns;
    Euler *variable = eqns->variables.data;

    long idx = find_entity(&mesh->entities, entity);
    long num = face_off[idx + 1] - face_off[idx];

    equations_boundary(eqns, eqns->variables.data, time);

    vector *center = arena_malloc(num, sizeof(*center));
    vector *normal = arena_malloc(num, sizeof(*normal));
    scalar *area = arena_malloc(num, sizeof(*area));
    scalar *pressure = arena_malloc(num, sizeof(*pressure));
    for (long j = 0, i = face_off[idx]; i < face_off[idx + 1]; i++, j++) {
        long left = mesh->faces.cell[i].left;
        long right = mesh->faces.cell[i].right;
        center[j] = mesh->faces.center[i];
        normal[j] = mesh->faces.basis[i].normal;
        area[j] = mesh->faces.area[i];
        pressure[j] = (variable[left].pressure + variable[right].pressure) / 2;
    }

    scalar *pressure_c = arena_malloc(num, sizeof(*pressure_c));
    vector coef = {0};
    scalar dynamic_pressure = reference->density * vector_norm2(reference->velocity) / 2;
    for (long i = 0; i < num; i++) {
        pressure_c[i] = (pressure[i] - reference->pressure) / dynamic_pressure;
        coef.x += (pressure[i] * area[i] * normal[i].x) / (dynamic_pressure * length);
        coef.y += (pressure[i] * area[i] * normal[i].y) / (dynamic_pressure * length);
    }
    coef = sync_vector_sum(coef);

    scalar alpha = atan2(reference->velocity.y, reference->velocity.x);
    scalar lift_c = (cos(alpha) * coef.y) - (sin(alpha) * coef.x);
    scalar drag_c = (cos(alpha) * coef.x) + (sin(alpha) * coef.y);

    println("Airfoil polar");
    println("\t lift coefficient : %g", lift_c);
    println("\t drag coefficient : %g", drag_c);

    if (sim->prefix) {
        char fname[128];
        sprintf(fname, "%s_polar.h5", sim->prefix);
        hid_t file = h5io_file_create(fname);
        bool root = (sync.rank == 0);
        h5io_dataset_write("center", center, num, 3, H5IO_SCALAR, file);
        h5io_dataset_write("pressure_c", pressure_c, num, 1, H5IO_SCALAR, file);
        h5io_dataset_write("lift_c", &lift_c, root, 1, H5IO_SCALAR, file);
        h5io_dataset_write("drag_c", &drag_c, root, 1, H5IO_SCALAR, file);
        h5io_file_close(file);
    }

    arena_load(save);
}
