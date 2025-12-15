#include <assert.h>
#include <math.h>
#include <string.h>

#include "equations.h"
#include "sync.h"
#include "teal/arena.h"
#include "teal/h5io.h"
#include "teal/sync.h"

// Write global field data (properties and time) to the VTKHDF file.
static void write_field_data(const Equations *eqns, scalar time, hid_t loc)
{
    hid_t group = h5io_group_create("FieldData", loc);

    bool root = (sync.rank == 0);

    h5io_dataset_write("time", &time, root, 1, H5IO_SCALAR, group);

    long num = eqns->properties.num;
    Name *name = eqns->properties.name;
    scalar *property = eqns->properties.data;
    for (long i = 0; i < num; i++) {
        h5io_dataset_write(name[i], &property[i], root, 1, H5IO_SCALAR, group);
    }

    h5io_group_close(group);
}

// Write a strided set of variables, splitting by component count.
static void write_variables(const long *dim, const Name *name, const void *variable_, long num,
                            long stride, long num_cells, hid_t loc)
{
    const scalar(*variable)[stride] = (void *)variable_;
    for (long j = 0, i = 0; i < num; j += dim[i++]) {
        Arena save = arena_save();

        scalar(*buf)[dim[i]] = arena_malloc(num_cells, sizeof(*buf));
        for (long k = 0; k < num_cells; k++) {
            memcpy(buf[k], &variable[k][j], sizeof(*buf));
        }
        h5io_dataset_write(name[i], buf, num_cells, dim[i], H5IO_SCALAR, loc);

        arena_load(save);
    }
}

// Compute and write user variables at the given time.
static void write_user_variables(const Equations *eqns, scalar time, long num_cells, hid_t loc)
{
    Arena save = arena_save();

    vector *center = eqns->mesh->cells.center;

    long stride_v = eqns->variables.stride;
    scalar(*variable)[stride_v] = eqns->variables.data;
    scalar *property = eqns->properties.data;

    long num = eqns->user.num;
    long stride_u = eqns->user.stride;
    long *dim = eqns->user.dim;
    Name *name = eqns->user.name;
    Compute *compute = eqns->user.compute;
    Update *conserved = eqns->user.conserved;

    scalar(*user)[stride_u] = arena_calloc(num_cells, sizeof(*user));

    if (conserved) {
        for (long i = 0; i < num_cells; i++) {
            compute(user[i], property, center[i], time, variable[i]);
            conserved(user[i], property);
        }
    }
    else {
        for (long i = 0; i < num_cells; i++) {
            compute(user[i], property, center[i], time, variable[i]);
        }
    }
    write_variables(dim, (void *)name, user, num, stride_u, num_cells, loc);

    arena_load(save);
}

// Write primary variables, time step, and optional user fields.
static void write_cell_data(const Equations *eqns, scalar time, hid_t loc)
{
    Arena save = arena_save();

    hid_t group = h5io_group_create("CellData", loc);

    long num_cells = eqns->mesh->cells.off_periodic;

    long num = eqns->variables.num;
    long stride = eqns->variables.stride;
    long *dim = eqns->variables.dim;
    Name *name = eqns->variables.name;
    scalar(*variable)[stride] = eqns->variables.data;

    Request req = sync_variables(eqns, variable, stride);
    equations_boundary(eqns, variable, time);
    sync_wait(eqns, req.send);
    write_variables(dim, (void *)name, variable, num, stride, num_cells, group);
    sync_wait(eqns, req.recv);

    if (eqns->user.num > 0) {
        write_user_variables(eqns, time, num_cells, group);
    }

    h5io_group_close(group);

    arena_load(save);
}

void equations_write(const Equations *eqns, const char *prefix, scalar time, long index)
{
    assert(eqns && prefix && isfinite(time) && time >= 0 && index >= 0);

    char fname[128];
    sprintf(fname, "%s_%05ld.vtkhdf", prefix, index);

    char lname[128];
    char *slash = strrchr(prefix, '/');
    sprintf(lname, "%s_mesh.vtkhdf", slash ? (slash + 1) : prefix);

    hid_t file = h5io_file_create(fname);
    hid_t vtkhdf = h5io_group_create("VTKHDF", file);

    h5io_attribute_write("Version", (long[]){1, 0}, 2, H5IO_LONG, vtkhdf);
    h5io_attribute_write("Type", "UnstructuredGrid", 1, H5IO_STRING, vtkhdf);

    h5io_link_create(lname, "/VTKHDF/NumberOfPoints", "NumberOfPoints", vtkhdf);
    h5io_link_create(lname, "/VTKHDF/Points", "Points", vtkhdf);

    h5io_link_create(lname, "/VTKHDF/NumberOfCells", "NumberOfCells", vtkhdf);
    h5io_link_create(lname, "/VTKHDF/NumberOfConnectivityIds", "NumberOfConnectivityIds", vtkhdf);
    h5io_link_create(lname, "/VTKHDF/Offsets", "Offsets", vtkhdf);
    h5io_link_create(lname, "/VTKHDF/Connectivity", "Connectivity", vtkhdf);
    h5io_link_create(lname, "/VTKHDF/Types", "Types", vtkhdf);

    write_field_data(eqns, time, vtkhdf);
    write_cell_data(eqns, time, vtkhdf);

    hid_t cell_data = h5io_group_open("CellData", vtkhdf);
    h5io_link_create(lname, "/VTKHDF/CellData/entity", "entity", cell_data);
    h5io_group_close(cell_data);

    h5io_group_close(vtkhdf);
    h5io_file_close(file);
}
