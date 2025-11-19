#include <assert.h>
#include <math.h>
#include <string.h>

#include "equations.h"
#include "teal/arena.h"
#include "teal/h5io.h"
#include "teal/sync.h"

static void write_field_data(const Equations *eqns, scalar time, hid_t loc)
{
    hid_t group = h5io_group_create("FieldData", loc);

    bool root = (sync.rank == 0);

    for (int i = 0; i < eqns->properties.num; i++) {
        h5io_dataset_write(eqns->properties.name[i], &eqns->properties.data[i], root, 1,
                           H5IO_SCALAR, group);
    }

    h5io_dataset_write("time", &time, root, 1, H5IO_SCALAR, group);

    h5io_group_close(group);
}

static void write_variables(const Name *name, const Type *type, const void *variable_, int num,
                            int stride, int num_cells, hid_t loc)
{
    const scalar(*variable)[stride] = variable_;
    for (int j = 0, i = 0; i < num; j += type[i++]) {
        Arena save = arena_save();

        scalar(*buf)[type[i]] = arena_malloc(num_cells, sizeof(*buf));
        for (int k = 0; k < num_cells; k++) {
            memcpy(buf[k], &variable[k][j], sizeof(*buf));
        }

        h5io_dataset_write(name[i], buf, num_cells, type[i], H5IO_SCALAR, loc);

        arena_load(save);
    }
}

static void write_user_variables(const Equations *eqns, scalar time, int num_cells, hid_t loc)
{
    Arena save = arena_save();

    vector *center = eqns->mesh->cells.center;
    scalar(*variable)[eqns->variables.stride] = eqns->variables.data;
    scalar *property = eqns->properties.data;

    scalar(*user)[eqns->user.stride] = arena_calloc(num_cells, sizeof(*user));
    for (int i = 0; i < num_cells; i++) {
        eqns->user.compute(user[i], property, center[i], time, variable[i]);
        if (eqns->user.conserved) {
            eqns->user.conserved(user[i], property);
        }
    }

    write_variables(eqns->user.name, eqns->user.type, user, eqns->user.num, eqns->user.stride,
                    num_cells, loc);

    arena_load(save);
}

static void write_cell_data(const Equations *eqns, scalar time, hid_t loc)
{
    Arena save = arena_save();

    hid_t group = h5io_group_create("CellData", loc);

    int num_cells = eqns->mesh->cells.off_periodic;
    equations_boundary(eqns, eqns->variables.data, time);
    write_variables(eqns->variables.name, eqns->variables.type, eqns->variables.data,
                    eqns->variables.num, eqns->variables.stride, num_cells, group);

    if (eqns->user.num > 0) {
        write_user_variables(eqns, time, num_cells, group);
    }

    scalar *step = arena_calloc(num_cells, sizeof(*step));
    equations_timestep(eqns, eqns->variables.data, step);
    h5io_dataset_write("step", step, num_cells, 1, H5IO_SCALAR, group);

    h5io_group_close(group);

    arena_load(save);
}

void equations_write(const Equations *eqns, scalar time, const char *prefix, int index)
{
    assert(eqns && isfinite(time) && time >= 0 && prefix && index >= 0);

    char fname[128];
    sprintf(fname, "%s_%05d.vtkhdf", prefix, index);

    char lname[128];
    char *slash = strrchr(prefix, '/');
    sprintf(lname, "%s_mesh.h5", slash ? (slash + 1) : prefix);

    hid_t file = h5io_file_create(fname);
    hid_t vtkhdf = h5io_group_create("VTKHDF", file);

    h5io_attribute_write("Version", (int[]){1, 0}, 2, H5IO_INT, vtkhdf);
    h5io_attribute_write("Type", "UnstructuredGrid", 1, H5IO_STRING, vtkhdf);

    h5io_link_create(lname, "/nodes/tot", "NumberOfPoints", vtkhdf);
    h5io_link_create(lname, "/nodes/coord", "Points", vtkhdf);

    h5io_link_create(lname, "/cells/tot", "NumberOfCells", vtkhdf);
    h5io_link_create(lname, "/cells/tot_idx", "NumberOfConnectivityIds", vtkhdf);
    h5io_link_create(lname, "/cells/node/off", "Offsets", vtkhdf);
    h5io_link_create(lname, "/cells/node/idx", "Connectivity", vtkhdf);
    h5io_link_create(lname, "/cells/type", "Types", vtkhdf);

    write_field_data(eqns, time, vtkhdf);
    write_cell_data(eqns, time, vtkhdf);

    h5io_group_close(vtkhdf);
    h5io_file_close(file);
}
