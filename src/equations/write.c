#include <math.h>
#include <string.h>

#include "equations.h"
#include "teal/arena.h"
#include "teal/assert.h"
#include "teal/h5io.h"
#include "teal/sync.h"

static void write_field_data(const Equations *eqns, scalar time, hid_t loc)
{
    hid_t group = h5io_group_create("FieldData", loc);

    bool root = (sync.rank == 0);

    for (number i = 0; i < eqns->properties.num; i++) {
        h5io_dataset_write(eqns->properties.name[i], &eqns->properties.property[i], root, 1,
                           H5IO_SCALAR, group);
    }

    h5io_dataset_write("time", &time, root, 1, H5IO_SCALAR, group);

    h5io_group_close(group);
}

static void write_variables(const Name *name, const Type *type, const void *variable_, number num,
                            number stride, number num_cells, hid_t loc)
{
    const scalar(*variable)[stride] = variable_;
    for (number j = 0, i = 0; i < num; j += type[i++]) {
        Arena save = arena_save();

        scalar(*buf)[type[i]] = arena_malloc(num_cells, sizeof(*buf));
        for (number k = 0; k < num_cells; k++) {
            memcpy(buf[k], &variable[k][j], sizeof(*buf));
        }

        h5io_dataset_write(name[i], buf, num_cells, type[i], H5IO_SCALAR, loc);

        arena_load(save);
    }
}

static void write_user_variables(const Equations *eqns, scalar time, number num_cells, hid_t loc)
{
    Arena save = arena_save();

    vector *center = eqns->mesh->cells.center;
    scalar(*variable)[eqns->variables.stride] = eqns->variables.variable;
    scalar *property = eqns->properties.property;

    scalar(*user)[eqns->user.stride] = arena_calloc(num_cells, sizeof(*user));
    for (number i = 0; i < num_cells; i++) {
        eqns->user.compute(user[i], variable[i], property, center[i], time);
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
    hid_t group = h5io_group_create("CellData", loc);

    number num_cells = eqns->mesh->cells.off_periodic;
    write_variables(eqns->variables.name, eqns->variables.type, eqns->variables.variable,
                    eqns->variables.num, eqns->variables.stride, num_cells, group);

    if (eqns->user.num > 0) {
        write_user_variables(eqns, time, num_cells, group);
    }

    h5io_group_close(group);
}

void equations_write(const Equations *eqns, scalar time, const char *prefix, number index)
{
    assert(eqns && isfinite(time) && time >= 0 && prefix && index >= 0);

    char fname[128];
    sprintf(fname, "%s_%05td.vtkhdf", prefix, index);

    char lname[128];
    char *slash = strrchr(prefix, '/');
    sprintf(lname, "%s_mesh.h5", slash ? (slash + 1) : prefix);

    hid_t file = h5io_file_create(fname);
    hid_t vtkhdf = h5io_group_create("VTKHDF", file);

    h5io_attribute_write("Version", (number[]){1, 0}, 2, H5IO_NUMBER, vtkhdf);
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
