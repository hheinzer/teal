#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "equations.h"
#include "exchange.h"
#include "h5io.h"
#include "sync.h"
#include "teal.h"

static void link_mesh_data(const char *name, hid_t loc)
{
    String lname;
    char *slash = strrchr(name, '/');
    sprintf(lname, "%s_mesh.hdf", slash ? (slash + 1) : name);

    h5io_link_create(lname, "/VTKHDF/NumberOfPoints", "NumberOfPoints", loc);
    h5io_link_create(lname, "/VTKHDF/NumberOfCells", "NumberOfCells", loc);
    h5io_link_create(lname, "/VTKHDF/NumberOfConnectivityIds", "NumberOfConnectivityIds", loc);

    h5io_link_create(lname, "/VTKHDF/Points", "Points", loc);
    h5io_link_create(lname, "/VTKHDF/Offsets", "Offsets", loc);
    h5io_link_create(lname, "/VTKHDF/Connectivity", "Connectivity", loc);
    h5io_link_create(lname, "/VTKHDF/Types", "Types", loc);

    hid_t group = h5io_group_create("CellData", loc);
    h5io_link_create(lname, "/VTKHDF/CellData/vtkGhostType", "vtkGhostType", group);
    h5io_link_create(lname, "/VTKHDF/CellData/entity", "entity", group);
    h5io_link_create(lname, "/VTKHDF/CellData/rank", "rank", group);
    h5io_link_create(lname, "/VTKHDF/CellData/volume", "volume", group);
    h5io_link_create(lname, "/VTKHDF/CellData/center", "center", group);
    h5io_group_close(group);
}

static void write_field_data(const Equations *eqns, double time, hid_t loc)
{
    int root = (sync.rank == 0);

    int num = eqns->properties.num;
    String *name = eqns->properties.name;
    double *property = eqns->properties.data;

    hid_t group = h5io_group_create("FieldData", loc);
    h5io_dataset_write("time", &time, root, 1, H5T_NATIVE_DOUBLE, group);
    for (int i = 0; i < num; i++) {
        h5io_dataset_write(name[i], &property[i], root, 1, H5T_NATIVE_DOUBLE, group);
    }
    h5io_group_close(group);
}

static void write_variables(const EquationsVariables *variables, int num_cells, hid_t loc)
{
    if (variables->num == 0) {
        return;
    }

    int num = variables->num;
    int stride = variables->stride;
    String *name = variables->name;
    int *dimension = variables->dimension;
    double (*data)[stride] = variables->data;

    int cap = 0;
    for (int i = 0; i < num; i++) {
        if (dimension[i] > cap) {
            cap = dimension[i];
        }
    }
    assert(cap > 0);

    void *buf_ = teal_calloc(num_cells, sizeof(double[cap]));

    int off = 0;
    for (int i = 0; i < num; i++) {
        if (name[i][0]) {
            double (*buf)[dimension[i]] = buf_;
            for (int j = 0; j < num_cells; j++) {
                memcpy(buf[j], &data[j][off], sizeof(*buf));
            }
            h5io_dataset_write(name[i], buf, num_cells, dimension[i], H5T_NATIVE_DOUBLE, loc);
        }
        off += dimension[i];
    }
    assert(off == stride);

    teal_free(buf_);
}

static void compute_reference(const Equations *eqns, double time, int num_cells)
{
    Vector *center = eqns->mesh->cells.center;

    int stride = eqns->primitive.stride;
    double (*primitive)[stride] = eqns->primitive.data;
    double (*reference)[stride] = eqns->reference.data;
    double *property = eqns->properties.data;
    Compute *compute = eqns->reference.func.compute;

    for (int i = 0; i < num_cells; i++) {
        compute(reference[i], property, center[i], time, primitive[i]);
    }
}

static void write_cell_data(const Equations *eqns, double time, int num_cells, hid_t loc)
{
    hid_t group = h5io_group_open("CellData", loc);

    int stride = eqns->primitive.stride;
    double (*primitive)[stride] = eqns->primitive.data;

    Exchange exchange = equations_exchange(eqns, primitive, stride);
    equations_boundary(eqns, primitive, time);
    equations_exchange_wait_recv(eqns, exchange);
    write_variables(&eqns->primitive, num_cells, group);
    equations_exchange_wait_send(eqns, exchange);

    if (eqns->reference.func.compute) {
        compute_reference(eqns, time, num_cells);
        write_variables(&eqns->reference, num_cells, group);
    }

    h5io_group_close(group);
}

void equations_write(const Equations *eqns, const char *name, double time, int index)
{
    assert(eqns && name);

    String fname;
    sprintf(fname, "%s_%05d.hdf", name, index);

    hid_t file = h5io_file_create(fname);
    hid_t vtkhdf = h5io_group_create("VTKHDF", file);

    h5io_attribute_write("Version", (int[]){1, 0}, 2, H5T_NATIVE_INT, vtkhdf);
    h5io_attribute_write("Type", "UnstructuredGrid", 16, H5T_C_S1, vtkhdf);

    link_mesh_data(name, vtkhdf);
    write_field_data(eqns, time, vtkhdf);

    int num_cells = eqns->mesh->cells.off_periodic;
    write_cell_data(eqns, time, num_cells, vtkhdf);

    h5io_group_close(vtkhdf);
    h5io_file_close(file);
}
