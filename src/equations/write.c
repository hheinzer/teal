#include <stdio.h>
#include <string.h>

#include "core/h5io.h"
#include "core/memory.h"
#include "core/utils.h"
#include "equations.h"
#include "teal.h"

void equations_write(const Equations *eqns, const char *prefix, long count, double time, long iter)
{
    // https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html#vtkhdf-file-format
    char fname[128], lname[128];
    sprintf(fname, "%s_%05ld.hdf", prefix, count);
    const char *slash = strrchr(prefix, '/');
    sprintf(lname, "%s_mesh.hdf", (slash ? slash + 1 : prefix));

    const int root = (teal.rank == 0);
    hid_t file = h5_file_create(fname);
    hid_t vtkhdf = h5_group_create("VTKHDF", file);

    const long version[] = {2, 2};
    h5_attribute_write("Version", version, 2, vtkhdf);
    h5_attribute_write("Type", "UnstructuredGrid", 0, vtkhdf);

    h5_link_create("NumberOfPoints", lname, "/VTKHDF/NumberOfPoints", vtkhdf);
    h5_link_create("NumberOfConnectivityIds", lname, "/VTKHDF/NumberOfConnectivityIds", vtkhdf);
    h5_link_create("NumberOfCells", lname, "/VTKHDF/NumberOfCells", vtkhdf);
    h5_link_create("Points", lname, "/VTKHDF/Points", vtkhdf);
    h5_link_create("Connectivity", lname, "/VTKHDF/Connectivity", vtkhdf);
    h5_link_create("Offsets", lname, "/VTKHDF/Offsets", vtkhdf);
    h5_link_create("Types", lname, "/VTKHDF/Types", vtkhdf);

    hid_t field = h5_group_create("FieldData", vtkhdf);
    h5_dataset_write("time", &time, H5DIMS(root), field);
    h5_dataset_write("iter", &iter, H5DIMS(root), field);
    h5_group_close(field);

    hid_t cell = h5_group_create("CellData", vtkhdf);
    h5_link_create("rank", lname, "/VTKHDF/CellData/rank", cell);
    h5_link_create("index", lname, "/VTKHDF/CellData/index", cell);
    h5_link_create("volume", lname, "/VTKHDF/CellData/volume", cell);
    h5_link_create("center", lname, "/VTKHDF/CellData/center", cell);
    h5_link_create("projection", lname, "/VTKHDF/CellData/projection", cell);

    const long n_inner_cells = eqns->mesh->n_inner_cells;
    smart double *buf = memory_calloc(n_inner_cells * N_DIMS, sizeof(*buf));

    const long n_vars = eqns->n_vars;
    const ALIAS(vdim, eqns->vars.dim);
    const ALIAS(vname, eqns->vars.name);
    const double(*vars)[n_vars] = (void *)eqns->vars.u;
    for (long n = 0, v = 0; v < n_vars; v += vdim[n], ++n) {
        for (long i = 0; i < n_inner_cells; ++i)
            for (long d = 0; d < vdim[n]; ++d) buf[i * vdim[n] + d] = vars[i][v + d];
        h5_dataset_write(vname[v], buf, H5DIMS(n_inner_cells, vdim[n]), cell);
    }

    const long n_user = eqns->n_user;
    const ALIAS(udim, eqns->user.dim);
    const ALIAS(uname, eqns->user.name);
    const ALIAS(x, eqns->mesh->cell.center);
    ALIAS(func, eqns->user.func);
    double user[n_user];
    for (long n = 0, v = 0; v < n_user; v += udim[n], ++n) {
        for (long i = 0; i < n_inner_cells; ++i) {
            func(user, vars[i], x[i], time);
            for (long d = 0; d < udim[n]; ++d) buf[i * udim[n] + d] = user[v + d];
        }
        h5_dataset_write(uname[v], buf, H5DIMS(n_inner_cells, udim[n]), cell);
    }

    h5_group_close(cell);
    h5_group_close(vtkhdf);
    h5_file_close(file);
}
