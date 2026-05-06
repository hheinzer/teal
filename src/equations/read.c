#include <assert.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>

#include "equations.h"
#include "h5io.h"
#include "kdtree.h"
#include "sync.h"
#include "teal.h"

static void read_field_data(const Equations *eqns, double *time, hid_t loc)
{
    int root = (sync.rank == 0);

    int num = eqns->properties.num;
    String *name = eqns->properties.name;
    double *property = eqns->properties.data;

    hid_t group = h5io_group_open("FieldData", loc);
    h5io_dataset_read("time", time, root, 1, H5T_NATIVE_DOUBLE, group);
    for (int i = 0; i < num; i++) {
        h5io_dataset_read(name[i], &property[i], root, 1, H5T_NATIVE_DOUBLE, group);
    }
    h5io_group_close(group);

    MPI_Bcast(time, 1, MPI_DOUBLE, 0, sync.comm);
    MPI_Bcast(property, num, MPI_DOUBLE, 0, sync.comm);
}

typedef struct {
    long beg, end;
} Slice;

static Slice slice_rank(long total)
{
    assert(sync.rank >= 0 && sync.size >= 0);
    long rank = sync.rank;
    long size = sync.size;
    long base = total / size;
    long extra = total % size;
    long beg = (rank * base) + ((rank < extra) ? rank : extra);
    long end = beg + base + (rank < extra);
    return (Slice){beg, end};
}

static void read_variables(const EquationsVariables *variables, int num_cells, void *data_,
                           hid_t loc)
{
    if (variables->num == 0) {
        return;
    }

    int num = variables->num;
    int stride = variables->stride;
    String *name = variables->name;
    int *dimension = variables->dimension;
    double (*data)[stride] = data_;

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
        double (*buf)[dimension[i]] = buf_;
        h5io_dataset_read(name[i], buf, num_cells, dimension[i], H5T_NATIVE_DOUBLE, loc);
        for (int j = 0; j < num_cells; j++) {
            memcpy(&data[j][off], buf[j], sizeof(*buf));
        }
        off += dimension[i];
    }
    assert(off == stride);

    teal_free(buf_);
}

typedef struct {
    Vector center;
    double best_dist2;
    long global;
} Query;

static MPI_Datatype datatype_query(void)
{
    MPI_Datatype datatype;
    int len[3] = {3, 1, 1};
    MPI_Aint off[3] = {offsetof(Query, center), offsetof(Query, best_dist2),
                       offsetof(Query, global)};
    MPI_Datatype type[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_LONG};
    MPI_Type_create_struct(3, len, off, type, &datatype);
    return sync_resized(datatype, sizeof(Query));
}

static long *compute_indices(const Equations *eqns, const Vector *center, int num_local,
                             int num_cells)
{
    Kdtree *tree = kdtree_init(center, num_local);

    int cap = num_cells;
    sync_max(&cap, 1, MPI_INT);

    Query *query = teal_calloc(cap, sizeof(*query));
    for (int i = 0; i < num_cells; i++) {
        query[i].center = eqns->mesh->cells.center[i];
        query[i].best_dist2 = DBL_MAX;
    }

    long prefix = num_local;
    sync_prefix(&prefix, 1, MPI_LONG);

    MPI_Datatype type = datatype_query();
    int num = num_cells;
    for (int rank = 0; rank < sync.size; rank++) {
        for (int i = 0; i < num; i++) {
            int index = -1;
            kdtree_nearest(tree, query[i].center, &index, 1);
            assert(index >= 0);
            double dist2 = vector_norm2(vector_sub(query[i].center, center[index]));
            if (dist2 < query[i].best_dist2) {
                query[i].best_dist2 = dist2;
                query[i].global = prefix + index;
            }
        }
        sync_rotate(query, &num, cap, type, 1);
    }
    assert(num == num_cells);
    MPI_Type_free(&type);

    long *global = teal_calloc(num_cells, sizeof(*global));
    for (int i = 0; i < num_cells; i++) {
        global[i] = query[i].global;
    }

    teal_free(query);
    kdtree_deinit(tree);
    return global;
}

static void read_cell_data(const Equations *eqns, int num_cells, hid_t loc)
{
    int root = (sync.rank == 0);

    long num_total;
    h5io_dataset_read("NumberOfCells", &num_total, root, 1, H5T_NATIVE_LONG, loc);
    MPI_Bcast(&num_total, 1, MPI_LONG, 0, sync.comm);

    Slice slice = slice_rank(num_total);
    assert(slice.end - slice.beg <= INT_MAX);
    int num_local = (int)(slice.end - slice.beg);
    assert(num_local > 0);

    hid_t group = h5io_group_open("CellData", loc);

    Vector *center = teal_calloc(num_local, sizeof(*center));
    h5io_dataset_read("center", center, num_local, 3, H5T_NATIVE_DOUBLE, group);

    int stride = eqns->primitive.stride;
    double (*data)[stride] = teal_calloc(num_local, sizeof(*data));
    read_variables(&eqns->primitive, num_local, data, group);

    h5io_group_close(group);

    long *global = compute_indices(eqns, center, num_local, num_cells);

    double (*primitive)[stride] = eqns->primitive.data;
    sync_collect(data, primitive, global, num_local, num_cells, MPI_DOUBLE, stride);

    teal_free(global);
    teal_free(data);
    teal_free(center);
}

void equations_read(const Equations *eqns, const char *fname, double *time, int *index)
{
    assert(eqns && fname && time && index);

    char *underscore = strrchr(fname, '_');
    *index = underscore ? ((int)strtol(underscore + 1, 0, 10) + 1) : 0;

    hid_t file = h5io_file_open(fname);
    hid_t vtkhdf = h5io_group_open("VTKHDF", file);

    read_field_data(eqns, time, vtkhdf);

    int num_cells = eqns->mesh->cells.num_inner;
    read_cell_data(eqns, num_cells, vtkhdf);

    h5io_group_close(vtkhdf);
    h5io_file_close(file);
}
