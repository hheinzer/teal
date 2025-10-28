#include "h5io.h"

#include <string.h>

#include "assert.h"
#include "sync.h"

hid_t h5io_file_create(const char *name)
{
    assert(name);

    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist, sync.comm, MPI_INFO_NULL);
    H5Pset_all_coll_metadata_ops(plist, true);
    H5Pset_coll_metadata_write(plist, true);

    hid_t file = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
    assert(file != H5I_INVALID_HID);

    H5Pclose(plist);
    return file;
}

hid_t h5io_file_open(const char *name)
{
    assert(name);

    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist, sync.comm, MPI_INFO_NULL);
    H5Pset_all_coll_metadata_ops(plist, true);

    hid_t file = H5Fopen(name, H5F_ACC_RDONLY, plist);
    assert(file != H5I_INVALID_HID);

    H5Pclose(plist);
    return file;
}

void h5io_file_close(hid_t file)
{
    H5Fclose(file);
}

hid_t h5io_group_create(const char *name, hid_t loc)
{
    assert(name);
    return H5Gcreate(loc, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

hid_t h5io_group_open(const char *name, hid_t loc)
{
    assert(name);
    return H5Gopen(loc, name, H5P_DEFAULT);
}

void h5io_group_close(hid_t group)
{
    H5Gclose(group);
}

void h5io_attribute_write(const char *object, const char *name, const void *buf, hsize_t dims,
                          hid_t type, hid_t loc)
{
    assert(name && buf && dims > 0);

    bool close_type = false;
    if (H5Tequal(type, H5T_C_S1) > 0) {
        type = H5Tcopy(H5T_C_S1);
        H5Tset_size(type, dims);
        H5Tset_strpad(type, H5T_STR_NULLPAD);
        close_type = true;
        dims = 1;
    }

    hid_t space = (dims == 1) ? H5Screate(H5S_SCALAR) : H5Screate_simple(1, &dims, 0);
    hid_t attr = object ? H5Aopen_by_name(loc, object, name, H5P_DEFAULT, H5P_DEFAULT)
                        : H5Acreate(loc, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, type, buf);

    H5Aclose(attr);
    H5Sclose(space);
    if (close_type) {
        H5Tclose(type);
    }
}

static bool attribute_dims_match(hid_t attr, hsize_t dims)
{
    hid_t space = H5Aget_space(attr);
    long ndims = H5Sget_simple_extent_ndims(space);

    if (ndims == 0) {
        H5Sclose(space);
        return dims == 1;
    }

    if (ndims != 1) {
        H5Sclose(space);
        return false;
    }

    hsize_t attr_dims[H5S_MAX_RANK];
    H5Sget_simple_extent_dims(space, attr_dims, 0);
    H5Sclose(space);

    return attr_dims[0] == dims;
}

void h5io_attribute_read(const char *object, const char *name, void *buf, hsize_t dims, hid_t type,
                         hid_t loc)
{
    assert(name && buf && dims > 0);

    bool close_type = false;
    if (H5Tequal(type, H5T_C_S1) > 0) {
        type = H5Tcopy(H5T_C_S1);
        H5Tset_size(type, dims);
        H5Tset_strpad(type, H5T_STR_NULLPAD);
        close_type = true;
        dims = 1;
    }

    hid_t attr = object ? H5Aopen_by_name(loc, object, name, H5P_DEFAULT, H5P_DEFAULT)
                        : H5Aopen(loc, name, H5P_DEFAULT);
    assert(attribute_dims_match(attr, dims));

    H5Aread(attr, type, buf);

    H5Aclose(attr);
    if (close_type) {
        H5Tclose(type);
    }
}

void h5io_dataset_write(const char *name, const void *buf, const hsize_t *count, long rank,
                        hid_t type, hid_t loc)
{
    assert(name && buf && count && (0 < rank && rank <= H5S_MAX_RANK));

    hsize_t dims[H5S_MAX_RANK];
    memcpy(dims, count, rank * sizeof(*count));
    MPI_Allreduce(count, dims, 1, HSIZE_AS_MPI_TYPE, MPI_SUM, sync.comm);

    hsize_t offset[H5S_MAX_RANK] = {0};
    MPI_Exscan(count, offset, 1, HSIZE_AS_MPI_TYPE, MPI_SUM, sync.comm);

    bool close_type = false;
    if (H5Tequal(type, H5T_C_S1) > 0) {
        assert(rank == 2);
        type = H5Tcopy(H5T_C_S1);
        H5Tset_size(type, count[1]);
        H5Tset_strpad(type, H5T_STR_NULLPAD);
        close_type = true;
        rank = 1;
    }

    hid_t file_space = H5Screate_simple(rank, dims, 0);
    hid_t dset = H5Dcreate(loc, name, type, file_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, 0, count, 0);

    hid_t plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);

    hid_t mem_space = H5Screate_simple(rank, count, 0);
    H5Dwrite(dset, type, mem_space, file_space, plist, buf);

    H5Sclose(mem_space);
    H5Pclose(plist);
    H5Dclose(dset);
    H5Sclose(file_space);
    if (close_type) {
        H5Tclose(type);
    }
}

static bool dataset_dims_match(hid_t dset, const hsize_t *count, long rank)
{
    hsize_t dims[H5S_MAX_RANK];
    memcpy(dims, count, rank * sizeof(*count));
    MPI_Allreduce(count, dims, 1, HSIZE_AS_MPI_TYPE, MPI_SUM, sync.comm);

    hsize_t dset_dims[H5S_MAX_RANK];
    hid_t space = H5Dget_space(dset);
    H5Sget_simple_extent_dims(space, dset_dims, 0);
    H5Sclose(space);

    for (long i = 0; i < rank; i++) {
        if (dims[i] != dset_dims[i]) {
            return false;
        }
    }
    return true;
}

void h5io_dataset_read(const char *name, void *buf, const hsize_t *count, long rank, hid_t type,
                       hid_t loc)
{
    assert(name && buf && count && (0 < rank && rank <= H5S_MAX_RANK));

    hsize_t offset[H5S_MAX_RANK] = {0};
    MPI_Exscan(count, offset, 1, HSIZE_AS_MPI_TYPE, MPI_SUM, sync.comm);

    bool close_type = false;
    if (H5Tequal(type, H5T_C_S1) > 0) {
        assert(rank == 2);
        type = H5Tcopy(H5T_C_S1);
        H5Tset_size(type, count[1]);
        H5Tset_strpad(type, H5T_STR_NULLPAD);
        close_type = true;
        rank = 1;
    }

    hid_t dset = H5Dopen(loc, name, H5P_DEFAULT);
    assert(dataset_dims_match(dset, count, rank));

    hid_t file_space = H5Dget_space(dset);
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, 0, count, 0);

    hid_t plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);

    hid_t mem_space = H5Screate_simple(rank, count, 0);
    H5Dread(dset, type, mem_space, file_space, plist, buf);

    H5Sclose(mem_space);
    H5Pclose(plist);
    H5Sclose(file_space);
    H5Dclose(dset);
    if (close_type) {
        H5Tclose(type);
    }
}

void h5io_link_create(const char *file, const char *object, const char *name, hid_t loc)
{
    assert(file && object && name);
    H5Lcreate_external(file, object, loc, name, H5P_DEFAULT, H5P_DEFAULT);
}
