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

void h5io_attribute_write(const char *name, const void *buf, long num, hid_t type, hid_t loc)
{
    assert(name && buf && num > 0);

    bool close_type = false;
    if (H5Tequal(type, H5T_C_S1) > 0) {
        assert(num == 1);
        type = H5Tcopy(H5T_C_S1);
        H5Tset_size(type, strlen(buf));
        H5Tset_strpad(type, H5T_STR_NULLPAD);
        close_type = true;
    }

    hsize_t dims = num;
    hid_t space = (num == 1) ? H5Screate(H5S_SCALAR) : H5Screate_simple(1, &dims, 0);
    hid_t attr = H5Acreate(loc, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, type, buf);

    H5Aclose(attr);
    H5Sclose(space);
    if (close_type) {
        H5Tclose(type);
    }
}

static bool attribute_num_matchs(hid_t attr, long num)
{
    hid_t space = H5Aget_space(attr);
    int ndims = H5Sget_simple_extent_ndims(space);

    if (ndims == 0) {
        H5Sclose(space);
        return num == 1;
    }

    if (ndims != 1) {
        H5Sclose(space);
        return false;
    }

    hsize_t attr_dims;
    H5Sget_simple_extent_dims(space, &attr_dims, 0);
    H5Sclose(space);

    hsize_t dims = num;
    return attr_dims == dims;
}

void h5io_attribute_read(const char *name, void *buf, long num, hid_t type, hid_t loc)
{
    assert(name && buf && num > 0);

    bool close_type = false;
    if (H5Tequal(type, H5T_C_S1) > 0) {
        type = H5Tcopy(H5T_C_S1);
        H5Tset_size(type, num);
        H5Tset_strpad(type, H5T_STR_NULLPAD);
        close_type = true;
        num = 1;
    }

    hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
    assert(attribute_num_matchs(attr, num));

    H5Aread(attr, type, buf);

    H5Aclose(attr);
    if (close_type) {
        H5Tclose(type);
    }
}

void h5io_dataset_write(const char *name, const void *buf, long num, long len, hid_t type,
                        hid_t loc)
{
    assert(name && buf && len > 0);

    int rank = (len == 1) ? 1 : 2;
    hsize_t count[2] = {num, len};
    hsize_t dims[2] = {sync_lsum(num), len};
    hsize_t offset[2] = {sync_lexsum(num), 0};

    bool close_type = false;
    if (H5Tequal(type, H5T_C_S1) > 0) {
        type = H5Tcopy(H5T_C_S1);
        H5Tset_size(type, len);
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

static bool dataset_dims_match(hid_t dset, long num, long len)
{
    hid_t space = H5Dget_space(dset);
    int ndims = H5Sget_simple_extent_ndims(space);

    if (ndims <= 0 || 2 < ndims) {
        H5Sclose(space);
        return false;
    }

    hsize_t dims[2] = {sync_lsum(num), len};

    hsize_t dset_dims[2];
    H5Sget_simple_extent_dims(space, dset_dims, 0);
    H5Sclose(space);

    for (long i = 0; i < ndims; i++) {
        if (dset_dims[i] != dims[i]) {
            return false;
        }
    }
    return true;
}

void h5io_dataset_read(const char *name, void *buf, long num, long len, hid_t type, hid_t loc)
{
    assert(name && buf && len > 0);

    int rank = (len == 1) ? 1 : 2;
    hsize_t count[2] = {num, len};
    hsize_t offset[2] = {sync_lexsum(num), 0};

    bool close_type = false;
    if (H5Tequal(type, H5T_C_S1) > 0) {
        type = H5Tcopy(H5T_C_S1);
        H5Tset_size(type, len);
        H5Tset_strpad(type, H5T_STR_NULLPAD);
        close_type = true;
        rank = 1;
    }

    hid_t dset = H5Dopen(loc, name, H5P_DEFAULT);
    assert(dataset_dims_match(dset, num, len));

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
