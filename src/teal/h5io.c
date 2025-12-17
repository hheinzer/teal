#include "h5io.h"

#include <assert.h>
#include <string.h>

#include "sync.h"
#include "utils.h"

hid_t h5io_file_open(const char *fname)
{
    assert(fname);

    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    assert(plist != H5I_INVALID_HID);
    H5Pset_fapl_mpio(plist, sync.comm, MPI_INFO_NULL);
    H5Pset_all_coll_metadata_ops(plist, true);

    hid_t file = H5Fopen(fname, H5F_ACC_RDONLY, plist);
    if (file == H5I_INVALID_HID) {
        error("could not open file (%s)", fname);
    }

    H5Pclose(plist);
    return file;
}

hid_t h5io_file_create(const char *fname)
{
    assert(fname);

    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    assert(plist != H5I_INVALID_HID);
    H5Pset_fapl_mpio(plist, sync.comm, MPI_INFO_NULL);
    H5Pset_all_coll_metadata_ops(plist, true);
    H5Pset_coll_metadata_write(plist, true);

    hid_t file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
    if (file == H5I_INVALID_HID) {
        error("could not create file (%s)", fname);
    }

    H5Pclose(plist);
    return file;
}

void h5io_file_close(hid_t file)
{
    assert(file != H5I_INVALID_HID);
    H5Fclose(file);
}

hid_t h5io_group_open(const char *name, hid_t loc)
{
    assert(name);
    hid_t group = H5Gopen(loc, name, H5P_DEFAULT);
    assert(group != H5I_INVALID_HID);
    return group;
}

hid_t h5io_group_create(const char *name, hid_t loc)
{
    assert(name);
    hid_t group = H5Gcreate(loc, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(group != H5I_INVALID_HID);
    return group;
}

void h5io_group_close(hid_t group)
{
    assert(group != H5I_INVALID_HID);
    H5Gclose(group);
}

long h5io_attribute_num(const char *name, hid_t loc)
{
    assert(name);

    hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
    assert(attr != H5I_INVALID_HID);

    hid_t space = H5Aget_space(attr);
    assert(space != H5I_INVALID_HID);

    long ndims = H5Sget_simple_extent_ndims(space);
    long num = 0;
    if (ndims == 0) {
        num = 1;
    }
    else if (ndims == 1) {
        hsize_t dims;
        H5Sget_simple_extent_dims(space, &dims, 0);
        num = dims;
    }

    H5Sclose(space);
    H5Aclose(attr);
    return num;
}

// Check whether an attribute has the expected length.
static bool attribute_num_matchs(hid_t attr, long num)
{
    bool match = false;

    hid_t space = H5Aget_space(attr);
    assert(space != H5I_INVALID_HID);

    long ndims = H5Sget_simple_extent_ndims(space);
    if (ndims == 0) {
        match = (num == 1);
    }
    else if (ndims != 1) {
        match = false;
    }
    else {
        hsize_t dims = num;
        hsize_t attr_dims;
        H5Sget_simple_extent_dims(space, &attr_dims, 0);
        match = (attr_dims == dims);
    }

    H5Sclose(space);
    return match;
}

void h5io_attribute_read(const char *name, void *buf, long num, hid_t type, hid_t loc)
{
    assert(name && buf && num > 0);

    bool close_type = false;
    if (H5Tequal(type, H5T_C_S1) > 0) {
        type = H5Tcopy(H5T_C_S1);
        assert(type != H5I_INVALID_HID);
        H5Tset_size(type, num);
        H5Tset_strpad(type, H5T_STR_NULLPAD);
        close_type = true;
        num = 1;
    }

    hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
    assert(attr != H5I_INVALID_HID);
    assert(attribute_num_matchs(attr, num));

    H5Aread(attr, type, buf);

    H5Aclose(attr);
    if (close_type) {
        H5Tclose(type);
    }
}

void h5io_attribute_write(const char *name, const void *buf, long num, hid_t type, hid_t loc)
{
    assert(name && buf && num > 0);

    bool close_type = false;
    if (H5Tequal(type, H5T_C_S1) > 0) {
        assert(num == 1);
        type = H5Tcopy(H5T_C_S1);
        assert(type != H5I_INVALID_HID);
        H5Tset_size(type, strlen(buf));
        H5Tset_strpad(type, H5T_STR_NULLPAD);
        close_type = true;
    }

    hsize_t dims = num;
    hid_t space = (num == 1) ? H5Screate(H5S_SCALAR) : H5Screate_simple(1, &dims, 0);
    assert(space != H5I_INVALID_HID);

    hid_t attr = H5Acreate(loc, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
    assert(attr != H5I_INVALID_HID);

    H5Awrite(attr, type, buf);

    H5Aclose(attr);
    H5Sclose(space);
    if (close_type) {
        H5Tclose(type);
    }
}

long h5io_dataset_num(const char *name, hid_t loc)
{
    assert(name);

    hid_t dset = H5Dopen(loc, name, H5P_DEFAULT);
    assert(dset != H5I_INVALID_HID);

    hid_t space = H5Dget_space(dset);
    assert(space != H5I_INVALID_HID);

    long ndims = H5Sget_simple_extent_ndims(space);
    long num = 0;
    if (ndims == 0) {
        num = 1;
    }
    else if (ndims == 1 || ndims == 2) {
        hsize_t dims[2];
        H5Sget_simple_extent_dims(space, dims, 0);
        num = dims[0];
    }

    H5Sclose(space);
    H5Dclose(dset);
    return num;
}

long h5io_dataset_len(const char *name, hid_t loc)
{
    assert(name);

    hid_t dset = H5Dopen(loc, name, H5P_DEFAULT);
    assert(dset != H5I_INVALID_HID);

    hid_t space = H5Dget_space(dset);
    assert(space != H5I_INVALID_HID);

    hid_t type = H5Dget_type(dset);
    assert(type != H5I_INVALID_HID);

    long len = 0;
    if (H5Tget_class(type) == H5T_STRING) {
        assert(H5Tis_variable_str(type) <= 0);
        len = H5Tget_size(type);
    }
    else {
        long ndims = H5Sget_simple_extent_ndims(space);
        if (ndims == 0 || ndims == 1) {
            len = 1;
        }
        else if (ndims == 2) {
            hsize_t dims[2];
            H5Sget_simple_extent_dims(space, dims, 0);
            len = dims[1];
        }
    }

    H5Tclose(type);
    H5Sclose(space);
    H5Dclose(dset);
    return len;
}

// Validate dataset shape against expected global rows and column length.
static bool dataset_dims_match(hid_t dset, long num, long len)
{
    bool match = false;

    hid_t space = H5Dget_space(dset);
    assert(space != H5I_INVALID_HID);

    long ndims = H5Sget_simple_extent_ndims(space);
    if (ndims <= 0 || 2 < ndims) {
        match = false;
    }
    else {
        hsize_t dims[2] = {sync_lsum(num), len};
        hsize_t dset_dims[2];
        H5Sget_simple_extent_dims(space, dset_dims, 0);
        if (ndims == 1) {
            match = (dset_dims[0] == dims[0]);
        }
        else {
            match = (dset_dims[0] == dims[0]) && (dset_dims[1] == dims[1]);
        }
    }

    H5Sclose(space);
    return match;
}

void h5io_dataset_read(const char *name, void *buf, long num, long len, hid_t type, hid_t loc)
{
    assert(name && buf && len > 0);

    long rank = (len == 1) ? 1 : 2;
    hsize_t count[2] = {num, len};
    hsize_t offset[2] = {sync_exsum(num), 0};

    bool close_type = false;
    if (H5Tequal(type, H5T_C_S1) > 0) {
        type = H5Tcopy(H5T_C_S1);
        assert(type != H5I_INVALID_HID);
        H5Tset_size(type, len);
        H5Tset_strpad(type, H5T_STR_NULLPAD);
        close_type = true;
        rank = 1;
    }

    hid_t dset = H5Dopen(loc, name, H5P_DEFAULT);
    assert(dset != H5I_INVALID_HID);
    assert(dataset_dims_match(dset, num, len));

    hid_t file_space = H5Dget_space(dset);
    assert(file_space != H5I_INVALID_HID);
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, 0, count, 0);

    hid_t plist = H5Pcreate(H5P_DATASET_XFER);
    assert(plist != H5I_INVALID_HID);
    H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);

    hid_t mem_space = H5Screate_simple(rank, count, 0);
    assert(mem_space != H5I_INVALID_HID);
    H5Dread(dset, type, mem_space, file_space, plist, buf);

    H5Sclose(mem_space);
    H5Pclose(plist);
    H5Sclose(file_space);
    H5Dclose(dset);
    if (close_type) {
        H5Tclose(type);
    }
}

void h5io_dataset_write(const char *name, const void *buf, long num, long len, hid_t type,
                        hid_t loc)
{
    assert(name && buf && len > 0);

    long rank = (len == 1) ? 1 : 2;
    hsize_t count[2] = {num, len};
    hsize_t dims[2] = {sync_lsum(num), len};
    hsize_t offset[2] = {sync_exsum(num), 0};

    bool close_type = false;
    if (H5Tequal(type, H5T_C_S1) > 0) {
        type = H5Tcopy(H5T_C_S1);
        assert(type != H5I_INVALID_HID);
        H5Tset_size(type, len);
        H5Tset_strpad(type, H5T_STR_NULLPAD);
        close_type = true;
        rank = 1;
    }

    hid_t file_space = H5Screate_simple(rank, dims, 0);
    assert(file_space != H5I_INVALID_HID);

    hid_t dset = H5Dcreate(loc, name, type, file_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(dset != H5I_INVALID_HID);
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, 0, count, 0);

    hid_t plist = H5Pcreate(H5P_DATASET_XFER);
    assert(plist != H5I_INVALID_HID);
    H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);

    hid_t mem_space = H5Screate_simple(rank, count, 0);
    assert(mem_space != H5I_INVALID_HID);
    H5Dwrite(dset, type, mem_space, file_space, plist, buf);

    H5Sclose(mem_space);
    H5Pclose(plist);
    H5Dclose(dset);
    H5Sclose(file_space);
    if (close_type) {
        H5Tclose(type);
    }
}

void h5io_link_create(const char *file, const char *object, const char *name, hid_t loc)
{
    assert(file && object && name);
    H5Lcreate_external(file, object, loc, name, H5P_DEFAULT, H5P_DEFAULT);
}
