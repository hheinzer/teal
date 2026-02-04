#include <assert.h>

#include "h5io2.h"
#include "sync2.h"
#include "teal2.h"

static herr_t h5io_error(hid_t estack, void *ctx)
{
    (void)ctx;
    if (sync2.rank == 0) {
        H5Eprint2(estack, stderr);
        fflush(stderr);
    }
    teal2_error("h5io error");
}

hid_t h5io2_file_open(const char *fname)
{
    assert(fname);

    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist, sync2.comm, MPI_INFO_NULL);
    H5Pset_all_coll_metadata_ops(plist, 1);

    hid_t file = H5Fopen(fname, H5F_ACC_RDONLY, plist);
    if (file < 0) {
        teal2_error("could not open file (%s)", fname);
    }

    if (H5Eset_auto2(H5E_DEFAULT, h5io_error, 0) < 0) {
        teal2_error("could not install error handler");
    }

    H5Pclose(plist);
    return file;
}

hid_t h5io2_file_create(const char *fname)
{
    assert(fname);

    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist, sync2.comm, MPI_INFO_NULL);
    H5Pset_all_coll_metadata_ops(plist, 1);
    H5Pset_coll_metadata_write(plist, 1);

    hid_t file = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
    if (file < 0) {
        teal2_error("could not create file (%s)", fname);
    }

    if (H5Eset_auto2(H5E_DEFAULT, h5io_error, 0) < 0) {
        teal2_error("could not install error handler");
    }

    H5Pclose(plist);
    return file;
}

void h5io2_file_close(hid_t file)
{
    assert(file >= 0);
    H5Fclose(file);

    if (H5Eset_auto2(H5E_DEFAULT, 0, 0) < 0) {
        teal2_error("could not restore error handler");
    }
}

hid_t h5io2_group_open(const char *name, hid_t loc)
{
    assert(name && loc >= 0);
    return H5Gopen2(loc, name, H5P_DEFAULT);
}

hid_t h5io2_group_create(const char *name, hid_t loc)
{
    assert(name && loc >= 0);
    return H5Gcreate2(loc, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

void h5io2_group_close(hid_t group)
{
    assert(group >= 0);
    H5Gclose(group);
}

static hid_t type_string(size_t size)
{
    hid_t type = H5Tcopy(H5T_C_S1);
    H5Tset_size(type, size);
    H5Tset_strpad(type, H5T_STR_NULLPAD);
    return type;
}

static int attribute_dims_matchs(hid_t attr, hsize_t dims)
{
    int match = 0;
    hid_t space = H5Aget_space(attr);
    int ndims = H5Sget_simple_extent_ndims(space);
    if (ndims == 0) {
        match = (dims == 1);
    }
    else if (ndims == 1) {
        hsize_t attr_dims;
        H5Sget_simple_extent_dims(space, &attr_dims, 0);
        match = (dims == attr_dims);
    }
    H5Sclose(space);
    return match;
}

void h5io2_attribute_read(const char *name, void *buf, int num, hid_t type, hid_t loc)
{
    assert(name && buf && num > 0 && type >= 0 && loc >= 0);

    int close_type = 0;
    if (H5Tequal(type, H5T_C_S1) > 0) {
        type = type_string(num);
        num = 1;
        close_type = 1;
    }

    hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
    if (!attribute_dims_matchs(attr, num)) {
        teal2_error("attribute dims mismatch (%s)", name);
    }

    H5Aread(attr, type, buf);

    H5Aclose(attr);
    if (close_type) {
        H5Tclose(type);
    }
}

void h5io2_attribute_write(const char *name, const void *buf, int num, hid_t type, hid_t loc)
{
    assert(name && buf && num > 0 && type >= 0 && loc >= 0);

    int close_type = 0;
    if (H5Tequal(type, H5T_C_S1) > 0) {
        type = type_string(num);
        num = 1;
        close_type = 1;
    }

    hsize_t dims = num;
    hid_t space = (num == 1) ? H5Screate(H5S_SCALAR) : H5Screate_simple(1, &dims, 0);

    hid_t attr = H5Acreate2(loc, name, type, space, H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(attr, type, buf);

    H5Aclose(attr);
    H5Sclose(space);
    if (close_type) {
        H5Tclose(type);
    }
}

static int dataset_dims_match(hid_t dset, int rank, hsize_t num, hsize_t len)
{
    int match = 0;
    hid_t space = H5Dget_space(dset);
    int ndims = H5Sget_simple_extent_ndims(space);
    if (ndims == rank) {
        sync2_sum(&num, 1, HSIZE_AS_MPI_TYPE);
        hsize_t dset_dims[2];
        H5Sget_simple_extent_dims(space, dset_dims, 0);
        match = (num == dset_dims[0]);
        if (ndims == 2) {
            match &= (len == dset_dims[1]);
        }
    }
    H5Sclose(space);
    return match;
}

void h5io2_dataset_read(const char *name, void *buf, int num, int len, hid_t type, hid_t loc)
{
    assert(name && (buf || num == 0) && num >= 0 && len > 0 && type >= 0 && loc >= 0);

    int rank = (len == 1) ? 1 : 2;

    int close_type = 0;
    if (H5Tequal(type, H5T_C_S1) > 0) {
        type = type_string(len);
        rank = 1;
        close_type = 1;
    }

    hid_t dset = H5Dopen2(loc, name, H5P_DEFAULT);
    if (!dataset_dims_match(dset, rank, num, len)) {
        teal2_error("dataset dims mismatch (%s)", name);
    }

    hid_t fspace = H5Dget_space(dset);

    hsize_t count[2] = {num, len};
    hid_t mspace = H5Screate_simple(rank, count, 0);

    hsize_t start[2] = {num, 0};
    sync2_prefix(&start[0], 1, HSIZE_AS_MPI_TYPE);
    H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, 0, count, 0);

    hid_t plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);

    H5Dread(dset, type, mspace, fspace, plist, buf);

    H5Pclose(plist);
    H5Sclose(mspace);
    H5Sclose(fspace);
    H5Dclose(dset);
    if (close_type) {
        H5Tclose(type);
    }
}

void h5io2_dataset_write(const char *name, const void *buf, int num, int len, hid_t type, hid_t loc)
{
    assert(name && (buf || num == 0) && num >= 0 && len > 0 && type >= 0 && loc >= 0);

    int rank = (len == 1) ? 1 : 2;

    int close_type = 0;
    if (H5Tequal(type, H5T_C_S1) > 0) {
        type = type_string(len);
        rank = 1;
        close_type = 1;
    }

    hsize_t dims[2] = {num, len};
    sync2_sum(&dims[0], 1, HSIZE_AS_MPI_TYPE);
    hid_t fspace = H5Screate_simple(rank, dims, 0);

    hid_t dset = H5Dcreate2(loc, name, type, fspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t count[2] = {num, len};
    hid_t mspace = H5Screate_simple(rank, count, 0);

    hsize_t start[2] = {num, 0};
    sync2_prefix(&start[0], 1, HSIZE_AS_MPI_TYPE);
    H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, 0, count, 0);

    hid_t plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);

    H5Dwrite(dset, type, mspace, fspace, plist, buf);

    H5Pclose(plist);
    H5Sclose(mspace);
    H5Dclose(dset);
    H5Sclose(fspace);
    if (close_type) {
        H5Tclose(type);
    }
}

void h5io2_link_create(const char *file, const char *object, const char *name, hid_t loc)
{
    assert(file && object && name && loc >= 0);
    H5Lcreate_external(file, object, loc, name, H5P_DEFAULT, H5P_DEFAULT);
}
