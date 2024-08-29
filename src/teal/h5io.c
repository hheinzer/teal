#include "h5io.h"

#include <string.h>

#include "sync.h"

hid_t h5_file_create(const char *name)
{
    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist, sync.comm, MPI_INFO_NULL);
    H5Pset_all_coll_metadata_ops(plist, true);
    H5Pset_coll_metadata_write(plist, true);
    hid_t file = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
    H5Pclose(plist);
    return file;
}

hid_t h5_file_open(const char *name)
{
    return H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);
}

void h5_file_close(hid_t file)
{
    H5Fclose(file);
}

hid_t h5_group_create(const char *name, hid_t loc)
{
    return H5Gcreate2(loc, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

hid_t h5_group_open(const char *name, hid_t loc)
{
    return H5Gopen2(loc, name, H5P_DEFAULT);
}

void h5_group_close(hid_t group)
{
    H5Gclose(group);
}

void x__h5_attribute_write_string(const char *name, const char *buf, hsize_t, hid_t, hid_t loc)
{
    hid_t type = H5Tcopy(H5T_C_S1);
    H5Tset_size(type, strlen(buf));
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate2(loc, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, type, buf);
    H5Aclose(attr);
    H5Sclose(space);
    H5Tclose(type);
}

void x__h5_attribute_write_type(const char *name, const void *buf, hsize_t dims, hid_t type,
                                hid_t loc)
{
    hid_t space = (dims == 1 ? H5Screate(H5S_SCALAR) : H5Screate_simple(1, &dims, 0));
    hid_t attr = H5Acreate2(loc, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, type, buf);
    H5Aclose(attr);
    H5Sclose(space);
}

void x__h5_attribute_read(const char *name, void *buf, hid_t loc)
{
    hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
    hid_t type = H5Aget_type(attr);
    H5Aread(attr, type, buf);
    H5Tclose(type);
    H5Aclose(attr);
}

void x__h5_dataset_write_type(const char *name, const void *buf, const hsize_t *dims, int rank,
                              hid_t type, hid_t loc)
{
    hsize_t fdims[rank], offset[rank];
    for (int i = 0; i < rank; ++i) fdims[i] = dims[i];
    for (int i = 0; i < rank; ++i) offset[i] = 0;
    MPI_Allreduce(&dims[0], &fdims[0], 1, MPI_UINT64_T, MPI_SUM, sync.comm);
    MPI_Exscan(&dims[0], &offset[0], 1, MPI_UINT64_T, MPI_SUM, sync.comm);
    hid_t mspace = H5Screate_simple(rank, dims, 0);
    hid_t fspace = H5Screate_simple(rank, fdims, 0);
    hid_t dset = H5Dcreate2(loc, name, type, fspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sselect_hyperslab(fspace, H5S_SELECT_SET, offset, 0, dims, 0);
    hid_t plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, H5FD_MPIO_INDEPENDENT);
    H5Dwrite(dset, type, mspace, fspace, plist, buf);
    H5Pclose(plist);
    H5Sclose(mspace);
    H5Sclose(fspace);
    H5Dclose(dset);
}

void x__h5_dataset_read_type(const char *name, void *buf, hid_t loc)
{
    hid_t dset = H5Dopen2(loc, name, H5P_DEFAULT);
    hid_t type = H5Dget_type(dset);
    H5Dread(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
    H5Tclose(type);
    H5Dclose(dset);
}

void h5_link_create(const char *name, const char *file, const char *obj, hid_t loc)
{
    H5Lcreate_external(file, obj, loc, name, H5P_DEFAULT, H5P_DEFAULT);
}
