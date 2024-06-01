#include "hdf5_io.h"

#include <hdf5.h>
#include <string.h>

hid_t hdf5_file_create(const char *name) {
    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL);
    H5Pset_all_coll_metadata_ops(plist, true);
    H5Pset_coll_metadata_write(plist, true);
    hid_t file = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist);
    H5Pclose(plist);
    return file;
}

hid_t hdf5_file_open(const char *name) {
    hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist, MPI_COMM_WORLD, MPI_INFO_NULL);
    H5Pset_all_coll_metadata_ops(plist, true);
    H5Pset_coll_metadata_write(plist, true);
    hid_t file = H5Fopen(name, H5F_ACC_RDONLY, H5P_DEFAULT);
    H5Pclose(plist);
    return file;
}

void hdf5_file_close(const hid_t file) {
    H5Fclose(file);
}

hid_t hdf5_group_create(const hid_t loc, const char *name) {
    return H5Gcreate(loc, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
}

hid_t hdf5_group_open(const hid_t loc, const char *name) {
    return H5Gopen(loc, name, H5P_DEFAULT);
}

void hdf5_group_close(const hid_t group) {
    H5Gclose(group);
}

void hdf5_write_attribute_str(const hid_t loc, const char *name, const char *buf) {
    hid_t type = H5Tcopy(H5T_C_S1);
    H5Tset_size(type, strlen(buf));
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate(loc, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, type, buf);
    H5Aclose(attr);
    H5Sclose(space);
    H5Tclose(type);
}

void hdf5_write_attribute_buf(const hid_t loc, const char *name, const void *buf, const int rank,
                              const hsize_t *dims, const hid_t type) {
    hid_t space = H5Screate_simple(rank, dims, H5P_DEFAULT);
    hid_t attr = H5Acreate(loc, name, type, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, type, buf);
    H5Aclose(attr);
    H5Sclose(space);
}

void hdf5_write_dataset_buf(const hid_t loc, const char *name, const void *buf, const int rank,
                            const hsize_t *count, const hid_t type) {
    hsize_t dims[rank] = {}, offset[rank] = {};
    for (long i = 0; i < rank; ++i) dims[i] = count[i];
    MPI_Allreduce(&count[0], &dims[0], 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Exscan(&count[0], &offset[0], 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

    hid_t filespace = H5Screate_simple(rank, dims, H5P_DEFAULT);
    hid_t dset = H5Dcreate(loc, name, type, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, 0, count, 0);
    hid_t memspace = H5Screate_simple(rank, count, H5P_DEFAULT);
    hid_t plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset, type, memspace, filespace, plist, buf);
    H5Pclose(plist);
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Dclose(dset);
}

void hdf5_link_create(const hid_t loc, const char *name, const char *file, const char *obj) {
    H5Lcreate_external(file, obj, loc, name, H5P_DEFAULT, H5P_DEFAULT);
}
