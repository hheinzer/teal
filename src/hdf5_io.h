#ifndef HDF5_IO_H
#define HDF5_IO_H

#include <hdf5.h>

#define HDF5_DIMS(...) \
    (hsize_t[]) { __VA_ARGS__ }

#define HDF5_TYPE(buf)                           \
    _Generic(buf,                                \
        unsigned char *: H5T_NATIVE_UCHAR,       \
        long *: H5T_NATIVE_LONG,                 \
        hsize_t *: H5T_NATIVE_HSIZE,             \
        double *: H5T_NATIVE_DOUBLE,             \
        const unsigned char *: H5T_NATIVE_UCHAR, \
        const long *: H5T_NATIVE_LONG,           \
        const hsize_t *: H5T_NATIVE_HSIZE,       \
        const double *: H5T_NATIVE_DOUBLE)

hid_t hdf5_file_create(const char *name);
hid_t hdf5_file_open(const char *name);
void hdf5_file_close(const hid_t file);

hid_t hdf5_group_create(const hid_t loc, const char *name);
hid_t hdf5_group_open(const hid_t loc, const char *name);
void hdf5_group_close(const hid_t group);

void hdf5_write_attribute_str(const hid_t loc, const char *name, const char *buf);
void hdf5_write_attribute_buf(const hid_t loc, const char *name, const void *buf, const int rank,
                              const hsize_t *dims, const hid_t type);
#define hdf5_write_attribute(loc, name, buf, ...)          \
    _Generic(buf,                                          \
        char *: hdf5_write_attribute_str,                  \
        long *: hdf5_write_attribute_buf,                  \
        const char *: hdf5_write_attribute_str,            \
        const long *: hdf5_write_attribute_buf)(loc, name, \
                                                buf __VA_OPT__(, __VA_ARGS__, HDF5_TYPE(buf)))

void hdf5_write_dataset_buf(const hid_t loc, const char *name, const void *buf, const int rank,
                            const hsize_t *count, const hid_t type);
#define hdf5_write_dataset(loc, name, buf, rank, count) \
    hdf5_write_dataset_buf(loc, name, buf, rank, count, HDF5_TYPE(buf))

void hdf5_link_create(const hid_t loc, const char *name, const char *file, const char *obj);

#endif
