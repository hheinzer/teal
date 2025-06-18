#pragma once

#include <hdf5.h>

#define H5IO_UCHAR H5T_NATIVE_UCHAR
#define H5IO_INT H5T_NATIVE_INT
#define H5IO_LONG H5T_NATIVE_LONG
#define H5IO_DOUBLE H5T_NATIVE_DOUBLE
#define H5IO_STRING H5T_C_S1

hid_t h5io_file_create(const char *name);
hid_t h5io_file_open(const char *name);
void h5io_file_close(hid_t file);

hid_t h5io_group_create(const char *name, hid_t loc);
hid_t h5io_group_open(const char *name, hid_t loc);
void h5io_group_close(hid_t group);

/* If `!object`, the attribute is attached to `loc`; otherwise it is attached to `loc/object`. */
void h5io_attribute_write(const char *object, const char *name, const void *buf, hsize_t dims,
                          hid_t type, hid_t loc);
void h5io_attribute_read(const char *object, const char *name, void *buf, hsize_t dims, hid_t type,
                         hid_t loc);

/* Dataset I/O assumes decomposition along the leading dimension only. */
void h5io_dataset_write(const char *name, const void *buf, const hsize_t *count, int rank,
                        hid_t type, hid_t loc);
void h5io_dataset_read(const char *name, void *buf, const hsize_t *count, int rank, hid_t type,
                       hid_t loc);

void h5io_link_create(const char *file, const char *object, const char *name, hid_t loc);
