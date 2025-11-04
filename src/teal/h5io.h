#pragma once

#include <hdf5.h>

#define H5IO_UCHAR H5T_NATIVE_UCHAR
#define H5IO_INT H5T_NATIVE_INT
#define H5IO_LONG H5T_NATIVE_LONG
#define H5IO_STRING H5T_C_S1
#define H5IO_SCALAR (sizeof(scalar) == sizeof(float) ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE)

hid_t h5io_file_open(const char *name);
hid_t h5io_file_create(const char *name);
void h5io_file_close(hid_t file);

hid_t h5io_group_open(const char *name, hid_t loc);
hid_t h5io_group_create(const char *name, hid_t loc);
void h5io_group_close(hid_t group);

long h5io_attribute_num(const char *name, hid_t loc);

void h5io_attribute_read(const char *name, void *buf, long num, hid_t type, hid_t loc);
void h5io_attribute_write(const char *name, const void *buf, long num, hid_t type, hid_t loc);

long h5io_dataset_num(const char *name, hid_t loc);
long h5io_dataset_len(const char *name, hid_t loc);

void h5io_dataset_read(const char *name, void *buf, long num, long len, hid_t type, hid_t loc);
void h5io_dataset_write(const char *name, const void *buf, long num, long len, hid_t type,
                        hid_t loc);

void h5io_link_create(const char *file, const char *object, const char *name, hid_t loc);
