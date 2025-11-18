#pragma once

#include <hdf5.h>
#include <stdint.h>

#include "teal.h"

#define H5IO_INT H5T_NATIVE_INT
#define H5IO_SCALAR (sizeof(scalar) == sizeof(float) ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE)
#define H5IO_STRING H5T_C_S1

hid_t h5io_file_open(const char *name);
hid_t h5io_file_create(const char *name);
void h5io_file_close(hid_t file);

hid_t h5io_group_open(const char *name, hid_t loc);
hid_t h5io_group_create(const char *name, hid_t loc);
void h5io_group_close(hid_t group);

int h5io_attribute_num(const char *name, hid_t loc);

void h5io_attribute_read(const char *name, void *buf, int num, hid_t type, hid_t loc);
void h5io_attribute_write(const char *name, const void *buf, int num, hid_t type, hid_t loc);

int h5io_dataset_num(const char *name, hid_t loc);
int h5io_dataset_len(const char *name, hid_t loc);

void h5io_dataset_read(const char *name, void *buf, int num, int len, hid_t type, hid_t loc);
void h5io_dataset_write(const char *name, const void *buf, int num, int len, hid_t type, hid_t loc);

void h5io_link_create(const char *file, const char *object, const char *name, hid_t loc);
