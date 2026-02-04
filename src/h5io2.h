#pragma once

#include <hdf5.h>

// Open an HDF5 file in read-only mode.
hid_t h5io2_file_open(const char *fname);

// Create or truncate an HDF5 file.
hid_t h5io2_file_create(const char *fname);

// Close an HDF5 file handle.
void h5io2_file_close(hid_t file);

// Open a group at `loc`.
hid_t h5io2_group_open(const char *name, hid_t loc);

// Create a group at `loc`.
hid_t h5io2_group_create(const char *name, hid_t loc);

// Close a group handle.
void h5io2_group_close(hid_t group);

// Read an attribute into `buf` (string attributes use `num` as buffer size).
void h5io2_attribute_read(const char *name, void *buf, int num, hid_t type, hid_t loc);

// Write an attribute from `buf` (string attributes use `num` as buffer size).
void h5io2_attribute_write(const char *name, const void *buf, int num, hid_t type, hid_t loc);

// Read a dataset slice of shape `[num, len]` (or `[num]` if `len == 1`).
void h5io2_dataset_read(const char *name, void *buf, int num, int len, hid_t type, hid_t loc);

// Write a dataset slice of shape `[num, len]` (or `[num]` if `len == 1`).
void h5io2_dataset_write(const char *name, const void *buf, int num, int len, hid_t type,
                         hid_t loc);

// Create an external link in `loc` named `name` that points to `file:object`.
void h5io2_link_create(const char *file, const char *object, const char *name, hid_t loc);
