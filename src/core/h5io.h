#ifndef H5IO_H
#define H5IO_H

#include <hdf5.h>

#define H5DIMS(...) ((hsize_t[]){__VA_ARGS__})
#define H5RANK(dims) (sizeof(dims) / sizeof(*dims))
#define H5TYPE(buf)                            \
    _Generic(*buf,                             \
        char: H5T_NATIVE_CHAR,                 \
        signed char: H5T_NATIVE_SCHAR,         \
        unsigned char: H5T_NATIVE_UCHAR,       \
        short: H5T_NATIVE_SHORT,               \
        unsigned short: H5T_NATIVE_USHORT,     \
        int: H5T_NATIVE_INT,                   \
        unsigned int: H5T_NATIVE_UINT,         \
        long: H5T_NATIVE_LONG,                 \
        unsigned long: H5T_NATIVE_ULONG,       \
        long long: H5T_NATIVE_LLONG,           \
        unsigned long long: H5T_NATIVE_ULLONG, \
        float: H5T_NATIVE_FLOAT,               \
        double: H5T_NATIVE_DOUBLE,             \
        long double: H5T_NATIVE_LDOUBLE)

hid_t h5_file_create(const char *name);

hid_t h5_file_open(const char *name);

void h5_file_close(hid_t file);

hid_t h5_group_create(const char *name, hid_t loc);

hid_t h5_group_open(const char *name, hid_t loc);

void h5_group_close(hid_t group);

#define h5_attribute_write(name, buf, dims, loc)                                              \
    _Generic(buf, char *: x__h5_attribute_write_string, default: x__h5_attribute_write_type)( \
        name, buf, dims, H5TYPE(buf), loc)
void x__h5_attribute_write_string(const char *name, const char *buf, hsize_t, hid_t, hid_t loc);
void x__h5_attribute_write_type(const char *name, const void *buf, hsize_t dims, hid_t type,
                                hid_t loc);

#define h5_attribute_read(name, buf, loc) x__h5_attribute_read(name, buf, loc)
void x__h5_attribute_read(const char *name, void *buf, hid_t loc);

#define h5_dataset_write(name, buf, dims, loc) \
    x__h5_dataset_write_type(name, buf, dims, H5RANK(dims), H5TYPE(buf), loc)
void x__h5_dataset_write_type(const char *name, const void *buf, const hsize_t *dims, int rank,
                              hid_t type, hid_t loc);

#define h5_dataset_read(name, buf, loc) x__h5_dataset_read_type(name, buf, loc)
void x__h5_dataset_read_type(const char *name, void *buf, hid_t loc);

void h5_link_create(const char *name, const char *file, const char *obj, hid_t loc);

#endif
