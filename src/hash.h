#ifndef HASH_H
#define HASH_H

#include <stdint.h>

uint64_t hash_fnv_64a(const void *ptr, const long nmemb, const long size);

#endif
