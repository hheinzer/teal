#pragma once

#include <stdint.h>

uint64_t hash_fnv_64a(const void *ptr, long nmemb, long size);
