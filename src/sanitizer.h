#pragma once

#ifndef __has_feature
#define __has_feature(x) 0
#endif

#if __has_feature(address_sanitizer) || defined(__SANITIZE_ADDRESS__)

#include <sanitizer/asan_interface.h>

#define MAKE_REGION_NOACCESS(addr, size) __asan_poison_memory_region(addr, size)
#define MAKE_REGION_ADDRESSABLE(addr, size) __asan_unpoison_memory_region(addr, size)

enum { REDZONE = 64 };

#else

#define MAKE_REGION_NOACCESS(addr, size) ((void)(addr), (void)(size))
#define MAKE_REGION_ADDRESSABLE(addr, size) ((void)(addr), (void)(size))

enum { REDZONE = 0 };

#endif
