#pragma once

#define defer(f) [[gnu::cleanup(f)]]

#define smart defer(memory_free)

void memory_sum_setzero(void);

long memory_sum_get(void);

void *memory_calloc(long nmemb, long size);

void *memory_realloc(void *ptr, long nmemb, long size);

void *memory_copy(void *dest, const void *src, long nmemb, long size);

void *memory_duplicate(const void *ptr, long nmemb, long size);

void *memory_setzero(void *ptr, long nmemb, long size);

long memory_size(const void *ptr);

void memory_free(void *ptr);
