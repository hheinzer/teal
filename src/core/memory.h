#ifndef MEMORY_H
#define MEMORY_H

#define cleanup __attribute__((__cleanup__(memory_cleanup)))
#define fcleanup(func) __attribute__((__cleanup__(func)))

void *memory_calloc(long nmemb, long size);

void *memory_realloc(void *ptr, long nmemb, long size);

void memory_cleanup(void *ptr);

void *memory_duplicate(const void *ptr, long nmemb, long size);

char *memory_strdup(const char *src);

void *memory_copy(void *dest, const void *src, long nmemb, long size);

void *memory_setzero(void *ptr, long nmemb, long size);

#endif
