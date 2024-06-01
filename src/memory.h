#ifndef MEMORY_H
#define MEMORY_H

#define cleanup __attribute__((__cleanup__(memory_cleanup)))
#define fcleanup(func) __attribute__((__cleanup__(func)))

void *memory_calloc(const long nmemb, const long size);

void *memory_realloc(void *ptr, const long nmemb, const long size);

void memory_cleanup(void *ptr);

#endif
