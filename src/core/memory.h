#ifndef MEMORY_H
#define MEMORY_H

#define defer(func) [[gnu::cleanup(func)]]
#define smart defer(memory_free)

void *memory_calloc(long nmemb, long size);

void *memory_realloc(void *ptr, long nmemb, long size);

/* Deallocates memory pointed to by 'ptr' and sets 'ptr' to 0. WARNING: Even though the argument is
 * of type 'void *', the function must be called with '&ptr'. */
void memory_free(void *ptr);

void *memory_duplicate(const void *ptr, long nmemb, long size);

char *memory_strdup(const char *src);

void *memory_copy(void *dest, const void *src, long nmemb, long size);

void *memory_setzero(void *ptr, long nmemb, long size);

#endif
