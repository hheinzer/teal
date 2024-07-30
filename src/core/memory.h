#pragma once

/* Run 'func' on variable that is declared immediately after defer, when it exits scope. 'func'
 * must take a single argument of type pointer to variable. */
#define defer(func) [[gnu::cleanup(func)]]

/* Create a smart pointer, which is automatically freed once it exits scope. */
#define smart defer(memory_free)

/* Allocate memory for 'nmemb' elements of specified 'size'. */
void *memory_calloc(long nmemb, long size);

/* Reallocate memory of 'ptr' for 'nmemb' elements of specified 'size'. */
void *memory_realloc(void *ptr, long nmemb, long size);

/* Deallocate memory pointed to by 'ptr' and set 'ptr' to 0. WARNING: Even though the argument is
 * of type 'void *', the function must be called with '&ptr'. */
void memory_free(void *ptr);

/* Duplicate memory of 'ptr', which points to 'nmemb' elements of specified 'size'. */
void *memory_duplicate(const void *ptr, long nmemb, long size);

/* Duplicate string 'src'. */
char *memory_strdup(const char *src);

/* Copy memory to 'dest' from 'src', which points to 'nmemb' elements of specified 'size'. */
void *memory_copy(void *dest, const void *src, long nmemb, long size);

/* Zero memory of 'ptr', which points to 'nmemb' elements of specified 'size'. */
void *memory_setzero(void *ptr, long nmemb, long size);
