#ifndef UTILS_H
#define UTILS_H

#define MIN(a, b) (a < b ? a : b)
#define MAX(a, b) (a > b ? a : b)

void *utils_memdup(const void *src, const long nmemb, const long size);

char *utils_strdup(const char *src);

const char *utils_extension(const char *fname);

const char *utils_basename(const char *fname);

#endif
