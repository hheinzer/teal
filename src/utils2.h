#pragma once

typedef int Compare(const void *, const void *);

// Return true if `lhs == rhs` with tolerance.
int isclose(double lhs, double rhs);

// Copy `num * size` bytes from `src` to `dst`.
void *copy(void *dst, const void *src, int num, int size);

// Sort array elements in-place using compare.
void sort(void *base, int num, int size, Compare compare);

// Sort and remove duplicate array elements; returns new element count.
int unique(void *base, int num, int size, Compare compare);

// Linear search for `key`; returns pointer to match or `0`.
void *find(const void *key, const void *base, int num, int size, Compare compare);

// Binary search for `key` in a sorted array; returns pointer to match or `0`.
void *search(const void *key, const void *base, int num, int size, Compare compare);

// Binary search for the bucket index of `key` in a sorted array; returns `i` such that
// `base[i] <= key < base[i + 1]`, `-1` if `key < base[0]`, `num` if `key >= base[num]`,
// or `-2` if `num == 0`. The `base` array must have a minimum length of `num + 1`.
int digitize(const void *key, const void *base, int num, int size, Compare compare);
