#pragma once

typedef int Compare(const void *, const void *);

// Return true if `lhs == rhs` with tolerance.
int isclose(double lhs, double rhs);

// Sort and remove duplicate array elements; returns new element count.
int unique(void *base, int num, int size, Compare compare);

// Binary search for the bucket index of `key` in a sorted array; returns `i` such that
// `base[i] <= key < base[i + 1]`, `-1` if `key < base[0]`, `num` if `key >= base[num]`,
// or `-2` if `num == 0`. The `base` array must have a minimum length of `num + 1`.
int digitize(const void *key, const void *base, int num, int size, Compare compare);
