#pragma once

typedef struct {
    double x, y, z;
} vector;

// Return the sum of two vectors.
vector vector_add(vector lhs, vector rhs);

// Return the difference of two vectors.
vector vector_sub(vector lhs, vector rhs);

// Return a vector multiplied by a scalar.
vector vector_mul(double val, vector vec);

// Return a vector divided by a scalar.
vector vector_div(vector vec, double val);

// Return the componentwise absolute value of a vector.
vector vector_abs(vector vec);

// Add rhs to lhs in place.
void vector_inc(vector *lhs, vector rhs);

// Subtract rhs from lhs in place.
void vector_dec(vector *lhs, vector rhs);

// Scale a vector in place by a scalar.
void vector_scale(vector *vec, double val);

// Return the dot product of two vectors.
double vector_dot(vector lhs, vector rhs);

// Return the Euclidean norm of a vector.
double vector_norm(vector vec);

// Return the squared Euclidean norm of a vector.
double vector_norm2(vector vec);

// Return the cross product of two vectors.
vector vector_cross(vector lhs, vector rhs);

// Return a normalized vector.
vector vector_normalize(vector vec);

// Return the distance between two vectors.
double vector_distance(vector lhs, vector rhs);

// Return the componentwise minimum of an array of vectors.
vector vector_min(const vector *arr, int num);

// Return the componentwise maximum of an array of vectors.
vector vector_max(const vector *arr, int num);

// Return the componentwise sum of an array of vectors.
vector vector_sum(const vector *arr, int num);
