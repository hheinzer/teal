#pragma once

typedef struct {
    double x, y, z;
} Vector;

typedef struct {
    Vector left, right;
} VectorPair;

typedef struct {
    Vector unit;
    double norm;
} VectorNorm;

// Return `lhs + rhs`.
Vector vector_add(Vector lhs, Vector rhs);

// Return `lhs - rhs`.
Vector vector_sub(Vector lhs, Vector rhs);

// Return `lhs * rhs`.
Vector vector_mul(double lhs, Vector rhs);

// Return `lhs / rhs`.
Vector vector_div(Vector lhs, double rhs);

// Compute `lhs += rhs` in-place.
void vector_iadd(Vector *lhs, Vector rhs);

// Compute `lhs -= rhs` in-place.
void vector_isub(Vector *lhs, Vector rhs);

// Compute `lhs *= rhs` in-place.
void vector_imul(Vector *lhs, double rhs);

// Compute `lhs /= rhs` in-place.
void vector_idiv(Vector *lhs, double rhs);

// Return the component-wise absolute value.
Vector vector_abs(Vector vec);

// Return the dot product.
double vector_dot(Vector lhs, Vector rhs);

// Return the 2-norm.
double vector_norm(Vector vec);

// Return the squared 2-norm.
double vector_norm2(Vector vec);

// Return the cross product.
Vector vector_cross(Vector lhs, Vector rhs);
