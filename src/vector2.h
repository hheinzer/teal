#pragma once

typedef struct {
    double x, y, z;
} Vector;

// Return `lhs + rhs`.
Vector vector2_add(Vector lhs, Vector rhs);

// Return `lhs - rhs`.
Vector vector2_sub(Vector lhs, Vector rhs);

// Return `lhs * rhs`.
Vector vector2_mul(double lhs, Vector rhs);

// Return `lhs / rhs`.
Vector vector2_div(Vector lhs, double rhs);

// Return the component-wise absolute value.
Vector vector2_abs(Vector vec);

// Add `rhs` to `lhs` in-place.
void vector2_inc(Vector *lhs, Vector rhs);

// Subtract `rhs` from `lhs` in-place.
void vector2_dec(Vector *lhs, Vector rhs);

// Scale `lhs` by `rhs` in-place.
void vector2_scale(Vector *lhs, double rhs);

// Return the dot product.
double vector2_dot(Vector lhs, Vector rhs);

// Return the 2-norm.
double vector2_norm(Vector vec);

// Return the squared 2-norm.
double vector2_norm2(Vector vec);

// Return the cross product.
Vector vector2_cross(Vector lhs, Vector rhs);
