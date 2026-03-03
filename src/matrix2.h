#pragma once

#include "vector2.h"

typedef struct {
    Vector x, y, z;  // row-major order
} Matrix;

// Return matrix-Vector multiplication.
Vector matrix_vector(Matrix mat, Vector vec);

// Return matrix-matrix multiplication.
Matrix matrix_matrix(Matrix lhs, Matrix rhs);

// Return the transpose.
Matrix matrix_transpose(Matrix mat);

// Return the determinant.
double matrix_determinant(Matrix mat);

// Return the inverse.
Matrix matrix_inverse(Matrix mat);
