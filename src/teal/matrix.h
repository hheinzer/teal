#pragma once

#include "teal.h"

matrix matrix_add(matrix lhs, matrix rhs);
matrix matrix_sub(matrix lhs, matrix rhs);
matrix matrix_mul(scalar lhs, matrix rhs);
matrix matrix_div(matrix lhs, scalar rhs);

void matrix_inc(matrix *lhs, matrix rhs);
void matrix_dec(matrix *lhs, matrix rhs);

matrix matrix_transpose(matrix mat);

vector matrix_matvec(matrix lhs, vector rhs);
matrix matrix_matmul(matrix lhs, matrix rhs);

scalar matrix_det(matrix mat);
matrix matrix_inv(matrix mat);
