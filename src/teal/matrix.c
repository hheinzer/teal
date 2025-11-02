#include "matrix.h"

#include "assert.h"
#include "utils.h"
#include "vector.h"

matrix matrix_add(matrix lhs, matrix rhs)
{
    return (matrix){
        vector_add(lhs.x, rhs.x),
        vector_add(lhs.y, rhs.y),
        vector_add(lhs.z, rhs.z),
    };
}

matrix matrix_sub(matrix lhs, matrix rhs)
{
    return (matrix){
        vector_sub(lhs.x, rhs.x),
        vector_sub(lhs.y, rhs.y),
        vector_sub(lhs.z, rhs.z),
    };
}

matrix matrix_mul(scalar lhs, matrix rhs)
{
    return (matrix){
        vector_mul(lhs, rhs.x),
        vector_mul(lhs, rhs.y),
        vector_mul(lhs, rhs.z),
    };
}

matrix matrix_div(matrix lhs, scalar rhs)
{
    assert(!isclose(rhs, 0.0));
    return matrix_mul(1.0 / rhs, lhs);
}

matrix matrix_transpose(matrix mat)
{
    return (matrix){
        {mat.x.x, mat.y.x, mat.z.x},
        {mat.x.y, mat.y.y, mat.z.y},
        {mat.x.z, mat.y.z, mat.z.z},
    };
}

vector matrix_matvec(matrix lhs, vector rhs)
{
    return (vector){
        vector_dot(lhs.x, rhs),
        vector_dot(lhs.y, rhs),
        vector_dot(lhs.z, rhs),
    };
}

matrix matrix_matmul(matrix lhs, matrix rhs)
{
    matrix rhs_t = matrix_transpose(rhs);
    return (matrix){
        matrix_matvec(rhs_t, lhs.x),
        matrix_matvec(rhs_t, lhs.y),
        matrix_matvec(rhs_t, lhs.z),
    };
}

scalar matrix_det(matrix mat)
{
    return vector_dot(mat.x, vector_cross(mat.y, mat.z));
}

matrix matrix_inv(matrix mat)
{
    matrix cof = {
        vector_cross(mat.y, mat.z),
        vector_cross(mat.z, mat.x),
        vector_cross(mat.x, mat.y),
    };
    scalar det = vector_dot(mat.x, cof.x);
    return matrix_div(matrix_transpose(cof), det);
}
