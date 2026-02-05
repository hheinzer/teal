#include <assert.h>

#include "matrix2.h"
#include "utils2.h"

Vector matrix_vector(Matrix mat, Vector vec)
{
    return (Vector){
        vector2_dot(mat.x, vec),
        vector2_dot(mat.y, vec),
        vector2_dot(mat.z, vec),
    };
}

Matrix matrix_matrix(Matrix lhs, Matrix rhs)
{
    Matrix rhs_t = matrix_transpose(rhs);
    return (Matrix){
        matrix_vector(rhs_t, lhs.x),
        matrix_vector(rhs_t, lhs.y),
        matrix_vector(rhs_t, lhs.z),
    };
}

Matrix matrix_transpose(Matrix mat)
{
    return (Matrix){
        {mat.x.x, mat.y.x, mat.z.x},
        {mat.x.y, mat.y.y, mat.z.y},
        {mat.x.z, mat.y.z, mat.z.z},
    };
}

double matrix_determinant(Matrix mat)
{
    return vector2_dot(mat.x, vector2_cross(mat.y, mat.z));
}

Matrix matrix_inverse(Matrix mat)
{
    Matrix adj = {
        vector2_cross(mat.y, mat.z),
        vector2_cross(mat.z, mat.x),
        vector2_cross(mat.x, mat.y),
    };

    double det = vector2_dot(mat.x, adj.x);
    assert(!isclose(det, 0));

    double inv_det = 1 / det;
    return (Matrix){
        vector2_mul(inv_det, adj.x),
        vector2_mul(inv_det, adj.y),
        vector2_mul(inv_det, adj.z),
    };
}
