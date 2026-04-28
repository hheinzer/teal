#include "matrix.h"

#include <assert.h>

_Static_assert(sizeof(Matrix) == sizeof(Vector[3]), "Matrix must be packed");

Vector matrix_vector(Matrix mat, Vector vec)
{
    return (Vector){
        vector_dot(mat.x, vec),
        vector_dot(mat.y, vec),
        vector_dot(mat.z, vec),
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
    return vector_dot(mat.x, vector_cross(mat.y, mat.z));
}

Matrix matrix_inverse(Matrix mat)
{
    Matrix adj = {
        vector_cross(mat.y, mat.z),
        vector_cross(mat.z, mat.x),
        vector_cross(mat.x, mat.y),
    };

    double det = vector_dot(mat.x, adj.x);
    return (Matrix){
        vector_div(adj.x, det),
        vector_div(adj.y, det),
        vector_div(adj.z, det),
    };
}
