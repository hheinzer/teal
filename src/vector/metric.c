#include <math.h>

#include "vector.h"

double vector_dot(vector lhs, vector rhs)
{
    return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z);
}

double vector_norm(vector vec)
{
    return sqrt(vector_dot(vec, vec));
}

double vector_norm2(vector vec)
{
    return vector_dot(vec, vec);
}

vector vector_normalize(vector vec)
{
    return vector_div(vec, vector_norm(vec));
}

double vector_distance(vector lhs, vector rhs)
{
    return vector_norm(vector_sub(lhs, rhs));
}
