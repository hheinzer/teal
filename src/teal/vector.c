#include "vector.h"

#include <math.h>

vector vector_min(vector lhs, vector rhs)
{
    return (vector){fmin(lhs.x, rhs.x), fmin(lhs.y, rhs.y), fmin(lhs.z, rhs.z)};
}

vector vector_max(vector lhs, vector rhs)
{
    return (vector){fmax(lhs.x, rhs.x), fmax(lhs.y, rhs.y), fmax(lhs.z, rhs.z)};
}

vector vector_add(vector lhs, vector rhs)
{
    return (vector){lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

vector vector_sub(vector lhs, vector rhs)
{
    return (vector){lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

vector vector_mul(vector lhs, scalar rhs)
{
    return (vector){lhs.x * rhs, lhs.y * rhs, lhs.z * rhs};
}

vector vector_div(vector lhs, scalar rhs)
{
    return vector_mul(lhs, 1.0 / rhs);
}

vector vector_abs(vector vec)
{
    return (vector){fabs(vec.x), fabs(vec.y), fabs(vec.z)};
}

scalar vector_dot(vector lhs, vector rhs)
{
    return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z);
}

scalar vector_len(vector vec)
{
    return sqrt(vector_dot(vec, vec));
}

vector vector_unit(vector vec)
{
    return vector_div(vec, vector_len(vec));
}

vector vector_cross(vector lhs, vector rhs)
{
    return (vector){
        (lhs.y * rhs.z) - (lhs.z * rhs.y),
        (lhs.z * rhs.x) - (lhs.x * rhs.z),
        (lhs.x * rhs.y) - (lhs.y * rhs.x),
    };
}
