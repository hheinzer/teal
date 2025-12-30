#include <math.h>

#include "vector.h"

vector vector_add(vector lhs, vector rhs)
{
    return (vector){
        .x = lhs.x + rhs.x,
        .y = lhs.y + rhs.y,
        .z = lhs.z + rhs.z,
    };
}

vector vector_sub(vector lhs, vector rhs)
{
    return (vector){
        .x = lhs.x - rhs.x,
        .y = lhs.y - rhs.y,
        .z = lhs.z - rhs.z,
    };
}

vector vector_mul(double val, vector vec)
{
    return (vector){
        .x = vec.x * val,
        .y = vec.y * val,
        .z = vec.z * val,
    };
}

vector vector_div(vector vec, double val)
{
    return (vector){
        .x = vec.x / val,
        .y = vec.y / val,
        .z = vec.z / val,
    };
}

vector vector_abs(vector vec)
{
    return (vector){
        .x = fabs(vec.x),
        .y = fabs(vec.y),
        .z = fabs(vec.z),
    };
}

void vector_inc(vector *lhs, vector rhs)
{
    lhs->x += rhs.x;
    lhs->y += rhs.y;
    lhs->z += rhs.z;
}

void vector_dec(vector *lhs, vector rhs)
{
    lhs->x -= rhs.x;
    lhs->y -= rhs.y;
    lhs->z -= rhs.z;
}

void vector_scale(vector *vec, double val)
{
    vec->x *= val;
    vec->y *= val;
    vec->z *= val;
}
