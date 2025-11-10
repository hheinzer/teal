#include "vector.h"

#include <math.h>

#include "assert.h"
#include "utils.h"

vector vector_add(vector lhs, vector rhs)
{
    return (vector){
        lhs.x + rhs.x,
        lhs.y + rhs.y,
        lhs.z + rhs.z,
    };
}

vector vector_sub(vector lhs, vector rhs)
{
    return (vector){
        lhs.x - rhs.x,
        lhs.y - rhs.y,
        lhs.z - rhs.z,
    };
}

vector vector_mul(scalar lhs, vector rhs)
{
    return (vector){
        lhs * rhs.x,
        lhs * rhs.y,
        lhs * rhs.z,
    };
}

vector vector_div(vector lhs, scalar rhs)
{
    assert(!isclose(rhs, 0));
    return vector_mul(1 / rhs, lhs);
}

void vector_inc(vector *lhs, vector rhs)
{
    *lhs = vector_add(*lhs, rhs);
}

void vector_dec(vector *lhs, vector rhs)
{
    *lhs = vector_sub(*lhs, rhs);
}

void vector_scale(vector *lhs, scalar rhs)
{
    *lhs = vector_mul(rhs, *lhs);
}

vector vector_abs(vector vec)
{
    return (vector){
        fabs(vec.x),
        fabs(vec.y),
        fabs(vec.z),
    };
}

vector vector_sum(const vector *vec, number num)
{
    assert(vec ? (num >= 0) : (num == 0));
    vector sum = {0};
    for (number i = 0; i < num; i++) {
        vector_inc(&sum, vec[i]);
    }
    return sum;
}

scalar vector_dot(vector lhs, vector rhs)
{
    return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z);
}

scalar vector_norm(vector vec)
{
    return sqrt(vector_dot(vec, vec));
}

vector vector_normalize(vector vec)
{
    return vector_div(vec, vector_norm(vec));
}

vector vector_cross(vector lhs, vector rhs)
{
    return (vector){
        (lhs.y * rhs.z) - (lhs.z * rhs.y),
        (lhs.z * rhs.x) - (lhs.x * rhs.z),
        (lhs.x * rhs.y) - (lhs.y * rhs.x),
    };
}
