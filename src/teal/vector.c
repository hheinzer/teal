#include "vector.h"

#include <assert.h>
#include <math.h>

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
    assert(!is_close(rhs, 0));
    return vector_mul(1 / rhs, lhs);
}

vector vector_abs(vector vec)
{
    return (vector){
        fabs(vec.x),
        fabs(vec.y),
        fabs(vec.z),
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

void vector_scale(vector *lhs, scalar rhs)
{
    lhs->x *= rhs;
    lhs->y *= rhs;
    lhs->z *= rhs;
}

vector vector_min(const vector *vec, long num)
{
    assert(vec ? (num >= 0) : (num == 0));
    vector min = {INFINITY, INFINITY, INFINITY};
    for (long i = 0; i < num; i++) {
        min.x = fmin(min.x, vec[i].x);
        min.y = fmin(min.y, vec[i].y);
        min.z = fmin(min.z, vec[i].z);
    }
    return min;
}

vector vector_max(const vector *vec, long num)
{
    assert(vec ? (num >= 0) : (num == 0));
    vector max = {-INFINITY, -INFINITY, -INFINITY};
    for (long i = 0; i < num; i++) {
        max.x = fmax(max.x, vec[i].x);
        max.y = fmax(max.y, vec[i].y);
        max.z = fmax(max.z, vec[i].z);
    }
    return max;
}

vector vector_sum(const vector *vec, long num)
{
    assert(vec ? (num >= 0) : (num == 0));
    vector sum = {0};
    for (long i = 0; i < num; i++) {
        sum.x += vec[i].x;
        sum.y += vec[i].y;
        sum.z += vec[i].z;
    }
    return sum;
}

scalar vector_dot(vector lhs, vector rhs)
{
    return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z);
}

scalar vector_subdot(vector lhs, vector rhs, vector dot)
{
    return ((lhs.x - rhs.x) * dot.x) + ((lhs.y - rhs.y) * dot.y) + ((lhs.z - rhs.z) * dot.z);
}

scalar vector_norm(vector vec)
{
    return sqrt(vector_dot(vec, vec));
}

vector vector_cross(vector lhs, vector rhs)
{
    return (vector){
        (lhs.y * rhs.z) - (lhs.z * rhs.y),
        (lhs.z * rhs.x) - (lhs.x * rhs.z),
        (lhs.x * rhs.y) - (lhs.y * rhs.x),
    };
}
