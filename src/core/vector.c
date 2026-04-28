#include "vector.h"

#include <assert.h>
#include <math.h>

#include "utils.h"

_Static_assert(sizeof(Vector) == sizeof(double[3]), "Vector must be packed");

Vector vector_add(Vector lhs, Vector rhs)
{
    return (Vector){
        lhs.x + rhs.x,
        lhs.y + rhs.y,
        lhs.z + rhs.z,
    };
}

Vector vector_sub(Vector lhs, Vector rhs)
{
    return (Vector){
        lhs.x - rhs.x,
        lhs.y - rhs.y,
        lhs.z - rhs.z,
    };
}

Vector vector_mul(double lhs, Vector rhs)
{
    return (Vector){
        lhs * rhs.x,
        lhs * rhs.y,
        lhs * rhs.z,
    };
}

Vector vector_div(Vector lhs, double rhs)
{
    assert(!isclose(rhs, 0));
    return (Vector){
        lhs.x / rhs,
        lhs.y / rhs,
        lhs.z / rhs,
    };
}

void vector_iadd(Vector *lhs, Vector rhs)
{
    assert(lhs);
    lhs->x += rhs.x;
    lhs->y += rhs.y;
    lhs->z += rhs.z;
}

void vector_isub(Vector *lhs, Vector rhs)
{
    assert(lhs);
    lhs->x -= rhs.x;
    lhs->y -= rhs.y;
    lhs->z -= rhs.z;
}

void vector_imul(Vector *lhs, double rhs)
{
    assert(lhs);
    lhs->x *= rhs;
    lhs->y *= rhs;
    lhs->z *= rhs;
}

void vector_idiv(Vector *lhs, double rhs)
{
    assert(lhs && !isclose(rhs, 0));
    lhs->x /= rhs;
    lhs->y /= rhs;
    lhs->z /= rhs;
}

Vector vector_abs(Vector vec)
{
    return (Vector){
        fabs(vec.x),
        fabs(vec.y),
        fabs(vec.z),
    };
}

double vector_dot(Vector lhs, Vector rhs)
{
    return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z);
}

double vector_norm(Vector vec)
{
    return sqrt(vector_dot(vec, vec));
}

double vector_norm2(Vector vec)
{
    return vector_dot(vec, vec);
}

Vector vector_cross(Vector lhs, Vector rhs)
{
    return (Vector){
        (lhs.y * rhs.z) - (lhs.z * rhs.y),
        (lhs.z * rhs.x) - (lhs.x * rhs.z),
        (lhs.x * rhs.y) - (lhs.y * rhs.x),
    };
}
