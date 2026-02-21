#include <assert.h>
#include <math.h>

#include "utils2.h"
#include "vector2.h"

_Static_assert(sizeof(Vector) == sizeof(double[3]), "Vector must be packed");

Vector vector2_add(Vector lhs, Vector rhs)
{
    return (Vector){
        lhs.x + rhs.x,
        lhs.y + rhs.y,
        lhs.z + rhs.z,
    };
}

Vector vector2_sub(Vector lhs, Vector rhs)
{
    return (Vector){
        lhs.x - rhs.x,
        lhs.y - rhs.y,
        lhs.z - rhs.z,
    };
}

Vector vector2_mul(double lhs, Vector rhs)
{
    return (Vector){
        lhs * rhs.x,
        lhs * rhs.y,
        lhs * rhs.z,
    };
}

Vector vector2_div(Vector lhs, double rhs)
{
    assert(!isclose(rhs, 0));
    return (Vector){
        lhs.x / rhs,
        lhs.y / rhs,
        lhs.z / rhs,
    };
}

void vector2_iadd(Vector *lhs, Vector rhs)
{
    assert(lhs);
    lhs->x += rhs.x;
    lhs->y += rhs.y;
    lhs->z += rhs.z;
}

void vector2_isub(Vector *lhs, Vector rhs)
{
    assert(lhs);
    lhs->x -= rhs.x;
    lhs->y -= rhs.y;
    lhs->z -= rhs.z;
}

void vector2_imul(Vector *lhs, double rhs)
{
    assert(lhs);
    lhs->x *= rhs;
    lhs->y *= rhs;
    lhs->z *= rhs;
}

void vector2_idiv(Vector *lhs, double rhs)
{
    assert(lhs && !isclose(rhs, 0));
    lhs->x /= rhs;
    lhs->y /= rhs;
    lhs->z /= rhs;
}

Vector vector2_abs(Vector vec)
{
    return (Vector){
        fabs(vec.x),
        fabs(vec.y),
        fabs(vec.z),
    };
}

double vector2_dot(Vector lhs, Vector rhs)
{
    return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z);
}

double vector2_norm(Vector vec)
{
    return sqrt(vector2_dot(vec, vec));
}

double vector2_norm2(Vector vec)
{
    return vector2_dot(vec, vec);
}

Vector vector2_cross(Vector lhs, Vector rhs)
{
    return (Vector){
        (lhs.y * rhs.z) - (lhs.z * rhs.y),
        (lhs.z * rhs.x) - (lhs.x * rhs.z),
        (lhs.x * rhs.y) - (lhs.y * rhs.x),
    };
}
