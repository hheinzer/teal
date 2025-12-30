#include "vector.h"

vector vector_cross(vector lhs, vector rhs)
{
    return (vector){
        .x = (lhs.y * rhs.z) - (lhs.z * rhs.y),
        .y = (lhs.z * rhs.x) - (lhs.x * rhs.z),
        .z = (lhs.x * rhs.y) - (lhs.y * rhs.x),
    };
}
