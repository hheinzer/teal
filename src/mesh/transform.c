#include "transform.h"

void transform_to_local(vector *res, const Basis *mat, const vector *vec)
{
    res->x = (mat->normal.x * vec->x) + (mat->normal.y * vec->y) + (mat->normal.z * vec->z);
    res->y = (mat->tangent1.x * vec->x) + (mat->tangent1.y * vec->y) + (mat->tangent1.z * vec->z);
    res->z = (mat->tangent2.x * vec->x) + (mat->tangent2.y * vec->y) + (mat->tangent2.z * vec->z);
}

void transform_to_global(vector *res, const Basis *mat, const vector *vec)
{
    res->x = (mat->normal.x * vec->x) + (mat->tangent1.x * vec->y) + (mat->tangent2.x * vec->z);
    res->y = (mat->normal.y * vec->x) + (mat->tangent1.y * vec->y) + (mat->tangent2.y * vec->z);
    res->z = (mat->normal.z * vec->x) + (mat->tangent1.z * vec->y) + (mat->tangent2.z * vec->z);
}
