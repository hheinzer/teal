#include "transform.h"

void transform_to_local(vector *res, const Basis *mat, const vector *vec)
{
    res->x = (mat->n.x * vec->x) + (mat->n.y * vec->y) + (mat->n.z * vec->z);
    res->y = (mat->s.x * vec->x) + (mat->s.y * vec->y) + (mat->s.z * vec->z);
    res->z = (mat->t.x * vec->x) + (mat->t.y * vec->y) + (mat->t.z * vec->z);
}

void transform_to_global(vector *res, const Basis *mat, const vector *vec)
{
    res->x = (mat->n.x * vec->x) + (mat->s.x * vec->y) + (mat->t.x * vec->z);
    res->y = (mat->n.y * vec->x) + (mat->s.y * vec->y) + (mat->t.y * vec->z);
    res->z = (mat->n.z * vec->x) + (mat->s.z * vec->y) + (mat->t.z * vec->z);
}
