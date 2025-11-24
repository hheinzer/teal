#pragma once

#include "mesh.h"

void transform_to_local(vector *res, const Basis *mat, const vector *vec);

void transform_to_global(vector *res, const Basis *mat, const vector *vec);
