#pragma once

#include "teal.h"

vector vector_min(vector lhs, vector rhs);
vector vector_max(vector lhs, vector rhs);

vector vector_add(vector lhs, vector rhs);
vector vector_sub(vector lhs, vector rhs);
vector vector_mul(vector lhs, scalar rhs);
vector vector_div(vector lhs, scalar rhs);

vector vector_abs(vector vec);

scalar vector_dot(vector lhs, vector rhs);
scalar vector_len(vector vec);

vector vector_unit(vector vec);
vector vector_cross(vector lhs, vector rhs);
