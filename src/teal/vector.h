#pragma once

#include "teal.h"

vector vector_add(vector lhs, vector rhs);
vector vector_sub(vector lhs, vector rhs);
vector vector_mul(scalar lhs, vector rhs);
vector vector_div(vector lhs, scalar rhs);
vector vector_abs(vector vec);

void vector_inc(vector *lhs, vector rhs);
void vector_dec(vector *lhs, vector rhs);
void vector_scale(vector *lhs, scalar rhs);

vector vector_min(const vector *vec, long num);
vector vector_max(const vector *vec, long num);
vector vector_sum(const vector *vec, long num);

scalar vector_dot(vector lhs, vector rhs);
scalar vector_norm(vector vec);
vector vector_cross(vector lhs, vector rhs);
