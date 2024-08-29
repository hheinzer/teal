#pragma once

#include "mesh.h"

void reorder_cells(Mesh *mesh, const long *old2new, const long *new2old);
