#pragma once

#include "mesh.h"

void sync_all(const Mesh *mesh, double *u, long ldu);

void sync_mesh(Mesh *mesh, int rank);
