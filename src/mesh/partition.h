#pragma once

#include "mesh.h"

/* Partition the 'mesh' and distribute the partitions among all MPI ranks. */
void partition(Mesh *mesh);
