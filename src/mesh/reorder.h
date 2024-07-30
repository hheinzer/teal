#pragma once

#include "mesh.h"

/* Reorder the cells of the 'mesh' according to 'old2new' and 'new2old' index maps. */
void reorder_cells(Mesh *mesh, const long *old2new, const long *new2old);
