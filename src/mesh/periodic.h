#pragma once

#include "mesh.h"
#include "teal/dict.h"

void periodic_cell_to_node(const Mesh *mesh, Dict *periodic);

void periodic_cell_to_cell(const Mesh *mesh, Dict *periodic);
