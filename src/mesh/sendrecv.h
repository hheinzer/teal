#pragma once

#include "mesh.h"

void send_mesh(Mesh *mesh, long rank);

void recv_mesh(Mesh *mesh, long root);
