#pragma once

#include "mesh.h"

/* Send a 'mesh' to the specified 'rank'. */
void send_mesh(Mesh *mesh, long rank);

/* Receive a 'mesh' from the specified 'root'. */
void recv_mesh(Mesh *mesh, long root);
