#ifndef SENDRECV_H
#define SENDRECV_H

#include "mesh.h"

void send_mesh(Mesh *mesh, long rank);

void recv_mesh(Mesh *mesh, long root);

#endif
