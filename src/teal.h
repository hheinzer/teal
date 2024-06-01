#ifndef TEAL_H
#define TEAL_H

#include "airfoil.h"
#include "euler.h"
#include "mesh.h"
#include "simulation.h"

#define PI 3.14159265358979323846

void teal_init(int argc, char **argv);

void teal_finalize(void);

#endif
