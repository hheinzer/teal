#include "boundary.h"

#include <string.h>

#include "euler/boundary.h"
#include "navierstokes.h"

static ApplyBC adiabatic_wall;

ApplyBC *navierstokes_bc(const char *name)
{
    if (!strcmp(name, "adiabatic wall") || !strcmp(name, "wall")) return adiabatic_wall;
    return euler_bc(name);
}

void adiabatic_wall(const Equations *, const Matrix3d, const double *, const double *ui, double *ug)
{
    // Blazek 2015, sec. 8.2.2
    ug[D] = ui[D];
    ug[U] = -ui[U];
    ug[V] = -ui[V];
    ug[W] = -ui[W];
    ug[P] = ui[P];
}
