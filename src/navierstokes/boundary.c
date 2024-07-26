#include "boundary.h"

#include "navierstokes.h"

void adiabatic_wall(const Equations *, const double *, const double *, const double *ui, double *ug)
{
    // Blazek 2015, sec. 8.2.2
    ug[D] = ui[D];
    ug[U] = -ui[U];
    ug[V] = -ui[V];
    ug[W] = -ui[W];
    ug[P] = ui[P];
}
