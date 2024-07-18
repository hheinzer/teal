#ifndef ROTATE_H
#define ROTATE_H

void global_to_local(const double *n, const double *g, double *l);

void local_to_global(const double *n, const double *l, double *g);

#endif
