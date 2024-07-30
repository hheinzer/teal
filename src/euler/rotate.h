#pragma once

/* Rotate the global state 'g' into the local state 'l', using the normal vector 'n'. */
void global_to_local(const double *n, const double *g, double *l);

/* Rotate the local state 'l' into the global state 'g', using the normal vector 'n'. */
void local_to_global(const double *n, const double *l, double *g);
