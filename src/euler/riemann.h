#pragma once

/* Solve the exact Riemann problem for the left and right density, velocity, and pressure 'dl',
 * 'ul', 'pl', 'dr', 'ur', 'pr', the interface coordinate 's', and the heat capacity ratio 'gamma'.
 * The resulting density, velocity, and pressure will be written to 'd', 'u', and 'p'. */
void riemann(double *d, double *u, double *p, double dl, double ul, double pl, double dr, double ur,
             double pr, double s, double gamma);
