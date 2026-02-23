#pragma once

#include "private.h"

typedef struct {
    int num;
    long *tag;
    Vector *coord;
} GridNodes;

typedef struct {
    int num;
    int num_inner;
    int off_boundary;
    int off_periodic;
    long *node_off;
    long *node_idx;
} GridCells;

typedef struct {
    int num;
    int num_inner;
    int off_boundary;
    Name *name;
    int *cell_off;
    int *periodic;
    Matrix *rotation;
    Vector *translation;
    long *node_off;
    long *node_idx;
} GridEntities;

typedef struct {
    GridNodes nodes;
    GridCells cells;
    GridEntities entities;
} Grid;

Grid *grid_init(const char *fname);

void grid_deinit(Grid *grid);

long *grid_tag_to_idx(const Grid *grid, const long *node_tag, int num_tags);

Grid *grid_gmsh(const char *fname);
