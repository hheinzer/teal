#pragma once

#include <parmetis.h>

#include "private.h"

typedef struct {
    int num;
    idx_t *tag;
    Vector *coord;
} GridNodes;

typedef struct {
    int num;
    int num_inner;
    int off_boundary;
    int off_periodic;
    idx_t *node_off;
    idx_t *node_idx;
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
    idx_t *node_off;
    idx_t *node_idx;
} GridEntities;

typedef struct {
    GridNodes nodes;
    GridCells cells;
    GridEntities entities;
} Grid;

Grid *grid_init(const char *fname);

void grid_deinit(Grid *grid);

idx_t *grid_tag_to_idx(const Grid *grid, const idx_t *node_tag, int num_tags);

Grid *grid_gmsh(const char *fname);
