#include <assert.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "mesh.h"
#include "parse.h"
#include "sync.h"
#include "teal.h"
#include "utils.h"

typedef struct {
    double version;
    int32_t file_type;
    int32_t data_size;
} Format;

typedef struct {
    int32_t dim;
    int32_t tag;
    char name[128];
} Physical;

typedef struct {
    int32_t num;
    Physical *physical;
} Physicals;

typedef struct {
    int32_t tag;
    double x, y, z;
    uint64_t num_physical_tags;
    int32_t *physical_tag;
} Point;

typedef struct {
    int32_t tag;
    double min_x, min_y, min_z, max_x, max_y, max_z;
    uint64_t num_physical_tags;
    int32_t *physical_tag;
    uint64_t num_point_tags;
    int32_t *point_tag;
} Curve;

typedef struct {
    int32_t tag;
    double min_x, min_y, min_z, max_x, max_y, max_z;
    uint64_t num_physical_tags;
    int32_t *physical_tag;
    uint64_t num_curve_tags;
    int32_t *curve_tag;
} Surface;

typedef struct {
    int32_t tag;
    double min_x, min_y, min_z, max_x, max_y, max_z;
    uint64_t num_physical_tags;
    int32_t *physical_tag;
    uint64_t num_surface_tags;
    int32_t *surface_tag;
} Volume;

typedef struct {
    uint64_t num_points;
    uint64_t num_curves;
    uint64_t num_surfaces;
    uint64_t num_volumes;
    Point *point;
    Curve *curve;
    Surface *surface;
    Volume *volume;
} Entities;

typedef struct {
    int32_t entity_dim;
    int32_t entity_tag;
    int32_t parametric;
    int len_coord;
    uint64_t num_nodes;
    uint64_t *tag;
    double *coord;
} NodeBlock;

typedef struct {
    uint64_t num_blocks;
    uint64_t tot_nodes;
    uint64_t min_tag;
    uint64_t max_tag;
    NodeBlock *block;
} Nodes;

typedef struct {
    int32_t entity_dim;
    int32_t entity_tag;
    int32_t element_type;
    int len_node_tag;
    uint64_t num_elements;
    uint64_t *tag;
    uint64_t *node_tag;
} ElementBlock;

typedef struct {
    uint64_t num_blocks;
    uint64_t tot_elements;
    uint64_t min_tag;
    uint64_t max_tag;
    ElementBlock *block;
} Elements;

typedef struct {
    int32_t entity_dim;
    int32_t entity_tag;
    int32_t entity_tag_master;
    uint64_t num_affine;
    double *affine;
    uint64_t num_nodes;
    uint64_t *node_tag;
    uint64_t *node_tag_master;
} Periodic;

typedef struct {
    uint64_t num;
    Periodic *periodic;
} Periodics;

typedef struct {
    Format format;
    Physicals physicals;
    Entities entities;
    Nodes nodes;
    Elements elements;
    Periodics periodics;
} Gmsh;

static int read_format(Gmsh *gmsh, Parse *file)
{
    parse_ascii(file, &gmsh->format.version, 1, MPI_DOUBLE);
    if (!isclose(gmsh->format.version, 4.1)) {
        teal_error("unsuppored version (%g)", gmsh->format.version);
    }

    parse_ascii(file, &gmsh->format.file_type, 1, MPI_INT32_T);
    int mode = (gmsh->format.file_type == 0) ? ASCII : BINARY;

    parse_ascii(file, &gmsh->format.data_size, 1, MPI_INT32_T);
    if (gmsh->format.data_size != 8) {
        teal_error("unsuppored data-size (%d)", gmsh->format.data_size);
    }

    if (gmsh->format.file_type == 1) {
        int one;
        parse_binary(file, &one, 1, MPI_INT32_T, 0);
        if (one != 1) {
            mode |= SWAP;
        }
    }

    return mode;
}

static void read_physicals(Gmsh *gmsh, Parse *file)
{
    parse_ascii(file, &gmsh->physicals.num, 1, MPI_INT32_T);

    Physical *physical = teal_alloc(gmsh->physicals.num, sizeof(*physical));
    for (int i = 0; i < gmsh->physicals.num; i++) {
        parse_ascii(file, &physical[i].dim, 1, MPI_INT32_T);
        parse_ascii(file, &physical[i].tag, 1, MPI_INT32_T);
        parse_string(file, physical[i].name, sizeof(physical->name));
    }
    gmsh->physicals.physical = physical;
}

static Point *read_points(int num, Parse *file, int mode)
{
    Point *point = teal_alloc(num, sizeof(*point));
    for (int i = 0; i < num; i++) {
        parse(file, &point[i].tag, 1, MPI_INT32_T, mode);
        parse(file, &point[i].x, 1, MPI_DOUBLE, mode);
        parse(file, &point[i].y, 1, MPI_DOUBLE, mode);
        parse(file, &point[i].z, 1, MPI_DOUBLE, mode);

        parse(file, &point[i].num_physical_tags, 1, MPI_UINT64_T, mode);
        assert(point[i].num_physical_tags <= INT_MAX);
        int num_physical_tags = (int)point[i].num_physical_tags;

        int32_t *physical_tag = teal_alloc(num_physical_tags, sizeof(*physical_tag));
        parse(file, physical_tag, num_physical_tags, MPI_INT32_T, mode);
        point[i].physical_tag = physical_tag;
    }
    return point;
}

static Curve *read_curves(int num, Parse *file, int mode)
{
    Curve *curve = teal_alloc(num, sizeof(*curve));
    for (int i = 0; i < num; i++) {
        parse(file, &curve[i].tag, 1, MPI_INT32_T, mode);
        parse(file, &curve[i].min_x, 1, MPI_DOUBLE, mode);
        parse(file, &curve[i].min_y, 1, MPI_DOUBLE, mode);
        parse(file, &curve[i].min_z, 1, MPI_DOUBLE, mode);
        parse(file, &curve[i].max_x, 1, MPI_DOUBLE, mode);
        parse(file, &curve[i].max_y, 1, MPI_DOUBLE, mode);
        parse(file, &curve[i].max_z, 1, MPI_DOUBLE, mode);

        parse(file, &curve[i].num_physical_tags, 1, MPI_UINT64_T, mode);
        assert(curve[i].num_physical_tags <= INT_MAX);
        int num_physical_tags = (int)curve[i].num_physical_tags;

        int32_t *physical_tag = teal_alloc(num_physical_tags, sizeof(*physical_tag));
        parse(file, physical_tag, num_physical_tags, MPI_INT32_T, mode);
        curve[i].physical_tag = physical_tag;

        parse(file, &curve[i].num_point_tags, 1, MPI_UINT64_T, mode);
        assert(curve[i].num_point_tags <= INT_MAX);
        int num_point_tags = (int)curve[i].num_point_tags;

        int32_t *point_tag = teal_alloc(num_point_tags, sizeof(*point_tag));
        parse(file, point_tag, num_point_tags, MPI_INT32_T, mode);
        curve[i].point_tag = point_tag;
    }
    return curve;
}

static Surface *read_surfaces(int num, Parse *file, int mode)
{
    Surface *surface = teal_alloc(num, sizeof(*surface));
    for (int i = 0; i < num; i++) {
        parse(file, &surface[i].tag, 1, MPI_INT32_T, mode);
        parse(file, &surface[i].min_x, 1, MPI_DOUBLE, mode);
        parse(file, &surface[i].min_y, 1, MPI_DOUBLE, mode);
        parse(file, &surface[i].min_z, 1, MPI_DOUBLE, mode);
        parse(file, &surface[i].max_x, 1, MPI_DOUBLE, mode);
        parse(file, &surface[i].max_y, 1, MPI_DOUBLE, mode);
        parse(file, &surface[i].max_z, 1, MPI_DOUBLE, mode);

        parse(file, &surface[i].num_physical_tags, 1, MPI_UINT64_T, mode);
        assert(surface[i].num_physical_tags <= INT_MAX);
        int num_physical_tags = (int)surface[i].num_physical_tags;

        int32_t *physical_tag = teal_alloc(num_physical_tags, sizeof(*physical_tag));
        parse(file, physical_tag, num_physical_tags, MPI_INT32_T, mode);
        surface[i].physical_tag = physical_tag;

        parse(file, &surface[i].num_curve_tags, 1, MPI_UINT64_T, mode);
        assert(surface[i].num_curve_tags <= INT_MAX);
        int num_curve_tags = (int)surface[i].num_curve_tags;

        int32_t *curve_tag = teal_alloc(num_curve_tags, sizeof(*curve_tag));
        parse(file, curve_tag, num_curve_tags, MPI_INT32_T, mode);
        surface[i].curve_tag = curve_tag;
    }
    return surface;
}

static Volume *read_volumes(int num, Parse *file, int mode)
{
    Volume *volume = teal_alloc(num, sizeof(*volume));
    for (int i = 0; i < num; i++) {
        parse(file, &volume[i].tag, 1, MPI_INT32_T, mode);
        parse(file, &volume[i].min_x, 1, MPI_DOUBLE, mode);
        parse(file, &volume[i].min_y, 1, MPI_DOUBLE, mode);
        parse(file, &volume[i].min_z, 1, MPI_DOUBLE, mode);
        parse(file, &volume[i].max_x, 1, MPI_DOUBLE, mode);
        parse(file, &volume[i].max_y, 1, MPI_DOUBLE, mode);
        parse(file, &volume[i].max_z, 1, MPI_DOUBLE, mode);

        parse(file, &volume[i].num_physical_tags, 1, MPI_UINT64_T, mode);
        assert(volume[i].num_physical_tags <= INT_MAX);
        int num_physical_tags = (int)volume[i].num_physical_tags;

        int32_t *physical_tag = teal_alloc(num_physical_tags, sizeof(*physical_tag));
        parse(file, physical_tag, num_physical_tags, MPI_INT32_T, mode);
        volume[i].physical_tag = physical_tag;

        parse(file, &volume[i].num_surface_tags, 1, MPI_UINT64_T, mode);
        assert(volume[i].num_surface_tags <= INT_MAX);
        int num_surface_tags = (int)volume[i].num_surface_tags;

        int32_t *surface_tag = teal_alloc(num_surface_tags, sizeof(*surface_tag));
        parse(file, surface_tag, num_surface_tags, MPI_INT32_T, mode);
        volume[i].surface_tag = surface_tag;
    }
    return volume;
}

static void read_entities(Gmsh *gmsh, Parse *file, int mode)
{
    parse(file, &gmsh->entities.num_points, 1, MPI_UINT64_T, mode);
    assert(gmsh->entities.num_points <= INT_MAX);
    int num_points = (int)gmsh->entities.num_points;

    parse(file, &gmsh->entities.num_curves, 1, MPI_UINT64_T, mode);
    assert(gmsh->entities.num_curves <= INT_MAX);
    int num_curves = (int)gmsh->entities.num_curves;

    parse(file, &gmsh->entities.num_surfaces, 1, MPI_UINT64_T, mode);
    assert(gmsh->entities.num_surfaces <= INT_MAX);
    int num_surfaces = (int)gmsh->entities.num_surfaces;

    parse(file, &gmsh->entities.num_volumes, 1, MPI_UINT64_T, mode);
    assert(gmsh->entities.num_volumes <= INT_MAX);
    int num_volumes = (int)gmsh->entities.num_volumes;

    gmsh->entities.point = read_points(num_points, file, mode);
    gmsh->entities.curve = read_curves(num_curves, file, mode);
    gmsh->entities.surface = read_surfaces(num_surfaces, file, mode);
    gmsh->entities.volume = read_volumes(num_volumes, file, mode);
}

static long read_node_block(NodeBlock *block, long beg, long end, long off, Parse *file, int mode)
{
    parse(file, &block->entity_dim, 1, MPI_INT32_T, mode);
    assert(inrange(0, block->entity_dim, 3));

    parse(file, &block->entity_tag, 1, MPI_INT32_T, mode);
    parse(file, &block->parametric, 1, MPI_INT32_T, mode);
    block->len_coord = 3 + (block->parametric ? block->entity_dim : 0);

    parse(file, &block->num_nodes, 1, MPI_UINT64_T, mode);
    assert(block->num_nodes <= LONG_MAX);
    long tot_nodes = (long)block->num_nodes;

    long beg_nodes = max(beg, off);
    long end_nodes = min(end, off + tot_nodes);
    long num_nodes = max(0, end_nodes - beg_nodes);
    assert(num_nodes <= INT_MAX && tot_nodes == sync_lsum((int)num_nodes));
    block->num_nodes = (uint64_t)num_nodes;

    block->tag = teal_alloc((int)num_nodes, sizeof(*block->tag));
    parse(file, block->tag, (int)num_nodes, MPI_UINT64_T, mode | SPLIT);

    assert(num_nodes < INT_MAX / block->len_coord);
    int tot_coords = (int)num_nodes * block->len_coord;
    block->coord = teal_alloc(tot_coords, sizeof(*block->coord));
    parse(file, block->coord, tot_coords, MPI_DOUBLE, mode | SPLIT);

    return tot_nodes;
}

static void read_nodes(Gmsh *gmsh, Parse *file, int mode)
{
    parse(file, &gmsh->nodes.num_blocks, 1, MPI_UINT64_T, mode);
    assert(gmsh->nodes.num_blocks <= INT_MAX);
    int num_blocks = (int)gmsh->nodes.num_blocks;

    parse(file, &gmsh->nodes.tot_nodes, 1, MPI_UINT64_T, mode);
    assert(gmsh->nodes.tot_nodes <= LONG_MAX);
    long tot_nodes = (long)gmsh->nodes.tot_nodes;

    parse(file, &gmsh->nodes.min_tag, 1, MPI_UINT64_T, mode);
    parse(file, &gmsh->nodes.max_tag, 1, MPI_UINT64_T, mode);

    long base = tot_nodes / sync.size;
    long extra = tot_nodes % sync.size;
    long beg = (sync.rank * base) + ((sync.rank < extra) ? sync.rank : extra);
    long end = beg + base + (sync.rank < extra);

    NodeBlock *block = teal_alloc(num_blocks, sizeof(*block));
    long off = 0;
    for (int i = 0; i < num_blocks; i++) {
        off += read_node_block(&block[i], beg, end, off, file, mode);
    }
    assert(off == tot_nodes);
    gmsh->nodes.block = block;
}

static int32_t num_node_tags(int32_t element_type)
{
    switch (element_type) {
        case 2: return 3;  // triangle
        case 3:            // quadrangle
        case 4: return 4;  // tetrahedron
        case 5: return 8;  // hexahedron
        case 6: return 6;  // prism
        case 7: return 5;  // pyramid
        default: teal_error("unsupported element type (%d)", element_type);
    }
}

static long read_element_block(ElementBlock *block, long beg, long end, long off, Parse *file,
                               int mode)
{
    parse(file, &block->entity_dim, 1, MPI_INT32_T, mode);
    assert(inrange(0, block->entity_dim, 3));

    parse(file, &block->entity_tag, 1, MPI_INT32_T, mode);
    parse(file, &block->element_type, 1, MPI_INT32_T, mode);
    block->len_node_tag = num_node_tags(block->element_type);

    parse(file, &block->num_elements, 1, MPI_UINT64_T, mode);
    assert(block->num_elements <= LONG_MAX);
    long tot_elements = (long)block->num_elements;

    long beg_elements = max(beg, off);
    long end_elements = min(end, off + tot_elements);
    long num_elements = max(0, end_elements - beg_elements);
    assert(num_elements <= INT_MAX && tot_elements == sync_lsum(num_elements));
    block->num_elements = (uint64_t)num_elements;

    block->tag = teal_alloc((int)num_elements, sizeof(*block->tag));

    assert(num_elements < INT_MAX / block->len_node_tag);
    int tot_node_tags = (int)num_elements * block->len_node_tag;
    block->node_tag = teal_alloc(tot_node_tags, sizeof(*block->node_tag));

    int stride = 1 + block->len_node_tag;
    assert(num_elements <= INT_MAX / stride);
    uint64_t (*tmp)[stride] = teal_alloc((int)num_elements, sizeof(*tmp));
    parse(file, tmp, (int)num_elements * stride, MPI_UINT64_T, mode | SPLIT);
    for (int i = 0; i < num_elements; i++) {
        block->tag[i] = tmp[i][0];
        for (int j = 0; j < block->len_node_tag; j++) {
            block->node_tag[(i * block->len_node_tag) + j] = tmp[i][1 + j];
        }
    }
    teal_free(tmp);

    return tot_elements;
}

static void read_elements(Gmsh *gmsh, Parse *file, int mode)
{
    parse(file, &gmsh->elements.num_blocks, 1, MPI_UINT64_T, mode);
    assert(gmsh->elements.num_blocks <= INT_MAX);
    int num_blocks = (int)gmsh->elements.num_blocks;

    parse(file, &gmsh->elements.tot_elements, 1, MPI_UINT64_T, mode);
    assert(gmsh->elements.tot_elements <= LONG_MAX);
    long tot_elements = (long)gmsh->elements.tot_elements;

    parse(file, &gmsh->elements.min_tag, 1, MPI_UINT64_T, mode);
    parse(file, &gmsh->elements.max_tag, 1, MPI_UINT64_T, mode);

    long base = tot_elements / sync.size;
    long extra = tot_elements % sync.size;
    long beg = (sync.rank * base) + ((sync.rank < extra) ? sync.rank : extra);
    long end = beg + base + (sync.rank < extra);

    ElementBlock *block = teal_alloc(num_blocks, sizeof(*block));
    long off = 0;
    for (int i = 0; i < num_blocks; i++) {
        off += read_element_block(&block[i], beg, end, off, file, mode);
    }
    assert(off == tot_elements);
    gmsh->elements.block = block;
}

static void read_link(Periodic *link, Parse *file, int mode)
{
    parse(file, &link->entity_dim, 1, MPI_INT32_T, mode);
    parse(file, &link->entity_tag, 1, MPI_INT32_T, mode);
    parse(file, &link->entity_tag_master, 1, MPI_INT32_T, mode);

    parse(file, &link->num_affine, 1, MPI_UINT64_T, mode);
    assert(link->num_affine <= INT_MAX);
    int num_affine = (int)link->num_affine;

    link->affine = teal_alloc(num_affine, sizeof(*link->affine));
    parse(file, link->affine, num_affine, MPI_DOUBLE, mode);

    parse(file, &link->num_nodes, 1, MPI_UINT64_T, mode);
    assert(link->num_nodes <= INT_MAX);
    int num_nodes = (int)link->num_nodes;

    link->node_tag = teal_alloc(num_nodes, sizeof(*link->node_tag));
    link->node_tag_master = teal_alloc(num_nodes, sizeof(*link->node_tag_master));

    assert(num_nodes <= INT_MAX / 2);
    uint64_t (*tmp)[2] = teal_alloc(num_nodes, sizeof(*tmp));
    parse(file, tmp, num_nodes * 2, MPI_UINT64_T, mode);
    for (int i = 0; i < num_nodes; i++) {
        link->node_tag[i] = tmp[i][0];
        link->node_tag_master[i] = tmp[i][1];
    }
    teal_free(tmp);
}

static void read_periodics(Gmsh *gmsh, Parse *file, int mode)
{
    parse(file, &gmsh->periodics.num, 1, MPI_UINT64_T, mode);
    assert(gmsh->periodics.num <= INT_MAX);
    int num = (int)gmsh->periodics.num;

    Periodic *link = teal_alloc(num, sizeof(*link));
    for (int i = 0; i < num; i++) {
        read_link(&link[i], file, mode);
    }
    gmsh->periodics.periodic = link;
}

static Gmsh *gmsh_init(const char *fname)
{
    Gmsh *gmsh = teal_alloc(1, sizeof(*gmsh));
    Parse *file = parse_init(fname);
    char section[128];
    int mode = 0;
    while (parse_string(file, section, sizeof(section))) {
        if (!strcmp(section, "$MeshFormat")) {
            mode = read_format(gmsh, file);
        }
        else if (!strcmp(section, "$PhysicalNames")) {
            read_physicals(gmsh, file);
        }
        else if (!strcmp(section, "$Entities")) {
            read_entities(gmsh, file, mode);
        }
        else if (!strcmp(section, "$Nodes")) {
            read_nodes(gmsh, file, mode);
        }
        else if (!strcmp(section, "$Elements")) {
            read_elements(gmsh, file, mode);
        }
        else if (!strcmp(section, "$Periodic")) {
            read_periodics(gmsh, file, mode);
        }
        else if (!strncmp(section, "$End", 4)) {
            continue;
        }
        else {
            teal_error("unexpected mesh section (%s)", section);
        }
    }
    parse_deinit(file);
    return gmsh;
}

static void create_nodes(Mesh *mesh, const Gmsh *gmsh)
{
    int num_nodes = 0;
    for (uint64_t i = 0; i < gmsh->nodes.num_blocks; i++) {
        const NodeBlock *block = &gmsh->nodes.block[i];
        assert(block->num_nodes <= INT_MAX);
        assert(num_nodes <= INT_MAX - (int)block->num_nodes);
        num_nodes += (int)block->num_nodes;
    }
    assert(num_nodes > 0);
    mesh->nodes.num = num_nodes;

    vector *coord = teal_alloc(num_nodes, sizeof(*coord));
    int num = 0;
    for (uint64_t i = 0; i < gmsh->nodes.num_blocks; i++) {
        NodeBlock *block = &gmsh->nodes.block[i];
        double (*block_coord)[block->len_coord] = (void *)block->coord;
        for (uint64_t j = 0; j < block->num_nodes; j++) {
            coord[num].x = block_coord[j][0];
            coord[num].y = block_coord[j][1];
            coord[num].z = block_coord[j][2];
            num += 1;
        }
    }
    assert(num == num_nodes);
    mesh->nodes.coord = coord;
}

static long *tag_to_idx(long *tag, const Mesh *mesh, const Gmsh *gmsh)
{
    typedef struct {
        long tag, idx;
    } Map;

    int cap = sync_max(mesh->nodes.num);
    Map *map = teal_alloc(cap, sizeof(*map));

    long prefix = sync_lexsum(mesh->nodes.num);
    int num = 0;
    for (uint64_t i = 0; i < gmsh->nodes.num_blocks; i++) {
        NodeBlock *block = &gmsh->nodes.block[i];
        for (uint64_t j = 0; j < block->num_nodes; j++) {
            assert(inrange(1, block->tag[j], LONG_MAX));
            map[num].tag = (long)block->tag[j];
            map[num].idx = prefix + num;
            num += 1;
        }
    }
    assert(num == mesh->nodes.num);
    qsort(map, (size_t)num, sizeof(*map), cmp_long);

    int num_tags = mesh->cells.node.off[mesh->cells.num];
    long *idx = teal_alloc(num_tags, sizeof(*idx));

    MPI_Datatype datatype;
    MPI_Type_contiguous(2, MPI_LONG, &datatype);
    MPI_Type_commit(&datatype);

    int next = (sync.rank + 1) % sync.size;
    int prev = (sync.rank - 1 + sync.size) % sync.size;
    for (int rank = 0; rank < sync.size; rank++) {
        for (int i = 0; i < num_tags; i++) {
            if (tag[i] > 0) {
                Map key = {.tag = tag[i]};
                const Map *val = bsearch(&key, map, (size_t)num, sizeof(*map), cmp_long);
                if (val) {
                    idx[i] = val->idx;
                    tag[i] = -tag[i];
                }
            }
        }
        MPI_Sendrecv_replace(&num, 1, MPI_INT, next, 0, prev, 0, sync.comm, MPI_STATUS_IGNORE);
        MPI_Sendrecv_replace(map, cap, datatype, next, 1, prev, 1, sync.comm, MPI_STATUS_IGNORE);
    }
    MPI_Type_free(&datatype);

    teal_free(map);

    for (int i = 0; i < num_tags; i++) {
        assert(tag[i] < 0);
    }
    teal_free(tag);

    return idx;
}

static void create_inner_cells(Mesh *mesh, const Gmsh *gmsh)
{
    int num_cells = 0;
    for (uint64_t i = 0; i < gmsh->elements.num_blocks; i++) {
        const ElementBlock *block = &gmsh->elements.block[i];
        if (block->entity_dim != 3) {
            continue;
        }
        assert(block->num_elements <= INT_MAX);
        assert(num_cells <= INT_MAX - (int)block->num_elements);
        num_cells += (int)block->num_elements;
    }
    assert(num_cells > 0);
    mesh->cells.num = num_cells;

    int *node_off = teal_alloc(num_cells + 1, sizeof(*node_off));
    long *node_tag = teal_alloc(num_cells * MAX_CELL_NODES, sizeof(*node_tag));
    long num = 0;
    for (uint64_t i = 0; i < gmsh->elements.num_blocks; i++) {
        ElementBlock *block = &gmsh->elements.block[i];
        if (block->entity_dim != 3) {
            continue;
        }
        int len = block->len_node_tag;
        uint64_t (*block_node_tag)[len] = (void *)block->node_tag;
        for (uint64_t j = 0; j < block->num_elements; j++) {
            node_off[num + 1] = node_off[num] + len;
            for (int k = 0; k < len; k++) {
                assert(inrange(1, block_node_tag[j][k], LONG_MAX));
                node_tag[node_off[num] + k] = (long)block_node_tag[j][k];
            }
            num += 1;
        }
    }
    assert(num == num_cells);
    mesh->cells.node.off = node_off;
    mesh->cells.node.idx = tag_to_idx(node_tag, mesh, gmsh);
}

static void create_boundary_faces(Mesh *mesh, const Gmsh *gmsh)
{
}

static int cmp_name(const void *lhs_, const void *rhs_)
{
    const Physical *lhs = lhs_;
    const Physical *rhs = rhs_;
    if (lhs->dim == rhs->dim) {
        return cmp_int(&lhs->tag, &rhs->tag);  // ascending
    }
    return -cmp_int(&lhs->dim, &rhs->dim);  // descending
}

static void create_entities(Mesh *mesh, const Gmsh *gmsh)
{
    int num_entities = gmsh->physicals.num;
    assert(num_entities > 0);
    mesh->entities.num = num_entities;

    Physical *physical = memdup(gmsh->physicals.physical, num_entities, sizeof(*physical));
    qsort(physical, (size_t)num_entities, sizeof(*physical), cmp_name);

    Name *name = teal_alloc(num_entities, sizeof(*name));
    int num_inner = 0;
    int num_boundary = 0;
    for (int i = 0; i < num_entities; i++) {
        strcpy(name[i], physical[i].name);
        if (physical[i].dim == 3) {
            num_inner += 1;
        }
        else if (physical[i].dim == 2) {
            num_boundary += 1;
        }
        else {
            teal_error("invalid physical name (%d, %d, %s)", physical[i].dim, physical[i].tag,
                       physical[i].name);
        }
    }
    assert(num_boundary == num_entities - num_inner);
    mesh->entities.num_inner = num_inner;
    mesh->entities.name = name;

    teal_free(physical);
}

static void create_periodics(Mesh *mesh, const Gmsh *gmsh)
{
}

static void reorder_cells(Mesh *mesh, const Gmsh *gmsh)
{
}

static void gmsh_deinit(Gmsh *gmsh)
{
    if (!gmsh) {
        return;
    }

    teal_free(gmsh->physicals.physical);

    if (gmsh->entities.point) {
        for (uint64_t i = 0; i < gmsh->entities.num_points; i++) {
            teal_free(gmsh->entities.point[i].physical_tag);
        }
        teal_free(gmsh->entities.point);
    }
    if (gmsh->entities.curve) {
        for (uint64_t i = 0; i < gmsh->entities.num_curves; i++) {
            teal_free(gmsh->entities.curve[i].physical_tag);
            teal_free(gmsh->entities.curve[i].point_tag);
        }
        teal_free(gmsh->entities.curve);
    }
    if (gmsh->entities.surface) {
        for (uint64_t i = 0; i < gmsh->entities.num_surfaces; i++) {
            teal_free(gmsh->entities.surface[i].physical_tag);
            teal_free(gmsh->entities.surface[i].curve_tag);
        }
        teal_free(gmsh->entities.surface);
    }
    if (gmsh->entities.volume) {
        for (uint64_t i = 0; i < gmsh->entities.num_volumes; i++) {
            teal_free(gmsh->entities.volume[i].physical_tag);
            teal_free(gmsh->entities.volume[i].surface_tag);
        }
        teal_free(gmsh->entities.volume);
    }

    if (gmsh->nodes.block) {
        for (uint64_t i = 0; i < gmsh->nodes.num_blocks; i++) {
            teal_free(gmsh->nodes.block[i].tag);
            teal_free(gmsh->nodes.block[i].coord);
        }
        teal_free(gmsh->nodes.block);
    }

    if (gmsh->elements.block) {
        for (uint64_t i = 0; i < gmsh->elements.num_blocks; i++) {
            teal_free(gmsh->elements.block[i].tag);
            teal_free(gmsh->elements.block[i].node_tag);
        }
        teal_free(gmsh->elements.block);
    }

    if (gmsh->periodics.periodic) {
        for (uint64_t i = 0; i < gmsh->periodics.num; i++) {
            teal_free(gmsh->periodics.periodic[i].affine);
            teal_free(gmsh->periodics.periodic[i].node_tag);
            teal_free(gmsh->periodics.periodic[i].node_tag_master);
        }
        teal_free(gmsh->periodics.periodic);
    }

    teal_free(gmsh);
}

static void read_gmsh(Mesh *mesh, const char *fname)
{
    Gmsh *gmsh = gmsh_init(fname);

    create_nodes(mesh, gmsh);
    create_inner_cells(mesh, gmsh);
    create_boundary_faces(mesh, gmsh);
    create_entities(mesh, gmsh);
    create_periodics(mesh, gmsh);

    reorder_cells(mesh, gmsh);

    gmsh_deinit(gmsh);
}
