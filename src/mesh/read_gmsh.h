#include <assert.h>
#include <inttypes.h>
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
    uint64_t num_nodes;
    uint64_t *tag;
    double *coord;
} NodeBlock;

typedef struct {
    uint64_t num_blocks;
    uint64_t num_nodes;
    uint64_t min_tag;
    uint64_t max_tag;
    NodeBlock *block;
} Nodes;

typedef struct {
    int32_t entity_dim;
    int32_t entity_tag;
    int32_t element_type;
    uint64_t num_elements;
    uint64_t *tag;
    uint64_t *node_tag;
} ElementBlock;

typedef struct {
    uint64_t num_blocks;
    uint64_t num_elements;
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

// Read the mesh format and return the parse mode.
static int read_format(Gmsh *gmsh, Parse *file)
{
    parse_ascii(file, &gmsh->format.version, 1, MPI_DOUBLE);
    if (!isclose(gmsh->format.version, 4.1)) {
        teal_error("unsupported version (%g)", gmsh->format.version);
    }

    parse_ascii(file, &gmsh->format.file_type, 1, MPI_INT32_T);
    int mode = (gmsh->format.file_type == 0) ? ASCII : BINARY;

    parse_ascii(file, &gmsh->format.data_size, 1, MPI_INT32_T);
    if (gmsh->format.data_size != 8) {
        teal_error("unsupported data-size (%d)", gmsh->format.data_size);
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

// Read physical names into the gmsh struct.
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

// Read point entities into a new array.
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

// Read curve entities into a new array.
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

// Read surface entities into a new array.
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

// Read volume entities into a new array.
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

// Read the entities section into gmsh.
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

// Return the coordinate stride for a node block.
static int len_coord(const NodeBlock *block)
{
    return 3 + (block->parametric ? block->entity_dim : 0);
}

// Read one node block and return its original size.
static long read_node_block(NodeBlock *block, long beg_nodes, long end_nodes, long off, Parse *file,
                            int mode)
{
    parse(file, &block->entity_dim, 1, MPI_INT32_T, mode);
    assert(0 <= block->entity_dim && block->entity_dim <= 3);

    parse(file, &block->entity_tag, 1, MPI_INT32_T, mode);
    parse(file, &block->parametric, 1, MPI_INT32_T, mode);

    parse(file, &block->num_nodes, 1, MPI_UINT64_T, mode);
    assert(block->num_nodes <= LONG_MAX);
    long tot_nodes = (long)block->num_nodes;

    long beg = lmax(beg_nodes, off);
    long end = lmin(end_nodes, off + tot_nodes);
    long num_nodes = lmax(0, end - beg);
    assert(tot_nodes == sync_lsum(num_nodes));
    block->num_nodes = (uint64_t)num_nodes;

    assert(num_nodes <= INT_MAX);
    block->tag = teal_alloc((int)num_nodes, sizeof(*block->tag));
    parse(file, block->tag, (int)num_nodes, MPI_UINT64_T, mode | SPLIT);

    int len = len_coord(block);
    assert(num_nodes < INT_MAX / len);
    int tot_coords = (int)num_nodes * len;
    block->coord = teal_alloc(tot_coords, sizeof(*block->coord));
    parse(file, block->coord, tot_coords, MPI_DOUBLE, mode | SPLIT);

    return tot_nodes;
}

// Read the nodes section into gmsh.
static void read_nodes(Gmsh *gmsh, Parse *file, int mode)
{
    parse(file, &gmsh->nodes.num_blocks, 1, MPI_UINT64_T, mode);
    assert(gmsh->nodes.num_blocks <= INT_MAX);
    int num_blocks = (int)gmsh->nodes.num_blocks;

    parse(file, &gmsh->nodes.num_nodes, 1, MPI_UINT64_T, mode);
    assert(gmsh->nodes.num_nodes <= LONG_MAX);
    long tot_nodes = (long)gmsh->nodes.num_nodes;

    parse(file, &gmsh->nodes.min_tag, 1, MPI_UINT64_T, mode);
    parse(file, &gmsh->nodes.max_tag, 1, MPI_UINT64_T, mode);

    long base = tot_nodes / sync.size;
    long extra = tot_nodes % sync.size;
    long beg_nodes = (sync.rank * base) + ((sync.rank < extra) ? sync.rank : extra);
    long end_nodes = beg_nodes + base + (sync.rank < extra);
    long num_nodes = end_nodes - beg_nodes;
    assert(num_nodes >= 0);
    gmsh->nodes.num_nodes = (uint64_t)num_nodes;

    NodeBlock *block = teal_alloc(num_blocks, sizeof(*block));
    long off = 0;
    for (int i = 0; i < num_blocks; i++) {
        off += read_node_block(&block[i], beg_nodes, end_nodes, off, file, mode);
    }
    assert(off == tot_nodes);
    gmsh->nodes.block = block;
}

// Return the node tag count for an element block.
static int len_node_tags(const ElementBlock *block)
{
    switch (block->element_type) {
        case 2: return 3;  // triangle
        case 3:            // quadrangle
        case 4: return 4;  // tetrahedron
        case 5: return 8;  // hexahedron
        case 6: return 6;  // prism
        case 7: return 5;  // pyramid
        default: teal_error("unsupported element type (%d)", block->element_type);
    }
}

// Read one element block and return its original size.
static long read_element_block(ElementBlock *block, long beg_elements, long end_elements, long off,
                               Parse *file, int mode)
{
    parse(file, &block->entity_dim, 1, MPI_INT32_T, mode);
    assert(0 <= block->entity_dim && block->entity_dim <= 3);

    parse(file, &block->entity_tag, 1, MPI_INT32_T, mode);
    parse(file, &block->element_type, 1, MPI_INT32_T, mode);

    parse(file, &block->num_elements, 1, MPI_UINT64_T, mode);
    assert(block->num_elements <= LONG_MAX);
    long tot_elements = (long)block->num_elements;

    long beg = lmax(beg_elements, off);
    long end = lmin(end_elements, off + tot_elements);
    long num_elements = lmax(0, end - beg);
    assert(tot_elements == sync_lsum(num_elements));
    block->num_elements = (uint64_t)num_elements;

    assert(num_elements <= INT_MAX);
    block->tag = teal_alloc((int)num_elements, sizeof(*block->tag));

    int len = len_node_tags(block);
    assert(num_elements < INT_MAX / len);
    int tot_node_tags = (int)num_elements * len;
    block->node_tag = teal_alloc(tot_node_tags, sizeof(*block->node_tag));

    int stride = 1 + len;
    assert(num_elements <= INT_MAX / stride);
    uint64_t (*tmp)[stride] = teal_alloc((int)num_elements, sizeof(*tmp));
    parse(file, tmp, (int)num_elements * stride, MPI_UINT64_T, mode | SPLIT);
    for (int i = 0; i < num_elements; i++) {
        block->tag[i] = tmp[i][0];
        for (int j = 0; j < len; j++) {
            block->node_tag[(i * len) + j] = tmp[i][1 + j];
        }
    }
    teal_free(tmp);

    return tot_elements;
}

// Read the elements section into gmsh.
static void read_elements(Gmsh *gmsh, Parse *file, int mode)
{
    parse(file, &gmsh->elements.num_blocks, 1, MPI_UINT64_T, mode);
    assert(gmsh->elements.num_blocks <= INT_MAX);
    int num_blocks = (int)gmsh->elements.num_blocks;

    parse(file, &gmsh->elements.num_elements, 1, MPI_UINT64_T, mode);
    assert(gmsh->elements.num_elements <= LONG_MAX);
    long tot_elements = (long)gmsh->elements.num_elements;

    parse(file, &gmsh->elements.min_tag, 1, MPI_UINT64_T, mode);
    parse(file, &gmsh->elements.max_tag, 1, MPI_UINT64_T, mode);

    long base = tot_elements / sync.size;
    long extra = tot_elements % sync.size;
    long beg_elements = (sync.rank * base) + ((sync.rank < extra) ? sync.rank : extra);
    long end_elements = beg_elements + base + (sync.rank < extra);
    long num_elements = end_elements - beg_elements;
    assert(num_elements >= 0);
    gmsh->elements.num_elements = (uint64_t)num_elements;

    ElementBlock *block = teal_alloc(num_blocks, sizeof(*block));
    long off = 0;
    for (int i = 0; i < num_blocks; i++) {
        off += read_element_block(&block[i], beg_elements, end_elements, off, file, mode);
    }
    assert(off == tot_elements);
    gmsh->elements.block = block;
}

// Read one periodic link.
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

// Read the periodic section into gmsh.
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

// Parse a Gmsh file into an in-memory struct.
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

// Compare physicals by dimension then tag.
static int cmp_physical(const void *lhs_, const void *rhs_)
{
    const Physical *lhs = lhs_;
    const Physical *rhs = rhs_;
    if (lhs->dim != rhs->dim) {
        assert(sizeof(((Physical *)0)->dim) == sizeof(int));
        return -cmp_int(&lhs->dim, &rhs->dim);
    }
    assert(sizeof(((Physical *)0)->tag) == sizeof(int));
    return cmp_int(&lhs->tag, &rhs->tag);
}

// Compare points by tag.
static int cmp_point(const void *lhs_, const void *rhs_)
{
    const Point *lhs = lhs_;
    const Point *rhs = rhs_;
    assert(sizeof(((Point *)0)->tag) == sizeof(int));
    return cmp_int(&lhs->tag, &rhs->tag);
}

// Compare curves by tag.
static int cmp_curve(const void *lhs_, const void *rhs_)
{
    const Curve *lhs = lhs_;
    const Curve *rhs = rhs_;
    assert(sizeof(((Curve *)0)->tag) == sizeof(int));
    return cmp_int(&lhs->tag, &rhs->tag);
}

// Compare surfaces by tag.
static int cmp_surface(const void *lhs_, const void *rhs_)
{
    const Surface *lhs = lhs_;
    const Surface *rhs = rhs_;
    assert(sizeof(((Surface *)0)->tag) == sizeof(int));
    return cmp_int(&lhs->tag, &rhs->tag);
}

// Compare volumes by tag.
static int cmp_volume(const void *lhs_, const void *rhs_)
{
    const Volume *lhs = lhs_;
    const Volume *rhs = rhs_;
    assert(sizeof(((Volume *)0)->tag) == sizeof(int));
    return cmp_int(&lhs->tag, &rhs->tag);
}

// Compare node blocks by dimension then tag.
static int cmp_node_block(const void *lhs_, const void *rhs_)
{
    const NodeBlock *lhs = lhs_;
    const NodeBlock *rhs = rhs_;
    if (lhs->entity_dim != rhs->entity_dim) {
        assert(sizeof(((NodeBlock *)0)->entity_dim) == sizeof(int));
        return -cmp_int(&lhs->entity_dim, &rhs->entity_dim);
    }
    assert(sizeof(((NodeBlock *)0)->entity_tag) == sizeof(int));
    return cmp_int(&lhs->entity_tag, &rhs->entity_tag);
}

// Compare element blocks by dimension then tag.
static int cmp_element_block(const void *lhs_, const void *rhs_)
{
    const ElementBlock *lhs = lhs_;
    const ElementBlock *rhs = rhs_;
    if (lhs->entity_dim != rhs->entity_dim) {
        assert(sizeof(((ElementBlock *)0)->entity_dim) == sizeof(int));
        return -cmp_int(&lhs->entity_dim, &rhs->entity_dim);
    }
    assert(sizeof(((ElementBlock *)0)->entity_tag) == sizeof(int));
    return cmp_int(&lhs->entity_tag, &rhs->entity_tag);
}

// Compare periodic links by dimension then tag.
static int cmp_periodic(const void *lhs_, const void *rhs_)
{
    const Periodic *lhs = lhs_;
    const Periodic *rhs = rhs_;
    if (lhs->entity_dim != rhs->entity_dim) {
        assert(sizeof(((Periodic *)0)->entity_dim) == sizeof(int));
        return -cmp_int(&lhs->entity_dim, &rhs->entity_dim);
    }
    assert(sizeof(((Periodic *)0)->entity_tag) == sizeof(int));
    return cmp_int(&lhs->entity_tag, &rhs->entity_tag);
}

// Sort gmsh arrays into a deterministic order.
static void reorder(Gmsh *gmsh)
{
    assert(gmsh->physicals.num >= 0);
    qsort(gmsh->physicals.physical, (size_t)gmsh->physicals.num, sizeof(*gmsh->physicals.physical),
          cmp_physical);

    qsort(gmsh->entities.point, gmsh->entities.num_points, sizeof(*gmsh->entities.point),
          cmp_point);

    qsort(gmsh->entities.curve, gmsh->entities.num_curves, sizeof(*gmsh->entities.curve),
          cmp_curve);

    qsort(gmsh->entities.surface, gmsh->entities.num_surfaces, sizeof(*gmsh->entities.surface),
          cmp_surface);

    qsort(gmsh->entities.volume, gmsh->entities.num_volumes, sizeof(*gmsh->entities.volume),
          cmp_volume);

    qsort(gmsh->nodes.block, gmsh->nodes.num_blocks, sizeof(*gmsh->nodes.block), cmp_node_block);

    qsort(gmsh->elements.block, gmsh->elements.num_blocks, sizeof(*gmsh->elements.block),
          cmp_element_block);

    qsort(gmsh->periodics.periodic, gmsh->periodics.num, sizeof(*gmsh->periodics.periodic),
          cmp_periodic);
}

// Populate mesh nodes from gmsh data.
static void create_nodes(Mesh *mesh, const Gmsh *gmsh)
{
    assert(gmsh->nodes.num_blocks <= INT_MAX);
    int num_blocks = (int)gmsh->nodes.num_blocks;
    const NodeBlock *block = gmsh->nodes.block;

    assert(gmsh->nodes.num_nodes <= INT_MAX);
    mesh->nodes.num = (int)gmsh->nodes.num_nodes;

    vector *coord = teal_alloc(mesh->nodes.num, sizeof(*coord));
    int num = 0;
    for (int i = 0; i < num_blocks; i++) {
        int len = len_coord(&block[i]);
        assert(block[i].num_nodes <= INT_MAX);
        for (int j = 0; j < (int)block[i].num_nodes; j++) {
            coord[num].x = block[i].coord[(j * len) + 0];
            coord[num].y = block[i].coord[(j * len) + 1];
            coord[num].z = block[i].coord[(j * len) + 2];
            num += 1;
        }
    }
    assert(num == mesh->nodes.num);
    mesh->nodes.coord = coord;
}

// Map node tags to local indices.
static long *tag_to_idx(long *tag, int num_tags, const Mesh *mesh, const Gmsh *gmsh)
{
    assert(gmsh->nodes.num_blocks <= INT_MAX);
    int num_blocks = (int)gmsh->nodes.num_blocks;
    const NodeBlock *block = gmsh->nodes.block;

    typedef struct {
        long tag, idx;
    } Map;

    int cap = sync_max(mesh->nodes.num);
    Map *map = teal_alloc(cap, sizeof(*map));

    long prefix = sync_lexsum(mesh->nodes.num);
    int num = 0;
    for (int i = 0; i < num_blocks; i++) {
        assert(block[i].num_nodes <= INT_MAX);
        for (int j = 0; j < (int)block[i].num_nodes; j++) {
            assert(0 < block[i].tag[j] && block[i].tag[j] <= LONG_MAX);
            map[num].tag = (long)block[i].tag[j];
            map[num].idx = prefix + num;
            num += 1;
        }
    }
    assert(num == mesh->nodes.num);
    qsort(map, (size_t)num, sizeof(*map), cmp_long);

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

// Build a node graph for element blocks of one dimension.
static void create_node_graph(Graph *node, int num_items, int max_len, int dim, const Mesh *mesh,
                              const Gmsh *gmsh)
{
    assert(gmsh->elements.num_blocks <= INT_MAX);
    int num_blocks = (int)gmsh->elements.num_blocks;
    const ElementBlock *block = gmsh->elements.block;

    int *off = teal_alloc(num_items + 1, sizeof(*off));
    long *tag = teal_alloc(num_items * max_len, sizeof(*tag));
    int num = 0;
    for (int i = 0; i < num_blocks; i++) {
        if (block[i].entity_dim != dim) {
            continue;
        }
        int len = len_node_tags(&block[i]);
        assert(len <= max_len);
        assert(block[i].num_elements <= INT_MAX);
        for (int j = 0; j < (int)block[i].num_elements; j++) {
            off[num + 1] = off[num] + len;
            for (int k = 0; k < len; k++) {
                assert(0 < block[i].node_tag[(j * len) + k] &&
                       block[i].node_tag[(j * len) + k] <= LONG_MAX);
                tag[off[num] + k] = (long)block[i].node_tag[(j * len) + k];
            }
            num += 1;
        }
    }
    assert(num == num_items);
    node->off = off;
    node->idx = tag_to_idx(tag, off[num], mesh, gmsh);
}

// Build inner cell connectivity from element blocks.
static void create_inner_cells(Mesh *mesh, const Gmsh *gmsh)
{
    assert(gmsh->elements.num_blocks <= INT_MAX);
    int num_blocks = (int)gmsh->elements.num_blocks;
    const ElementBlock *block = gmsh->elements.block;

    int num_cells = 0;
    for (int i = 0; i < num_blocks; i++) {
        if (block[i].entity_dim != 3) {
            continue;
        }
        assert(block[i].num_elements <= INT_MAX);
        assert(num_cells <= INT_MAX - (int)block[i].num_elements);
        num_cells += (int)block[i].num_elements;
    }
    assert(num_cells > 0);
    mesh->cells.num = num_cells;

    create_node_graph(&mesh->cells.node, num_cells, MAX_CELL_NODES, 3, mesh, gmsh);
}

// Build boundary face connectivity from element blocks.
static void create_boundary_faces(Mesh *mesh, const Gmsh *gmsh)
{
    assert(gmsh->elements.num_blocks <= INT_MAX);
    int num_blocks = (int)gmsh->elements.num_blocks;
    const ElementBlock *block = gmsh->elements.block;

    int num_faces = 0;
    for (int i = 0; i < num_blocks; i++) {
        if (block[i].entity_dim != 2) {
            continue;
        }
        assert(block[i].num_elements <= INT_MAX);
        assert(num_faces <= INT_MAX - (int)block[i].num_elements);
        num_faces += (int)block[i].num_elements;
    }
    assert(num_faces >= 0);
    mesh->faces.num = num_faces;

    create_node_graph(&mesh->faces.node, num_faces, MAX_FACE_NODES, 2, mesh, gmsh);
}

// Return the physical index for a dimension and tag.
static int physical_index(int dim, int tag, const Gmsh *gmsh)
{
    Physical key = {.dim = dim, .tag = tag};
    assert(gmsh->physicals.num >= 0);
    const Physical *val = bsearch(&key, gmsh->physicals.physical, (size_t)gmsh->physicals.num,
                                  sizeof(*gmsh->physicals.physical), cmp_physical);
    if (!val) {
        teal_error("could not find physical (%d, %d)", dim, tag);
    }
    ptrdiff_t idx = val - gmsh->physicals.physical;
    assert(idx <= INT_MAX);
    return (int)idx;
}

// Find a volume entity by tag.
static const Volume *find_volume(int tag, const Gmsh *gmsh)
{
    Volume key = {.tag = tag};
    const Volume *val = bsearch(&key, gmsh->entities.volume, gmsh->entities.num_volumes,
                                sizeof(*gmsh->entities.volume), cmp_volume);
    if (!val) {
        teal_error("could not find volume entity (%d)", tag);
    }
    return val;
}

// Find a surface entity by tag.
static const Surface *find_surface(int tag, const Gmsh *gmsh)
{
    Surface key = {.tag = tag};
    const Surface *val = bsearch(&key, gmsh->entities.surface, gmsh->entities.num_surfaces,
                                 sizeof(*gmsh->entities.surface), cmp_surface);
    if (!val) {
        teal_error("could not find surface entity (%d)", tag);
    }
    return val;
}

// Return the physical entity index for a mesh entity tag.
static int entity_index(int dim, int tag, const Gmsh *gmsh)
{
    switch (dim) {
        case 3: {
            const Volume *volume = find_volume(tag, gmsh);
            if (volume->num_physical_tags != 1) {
                teal_error("unsupported number of physical tags (%" PRIu64 ")",
                           volume->num_physical_tags);
            }
            return physical_index(dim, volume->physical_tag[0], gmsh);
        }
        case 2: {
            const Surface *surface = find_surface(tag, gmsh);
            if (surface->num_physical_tags != 1) {
                teal_error("unsupported number of physical tags (%" PRIu64 ")",
                           surface->num_physical_tags);
            }
            return physical_index(dim, surface->physical_tag[0], gmsh);
        }
        default: teal_error("unsupported entity dimension (%d)", dim);
    }
}

// Create mesh entities from physical names.
static void create_entities(Mesh *mesh, const Gmsh *gmsh)
{
    int num_entities = gmsh->physicals.num;
    assert(num_entities > 0);
    mesh->entities.num = num_entities;

    const Physical *physical = gmsh->physicals.physical;

    Name *name = teal_alloc(num_entities, sizeof(*name));
    int num_inner = 0;
    int num_boundary = 0;
    for (int i = 0; i < num_entities; i++) {
        strcpy(name[i], physical[i].name);
        switch (physical[i].dim) {
            case 3: num_inner += 1; break;
            case 2: num_boundary += 1; break;
            default: teal_error("unsupported physical dimension (%d)", physical[i].dim);
        }
    }
    assert(num_boundary == num_entities - num_inner);
    mesh->entities.num_inner = num_inner;
    mesh->entities.name = name;

    assert(gmsh->elements.num_blocks <= INT_MAX);
    int num_blocks = (int)gmsh->elements.num_blocks;
    const ElementBlock *block = gmsh->elements.block;

    int *num_volumes = teal_alloc(num_entities, sizeof(*num_volumes));
    int *num_surfaces = teal_alloc(num_entities, sizeof(*num_surfaces));
    for (int i = 0; i < num_blocks; i++) {
        int idx = entity_index(block[i].entity_dim, block[i].entity_tag, gmsh);
        switch (block[i].entity_dim) {
            case 3:
                assert(block[i].num_elements <= INT_MAX);
                assert(num_volumes[idx] <= INT_MAX - (int)block[i].num_elements);
                num_volumes[idx] += (int)block[i].num_elements;
                break;
            case 2:
                assert(block[i].num_elements <= INT_MAX);
                assert(num_surfaces[idx] <= INT_MAX - (int)block[i].num_elements);
                num_surfaces[idx] += (int)block[i].num_elements;
                break;
            default: teal_error("unsupported entity dimension (%d)", block[i].entity_dim);
        }
    }

    int *cell_off = teal_alloc(num_entities + 1, sizeof(*cell_off));
    int *face_off = teal_alloc(num_entities + 1, sizeof(*face_off));
    for (int i = 0; i < num_entities; i++) {
        cell_off[i + 1] = cell_off[i] + num_volumes[i];
        face_off[i + 1] = face_off[i] + num_surfaces[i];
    }
    mesh->entities.cell_off = cell_off;
    mesh->entities.face_off = face_off;

    teal_free(num_volumes);
    teal_free(num_surfaces);
}

// Create mesh periodic data from gmsh.
static void create_periodics(Mesh *mesh, const Gmsh *gmsh)
{
}

// Release all gmsh allocations.
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

// Read a Gmsh file into a mesh.
static void read_gmsh(Mesh *mesh, const char *fname)
{
    Gmsh *gmsh = gmsh_init(fname);

    reorder(gmsh);

    create_nodes(mesh, gmsh);
    create_inner_cells(mesh, gmsh);
    create_boundary_faces(mesh, gmsh);
    create_entities(mesh, gmsh);
    create_periodics(mesh, gmsh);

    gmsh_deinit(gmsh);
}
