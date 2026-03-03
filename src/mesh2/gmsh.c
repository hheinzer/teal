#include <assert.h>
#include <inttypes.h>
#include <limits.h>
#include <stdint.h>
#include <string.h>

#include "arena2.h"
#include "grid.h"
#include "mesh2.h"
#include "parse2.h"
#include "sync2.h"
#include "teal2.h"
#include "utils2.h"

static Arena *arena;

typedef struct {
    double version;
    int32_t file_type;
    int32_t data_size;
} GmshFormat;

typedef struct {
    int32_t dim;
    int32_t tag;
    char name[128];
} Physical;

typedef struct {
    int32_t num;
    Physical *physical;
} GmshPhysicals;

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
} GmshEntities;

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
} GmshNodes;

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
} GmshElements;

typedef struct {
    int32_t entity_dim;
    int32_t entity_tag;
    int32_t entity_tag_master;
    uint64_t num_affine;
    double *affine;
    uint64_t num_nodes;
    uint64_t *node_tag;
    uint64_t *node_tag_master;
} Link;

typedef struct {
    uint64_t num_links;
    Link *link;
} GmshPeriodic;

typedef struct {
    GmshFormat format;
    GmshPhysicals physicals;
    GmshEntities entities;
    GmshNodes nodes;
    GmshElements elements;
    GmshPeriodic periodic;
} Gmsh;

static int read_format(Gmsh *gmsh, MPI_File file)
{
    parse2_ascii(file, &gmsh->format.version, 1, MPI_DOUBLE);
    if (!isclose(gmsh->format.version, 4.1)) {
        teal2_error("unsupported version (%g)", gmsh->format.version);
    }

    parse2_ascii(file, &gmsh->format.file_type, 1, MPI_INT32_T);

    parse2_ascii(file, &gmsh->format.data_size, 1, MPI_INT32_T);
    if (gmsh->format.data_size != 8) {
        teal2_error("unsupported data-size (%d)", gmsh->format.data_size);
    }

    int mode = 0;
    if (gmsh->format.file_type == 1) {
        int32_t one;
        parse2_binary(file, &one, 1, MPI_INT32_T, 0);
        mode = BINARY | ((one != 1) ? SWAP : 0);
    }
    return mode;
}

static void read_physicals(Gmsh *gmsh, MPI_File file)
{
    parse2_ascii(file, &gmsh->physicals.num, 1, MPI_INT32_T);

    Physical *physical = arena2_calloc(arena, gmsh->physicals.num, sizeof(*physical));
    for (int32_t i = 0; i < gmsh->physicals.num; i++) {
        parse2_ascii(file, &physical[i].dim, 1, MPI_INT32_T);
        parse2_ascii(file, &physical[i].tag, 1, MPI_INT32_T);
        parse2_string(file, physical[i].name, sizeof(physical[i].name));
    }
    gmsh->physicals.physical = physical;
}

static void read_points(Gmsh *gmsh, MPI_File file, int mode)
{
    assert(gmsh->entities.num_points <= INT_MAX);
    Point *point = arena2_calloc(arena, (int)gmsh->entities.num_points, sizeof(*point));
    for (uint64_t i = 0; i < gmsh->entities.num_points; i++) {
        parse2(file, &point[i].tag, 1, MPI_INT32_T, mode);
        parse2(file, &point[i].x, 1, MPI_DOUBLE, mode);
        parse2(file, &point[i].y, 1, MPI_DOUBLE, mode);
        parse2(file, &point[i].z, 1, MPI_DOUBLE, mode);

        parse2(file, &point[i].num_physical_tags, 1, MPI_UINT64_T, mode);
        assert(point[i].num_physical_tags <= INT_MAX);
        int num_physical_tags = (int)point[i].num_physical_tags;

        int32_t *physical_tag = arena2_calloc(arena, num_physical_tags, sizeof(*physical_tag));
        parse2(file, physical_tag, num_physical_tags, MPI_INT32_T, mode);
        point[i].physical_tag = physical_tag;
    }
    gmsh->entities.point = point;
}

static void read_curves(Gmsh *gmsh, MPI_File file, int mode)
{
    assert(gmsh->entities.num_curves <= INT_MAX);
    Curve *curve = arena2_calloc(arena, (int)gmsh->entities.num_curves, sizeof(*curve));
    for (uint64_t i = 0; i < gmsh->entities.num_curves; i++) {
        parse2(file, &curve[i].tag, 1, MPI_INT32_T, mode);
        parse2(file, &curve[i].min_x, 1, MPI_DOUBLE, mode);
        parse2(file, &curve[i].min_y, 1, MPI_DOUBLE, mode);
        parse2(file, &curve[i].min_z, 1, MPI_DOUBLE, mode);
        parse2(file, &curve[i].max_x, 1, MPI_DOUBLE, mode);
        parse2(file, &curve[i].max_y, 1, MPI_DOUBLE, mode);
        parse2(file, &curve[i].max_z, 1, MPI_DOUBLE, mode);

        parse2(file, &curve[i].num_physical_tags, 1, MPI_UINT64_T, mode);
        assert(curve[i].num_physical_tags <= INT_MAX);
        int num_physical_tags = (int)curve[i].num_physical_tags;

        int32_t *physical_tag = arena2_calloc(arena, num_physical_tags, sizeof(*physical_tag));
        parse2(file, physical_tag, num_physical_tags, MPI_INT32_T, mode);
        curve[i].physical_tag = physical_tag;

        parse2(file, &curve[i].num_point_tags, 1, MPI_UINT64_T, mode);
        assert(curve[i].num_point_tags <= INT_MAX);
        int num_point_tags = (int)curve[i].num_point_tags;

        int32_t *point_tag = arena2_calloc(arena, num_point_tags, sizeof(*point_tag));
        parse2(file, point_tag, num_point_tags, MPI_INT32_T, mode);
        curve[i].point_tag = point_tag;
    }
    gmsh->entities.curve = curve;
}

static void read_surfaces(Gmsh *gmsh, MPI_File file, int mode)
{
    assert(gmsh->entities.num_surfaces <= INT_MAX);
    Surface *surface = arena2_calloc(arena, (int)gmsh->entities.num_surfaces, sizeof(*surface));
    for (uint64_t i = 0; i < gmsh->entities.num_surfaces; i++) {
        parse2(file, &surface[i].tag, 1, MPI_INT32_T, mode);
        parse2(file, &surface[i].min_x, 1, MPI_DOUBLE, mode);
        parse2(file, &surface[i].min_y, 1, MPI_DOUBLE, mode);
        parse2(file, &surface[i].min_z, 1, MPI_DOUBLE, mode);
        parse2(file, &surface[i].max_x, 1, MPI_DOUBLE, mode);
        parse2(file, &surface[i].max_y, 1, MPI_DOUBLE, mode);
        parse2(file, &surface[i].max_z, 1, MPI_DOUBLE, mode);

        parse2(file, &surface[i].num_physical_tags, 1, MPI_UINT64_T, mode);
        assert(surface[i].num_physical_tags <= INT_MAX);
        int num_physical_tags = (int)surface[i].num_physical_tags;

        int32_t *physical_tag = arena2_calloc(arena, num_physical_tags, sizeof(*physical_tag));
        parse2(file, physical_tag, num_physical_tags, MPI_INT32_T, mode);
        surface[i].physical_tag = physical_tag;

        parse2(file, &surface[i].num_curve_tags, 1, MPI_UINT64_T, mode);
        assert(surface[i].num_curve_tags <= INT_MAX);
        int num_curve_tags = (int)surface[i].num_curve_tags;

        int32_t *curve_tag = arena2_calloc(arena, num_curve_tags, sizeof(*curve_tag));
        parse2(file, curve_tag, num_curve_tags, MPI_INT32_T, mode);
        surface[i].curve_tag = curve_tag;
    }
    gmsh->entities.surface = surface;
}

static void read_volumes(Gmsh *gmsh, MPI_File file, int mode)
{
    assert(gmsh->entities.num_volumes <= INT_MAX);
    Volume *volume = arena2_calloc(arena, (int)gmsh->entities.num_volumes, sizeof(*volume));
    for (uint64_t i = 0; i < gmsh->entities.num_volumes; i++) {
        parse2(file, &volume[i].tag, 1, MPI_INT32_T, mode);
        parse2(file, &volume[i].min_x, 1, MPI_DOUBLE, mode);
        parse2(file, &volume[i].min_y, 1, MPI_DOUBLE, mode);
        parse2(file, &volume[i].min_z, 1, MPI_DOUBLE, mode);
        parse2(file, &volume[i].max_x, 1, MPI_DOUBLE, mode);
        parse2(file, &volume[i].max_y, 1, MPI_DOUBLE, mode);
        parse2(file, &volume[i].max_z, 1, MPI_DOUBLE, mode);

        parse2(file, &volume[i].num_physical_tags, 1, MPI_UINT64_T, mode);
        assert(volume[i].num_physical_tags <= INT_MAX);
        int num_physical_tags = (int)volume[i].num_physical_tags;

        int32_t *physical_tag = arena2_calloc(arena, num_physical_tags, sizeof(*physical_tag));
        parse2(file, physical_tag, num_physical_tags, MPI_INT32_T, mode);
        volume[i].physical_tag = physical_tag;

        parse2(file, &volume[i].num_surface_tags, 1, MPI_UINT64_T, mode);
        assert(volume[i].num_surface_tags <= INT_MAX);
        int num_surface_tags = (int)volume[i].num_surface_tags;

        int32_t *surface_tag = arena2_calloc(arena, num_surface_tags, sizeof(*surface_tag));
        parse2(file, surface_tag, num_surface_tags, MPI_INT32_T, mode);
        volume[i].surface_tag = surface_tag;
    }
    gmsh->entities.volume = volume;
}

static void read_entities(Gmsh *gmsh, MPI_File file, int mode)
{
    parse2(file, &gmsh->entities.num_points, 1, MPI_UINT64_T, mode);
    parse2(file, &gmsh->entities.num_curves, 1, MPI_UINT64_T, mode);
    parse2(file, &gmsh->entities.num_surfaces, 1, MPI_UINT64_T, mode);
    parse2(file, &gmsh->entities.num_volumes, 1, MPI_UINT64_T, mode);

    read_points(gmsh, file, mode);
    read_curves(gmsh, file, mode);
    read_surfaces(gmsh, file, mode);
    read_volumes(gmsh, file, mode);
}

typedef struct {
    uint64_t beg, end;
} Slice;

static Slice slice_rank(uint64_t total)
{
    assert(sync2.rank >= 0 && sync2.size >= 0);
    uint64_t rank = (uint64_t)sync2.rank;
    uint64_t size = (uint64_t)sync2.size;
    uint64_t base = total / size;
    uint64_t extra = total % size;
    uint64_t beg = (rank * base) + ((rank < extra) ? rank : extra);
    uint64_t end = beg + base + (rank < extra);
    return (Slice){beg, end};
}

static int num_overlap(uint64_t beg, uint64_t end, uint64_t offset, uint64_t tot)
{
    uint64_t min = (end < offset + tot) ? end : offset + tot;
    uint64_t max = (beg > offset) ? beg : offset;
    uint64_t num = (min > max) ? min - max : 0;
    assert(num <= INT_MAX);
    return (int)num;
}

static int len_coord(const NodeBlock *block)
{
    return 3 + (block->parametric ? block->entity_dim : 0);
}

static void read_nodes(Gmsh *gmsh, MPI_File file, int mode)
{
    parse2(file, &gmsh->nodes.num_blocks, 1, MPI_UINT64_T, mode);
    parse2(file, &gmsh->nodes.num_nodes, 1, MPI_UINT64_T, mode);
    parse2(file, &gmsh->nodes.min_tag, 1, MPI_UINT64_T, mode);
    parse2(file, &gmsh->nodes.max_tag, 1, MPI_UINT64_T, mode);

    uint64_t total = gmsh->nodes.num_nodes;
    Slice slice = slice_rank(total);
    gmsh->nodes.num_nodes = slice.end - slice.beg;

    assert(gmsh->nodes.num_blocks <= INT_MAX);
    NodeBlock *block = arena2_calloc(arena, (int)gmsh->nodes.num_blocks, sizeof(*block));
    uint64_t offset = 0;
    for (uint64_t i = 0; i < gmsh->nodes.num_blocks; i++) {
        parse2(file, &block[i].entity_dim, 1, MPI_INT32_T, mode);
        parse2(file, &block[i].entity_tag, 1, MPI_INT32_T, mode);
        parse2(file, &block[i].parametric, 1, MPI_INT32_T, mode);
        parse2(file, &block[i].num_nodes, 1, MPI_UINT64_T, mode);

        uint64_t tot = block[i].num_nodes;
        int num = num_overlap(slice.beg, slice.end, offset, tot);
        block[i].num_nodes = num;

        block[i].tag = arena2_calloc(arena, num, sizeof(*block[i].tag));
        parse2(file, block[i].tag, num, MPI_UINT64_T, mode | SPLIT);

        int len = len_coord(&block[i]);
        block[i].coord = arena2_calloc(arena, num * len, sizeof(*block[i].coord));
        parse2(file, block[i].coord, num * len, MPI_DOUBLE, mode | SPLIT);

        offset += tot;
    }
    assert(offset == total);

    gmsh->nodes.block = block;
}

static int len_node_tags(const ElementBlock *block)
{
    switch (block->element_type) {
        case 2: return 3;  // triangle
        case 3:            // quadrangle
        case 4: return 4;  // tetrahedron
        case 5: return 8;  // hexahedron
        case 6: return 6;  // prism
        case 7: return 5;  // pyramid
        default: teal2_error("unsupported element type (%d)", block->element_type);
    }
}

static void read_elements(Gmsh *gmsh, MPI_File file, int mode)
{
    parse2(file, &gmsh->elements.num_blocks, 1, MPI_UINT64_T, mode);
    parse2(file, &gmsh->elements.num_elements, 1, MPI_UINT64_T, mode);
    parse2(file, &gmsh->elements.min_tag, 1, MPI_UINT64_T, mode);
    parse2(file, &gmsh->elements.max_tag, 1, MPI_UINT64_T, mode);

    uint64_t total = gmsh->elements.num_elements;
    Slice slice = slice_rank(total);
    gmsh->elements.num_elements = slice.end - slice.beg;

    assert(gmsh->elements.num_blocks <= INT_MAX);
    ElementBlock *block = arena2_calloc(arena, (int)gmsh->elements.num_blocks, sizeof(*block));
    uint64_t offset = 0;
    for (uint64_t i = 0; i < gmsh->elements.num_blocks; i++) {
        parse2(file, &block[i].entity_dim, 1, MPI_INT32_T, mode);
        parse2(file, &block[i].entity_tag, 1, MPI_INT32_T, mode);
        parse2(file, &block[i].element_type, 1, MPI_INT32_T, mode);
        parse2(file, &block[i].num_elements, 1, MPI_UINT64_T, mode);

        uint64_t tot = block[i].num_elements;
        int num = num_overlap(slice.beg, slice.end, offset, tot);
        block[i].num_elements = num;

        block[i].tag = arena2_calloc(arena, num, sizeof(*block[i].tag));

        int len = len_node_tags(&block[i]);
        block[i].node_tag = arena2_calloc(arena, num * len, sizeof(*block[i].node_tag));

        Save *save = arena2_save(arena);

        uint64_t (*buf)[1 + len] = arena2_calloc(arena, num, (int)sizeof(*buf));
        parse2(file, buf, num * (1 + len), MPI_UINT64_T, mode | SPLIT);

        for (int j = 0; j < num; j++) {
            block[i].tag[j] = buf[j][0];
            for (int k = 0; k < len; k++) {
                block[i].node_tag[(j * len) + k] = buf[j][1 + k];
            }
        }

        arena2_load(arena, save);
        offset += tot;
    }
    assert(offset == total);

    gmsh->elements.block = block;
}

static void read_periodic(Gmsh *gmsh, MPI_File file, int mode)
{
    parse2(file, &gmsh->periodic.num_links, 1, MPI_UINT64_T, mode);

    assert(gmsh->periodic.num_links <= INT_MAX);
    Link *link = arena2_calloc(arena, (int)gmsh->periodic.num_links, sizeof(*link));
    for (uint64_t i = 0; i < gmsh->periodic.num_links; i++) {
        parse2(file, &link[i].entity_dim, 1, MPI_INT32_T, mode);
        parse2(file, &link[i].entity_tag, 1, MPI_INT32_T, mode);
        parse2(file, &link[i].entity_tag_master, 1, MPI_INT32_T, mode);
        parse2(file, &link[i].num_affine, 1, MPI_UINT64_T, mode);

        assert(link[i].num_affine <= INT_MAX);
        link[i].affine = arena2_calloc(arena, (int)link[i].num_affine, sizeof(*link[i].affine));
        parse2(file, link[i].affine, (int)link[i].num_affine, MPI_DOUBLE, mode);

        parse2(file, &link[i].num_nodes, 1, MPI_UINT64_T, mode);
        uint64_t total = link[i].num_nodes;
        uint64_t base = total / sync2.size;
        uint64_t extra = total % sync2.size;

        assert(base <= INT_MAX && sync2.rank >= 0);
        int num = (int)base + ((uint64_t)sync2.rank < extra);
        link[i].num_nodes = num;

        link[i].node_tag = arena2_calloc(arena, num, sizeof(*link[i].node_tag));
        link[i].node_tag_master = arena2_calloc(arena, num, sizeof(*link[i].node_tag_master));

        Save *save = arena2_save(arena);

        uint64_t (*buf)[2] = arena2_calloc(arena, num, sizeof(*buf));
        parse2(file, buf, num * 2, MPI_UINT64_T, mode | SPLIT);

        for (int j = 0; j < num; j++) {
            link[i].node_tag[j] = buf[j][0];
            link[i].node_tag_master[j] = buf[j][1];
        }

        arena2_load(arena, save);
    }
    gmsh->periodic.link = link;
}

static Gmsh *read(const char *fname)
{
    MPI_File file = parse2_open(fname);

    Gmsh *gmsh = arena2_calloc(arena, 1, sizeof(*gmsh));
    int mode = 0;
    char section[128];
    while (parse2_string(file, section, sizeof(section))) {
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
            read_periodic(gmsh, file, mode);
        }
        else if (!strncmp(section, "$End", 4)) {
            continue;
        }
        else {
            teal2_error("unsupported mesh section (%s)", section);
        }
    }

    parse2_close(file);
    return gmsh;
}

static int32_t *get_physical_tags(const ElementBlock *block, const Gmsh *gmsh,
                                  uint64_t *num_physical_tags)
{
    switch (block->entity_dim) {
        case 2:
            for (uint64_t i = 0; i < gmsh->entities.num_surfaces; i++) {
                if (gmsh->entities.surface[i].tag == block->entity_tag) {
                    *num_physical_tags = gmsh->entities.surface[i].num_physical_tags;
                    return gmsh->entities.surface[i].physical_tag;
                }
            }
            teal2_error("could not find surface entity (%d)", block->entity_tag);
        case 3:
            for (uint64_t i = 0; i < gmsh->entities.num_volumes; i++) {
                if (gmsh->entities.volume[i].tag == block->entity_tag) {
                    *num_physical_tags = gmsh->entities.volume[i].num_physical_tags;
                    return gmsh->entities.volume[i].physical_tag;
                }
            }
            teal2_error("could not find volume entity (%d)", block->entity_tag);
        default: teal2_error("unsupported element block dimension (%d)", block->entity_dim);
    }
}

static Physical *get_physical(const ElementBlock *block, const Gmsh *gmsh)
{
    uint64_t num_physical_tags;
    int32_t *physical_tag = get_physical_tags(block, gmsh, &num_physical_tags);
    if (num_physical_tags != 1) {
        teal2_error("unsupported number of physical tags (%" PRIu64 ")", num_physical_tags);
    }
    for (int32_t i = 0; i < gmsh->physicals.num; i++) {
        if (gmsh->physicals.physical[i].dim == block->entity_dim &&
            gmsh->physicals.physical[i].tag == physical_tag[0]) {
            return &gmsh->physicals.physical[i];
        }
    }
    teal2_error("could not find physical name (%d, %d)", block->entity_dim, physical_tag[0]);
}

static Link *get_link(const ElementBlock *block, const Gmsh *gmsh)
{
    for (uint64_t i = 0; i < gmsh->periodic.num_links; i++) {
        if (gmsh->periodic.link[i].entity_dim != block->entity_dim) {
            continue;
        }
        if (gmsh->periodic.link[i].entity_tag == block->entity_tag) {
            return &gmsh->periodic.link[i];
        }
        if (gmsh->periodic.link[i].entity_tag_master == block->entity_tag) {
            return &gmsh->periodic.link[i];
        }
    }
    return 0;
}

typedef struct {
    ElementBlock *block;
    Physical *physical;
    Link *link;
} Entity;

static int compare_entity(const void *lhs_, const void *rhs_)
{
    const Entity *lhs = lhs_;
    const Entity *rhs = rhs_;
    if (!lhs->link != !rhs->link) {
        return (!lhs->link && rhs->link) ? -1 : +1;
    }
    if (lhs->physical->dim != rhs->physical->dim) {
        return (lhs->physical->dim > rhs->physical->dim) ? -1 : +1;
    }
    if (lhs->physical->tag != rhs->physical->tag) {
        return (lhs->physical->tag < rhs->physical->tag) ? -1 : +1;
    }
    if (lhs->block->entity_tag != rhs->block->entity_tag) {
        return (lhs->block->entity_tag < rhs->block->entity_tag) ? -1 : +1;
    }
    return 0;
}

static Entity *extract(const Gmsh *gmsh)
{
    assert(gmsh->elements.num_blocks <= INT_MAX);
    int num_blocks = (int)gmsh->elements.num_blocks;

    Entity *entity = arena2_calloc(arena, num_blocks, sizeof(*entity));
    for (uint64_t i = 0; i < gmsh->elements.num_blocks; i++) {
        entity[i].block = &gmsh->elements.block[i];
        entity[i].physical = get_physical(&gmsh->elements.block[i], gmsh);
        entity[i].link = get_link(&gmsh->elements.block[i], gmsh);
    }

    sort(entity, num_blocks, sizeof(*entity), compare_entity);
    return entity;
}

static void create_nodes(Grid *grid, const Gmsh *gmsh)
{
    if (gmsh->nodes.num_blocks == 0 || gmsh->nodes.num_nodes == 0) {
        teal2_error("mesh does not contain any nodes");
    }

    assert(gmsh->nodes.num_nodes <= INT_MAX);
    int num_nodes = (int)gmsh->nodes.num_nodes;

    long *tag = teal2_calloc(num_nodes, sizeof(*tag));
    Vector *coord = teal2_calloc(num_nodes, sizeof(*coord));

    int num = 0;
    for (uint64_t i = 0; i < gmsh->nodes.num_blocks; i++) {
        int len = len_coord(&gmsh->nodes.block[i]);
        for (uint64_t j = 0; j < gmsh->nodes.block[i].num_nodes; j++) {
            assert(gmsh->nodes.block[i].tag[j] <= LONG_MAX);
            tag[num] = (long)gmsh->nodes.block[i].tag[j];
            coord[num].x = gmsh->nodes.block[i].coord[(j * len) + 0];
            coord[num].y = gmsh->nodes.block[i].coord[(j * len) + 1];
            coord[num].z = gmsh->nodes.block[i].coord[(j * len) + 2];
            num += 1;
        }
    }
    assert(num == num_nodes);

    grid->nodes.num = num_nodes;
    grid->nodes.tag = tag;
    grid->nodes.coord = coord;
}

static void create_cells(Grid *grid, const Entity *entity, const Gmsh *gmsh)
{
    if (gmsh->elements.num_blocks == 0 || gmsh->elements.num_elements == 0) {
        teal2_error("mesh does not contain any elements");
    }

    assert(gmsh->elements.num_elements <= INT_MAX);
    int num_cells = (int)gmsh->elements.num_elements;

    long *node_off = teal2_calloc(num_cells + 1, sizeof(*node_off));
    long *node_tag = teal2_calloc(num_cells * MAX_CELL_NODES, sizeof(*node_tag));

    int num = 0;
    for (uint64_t i = 0; i < gmsh->elements.num_blocks; i++) {
        int len = len_node_tags(entity[i].block);
        for (uint64_t j = 0; j < entity[i].block->num_elements; j++) {
            node_off[num + 1] = node_off[num];
            for (int k = 0; k < len; k++) {
                assert(node_off[num + 1] - node_off[num] < MAX_CELL_NODES);
                assert(entity[i].block->node_tag[(j * len) + k] <= LONG_MAX);
                node_tag[node_off[num + 1]++] = (long)entity[i].block->node_tag[(j * len) + k];
            }
            num += 1;
        }
    }
    assert(num == num_cells);

    assert(node_off[num_cells] <= INT_MAX);
    int num_tags = (int)node_off[num_cells];

    grid->cells.num = num_cells;
    grid->cells.node_off = node_off;
    grid->cells.node_idx = grid_tag_to_idx(grid, node_tag, num_tags);

    teal2_free(node_tag);
}

static int count_entities(Grid *grid, const Entity *entity, const Gmsh *gmsh)
{
    if (entity[0].physical->dim != 3) {
        teal2_error("mesh does not contain any volumes");
    }

    int num_inner = 1;
    int num_boundary = 0;
    int num_periodic = 0;
    int num_nodes = 0;
    for (uint64_t i = 1; i < gmsh->elements.num_blocks; i++) {
        if (!entity[i].link) {
            if (entity[i].physical != entity[i - 1].physical) {
                if (entity[i].physical->dim == 3) {
                    num_inner += 1;
                }
                else if (entity[i].physical->dim == 2) {
                    num_boundary += 1;
                }
                else {
                    teal2_error("unsupported physical dimension (%d)", entity[i].physical->dim);
                }
            }
        }
        else {
            if (entity[i].link->entity_dim != 2) {
                teal2_error("unsupported periodic dimension (%d)", entity[i].physical->dim);
            }
            num_periodic += 1;
            assert(entity[i].link->num_nodes <= INT_MAX);
            num_nodes += (int)entity[i].link->num_nodes;
        }
    }

    grid->entities.num = num_inner + num_boundary + num_periodic;
    grid->entities.num_inner = num_inner;
    grid->entities.off_boundary = num_inner + num_boundary;
    return num_nodes;
}

static int get_periodic(const Entity *key, const Entity *entity, const Gmsh *gmsh)
{
    int entity_tag = (key->link->entity_tag == key->block->entity_tag)
                         ? key->link->entity_tag_master
                         : key->link->entity_tag;
    int num = 1;
    for (uint64_t i = 1; i < gmsh->elements.num_blocks; i++) {
        if (!entity[i].link) {
            if (entity[i].physical != entity[i - 1].physical) {
                num += 1;
            }
        }
        else {
            if (entity[i].block->entity_tag == entity_tag) {
                return num;
            }
            num += 1;
        }
    }
    teal2_error("could not find periodic link (%d)", entity_tag);
}

static void get_affine_transformation(Matrix *rotation, Vector *translation, const Entity *entity)
{
    if (entity->link->num_affine != 16) {
        teal2_error("unsupported affine transformation size (%" PRIu64 ")",
                    entity->link->num_affine);
    }

    int len = 4;

    rotation->x.x = entity->link->affine[(0 * len) + 0];
    rotation->x.y = entity->link->affine[(0 * len) + 1];
    rotation->x.z = entity->link->affine[(0 * len) + 2];

    rotation->y.x = entity->link->affine[(1 * len) + 0];
    rotation->y.y = entity->link->affine[(1 * len) + 1];
    rotation->y.z = entity->link->affine[(1 * len) + 2];

    rotation->z.x = entity->link->affine[(2 * len) + 0];
    rotation->z.y = entity->link->affine[(2 * len) + 1];
    rotation->z.z = entity->link->affine[(2 * len) + 2];

    translation->x = entity->link->affine[(0 * len) + 3];
    translation->y = entity->link->affine[(1 * len) + 3];
    translation->z = entity->link->affine[(2 * len) + 3];

    if (!isclose(matrix_determinant(*rotation), 1)) {
        teal2_error("invalid affine transformation determinant (%g)",
                    matrix_determinant(*rotation));
    }

    if (entity->block->entity_tag == entity->link->entity_tag_master) {
        *rotation = matrix_inverse(*rotation);
        *translation = vector2_mul(-1, matrix_vector(*rotation, *translation));
    }
}

static void get_node_mapping(long *node_off, long *node_tag, const Entity *entity)
{
    assert(entity->link->num_nodes <= LONG_MAX);
    *node_off += (long)entity->link->num_nodes;
    if (entity->block->entity_tag == entity->link->entity_tag) {
        for (uint64_t i = 0; i < entity->link->num_nodes; i++) {
            assert(entity->link->node_tag[i] <= LONG_MAX);
            node_tag[i] = (long)entity->link->node_tag[i];
        }
    }
    else {
        for (uint64_t i = 0; i < entity->link->num_nodes; i++) {
            assert(entity->link->node_tag_master[i] <= LONG_MAX);
            node_tag[i] = (long)entity->link->node_tag_master[i];
        }
    }
}

static void create_entities(Grid *grid, const Entity *entity, const Gmsh *gmsh)
{
    assert(gmsh->elements.num_blocks <= INT_MAX);
    int num_blocks = (int)gmsh->elements.num_blocks;

    int num_nodes = count_entities(grid, entity, gmsh);

    String *name = teal2_calloc(grid->entities.num, sizeof(*name));
    int *cell_off = teal2_calloc(grid->entities.num + 1, sizeof(*cell_off));
    int *periodic = teal2_calloc(grid->entities.num, sizeof(*periodic));
    Matrix *rotation = teal2_calloc(grid->entities.num, sizeof(*rotation));
    Vector *translation = teal2_calloc(grid->entities.num, sizeof(*translation));
    long *node_off = teal2_calloc(grid->entities.num + 1, sizeof(*node_off));
    long *node_tag = teal2_calloc(num_nodes, sizeof(*node_tag));

    int num = 0;
    for (int i = 0; i < grid->entities.num; i++) {
        cell_off[i + 1] = cell_off[i];
        node_off[i + 1] = node_off[i];
        strcpy(name[i], entity[num].physical->name);
        if (!entity[num].link) {
            do {
                assert(entity[num].block->num_elements <= INT_MAX);
                cell_off[i + 1] += (int)entity[num].block->num_elements;
                num += 1;
            } while (num < num_blocks && entity[num].physical == entity[num - 1].physical);
        }
        else {
            assert(entity[num].block->num_elements <= INT_MAX);
            cell_off[i + 1] += (int)entity[num].block->num_elements;
            periodic[i] = get_periodic(&entity[num], entity, gmsh);
            get_affine_transformation(&rotation[i], &translation[i], &entity[num]);
            get_node_mapping(&node_off[i + 1], &node_tag[node_off[i]], &entity[num]);
            num += 1;
        }
    }
    assert(num == num_blocks);

    grid->cells.num_inner = cell_off[grid->entities.num_inner];
    grid->cells.off_boundary = cell_off[grid->entities.off_boundary];
    grid->cells.off_periodic = cell_off[grid->entities.num];
    assert(grid->cells.off_periodic == grid->cells.num);

    grid->entities.name = name;
    grid->entities.cell_off = cell_off;
    grid->entities.periodic = periodic;
    grid->entities.rotation = rotation;
    grid->entities.translation = translation;
    grid->entities.node_off = node_off;
    grid->entities.node_idx = grid_tag_to_idx(grid, node_tag, num_nodes);

    teal2_free(node_tag);
}

Grid *grid_gmsh(const char *fname)
{
    arena = arena2_init(0);

    Gmsh *gmsh = read(fname);
    Entity *entity = extract(gmsh);

    Grid *grid = teal2_calloc(1, sizeof(*grid));
    create_nodes(grid, gmsh);
    create_cells(grid, entity, gmsh);
    create_entities(grid, entity, gmsh);

    arena2_deinit(arena);
    return grid;
}
