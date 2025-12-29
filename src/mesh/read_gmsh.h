#include <assert.h>
#include <stdint.h>
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
} Name;

typedef struct {
    int32_t num;
    Name *name;
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
    uint64_t tot_nodes;
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
} Link;

typedef struct {
    uint64_t num;
    Link *link;
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
    Name *name = teal_alloc(gmsh->physicals.num, sizeof(*name));
    for (int i = 0; i < gmsh->physicals.num; i++) {
        parse_ascii(file, &name[i].dim, 1, MPI_INT32_T);
        parse_ascii(file, &name[i].tag, 1, MPI_INT32_T);
        parse_string(file, name[i].name, sizeof(name->name));
    }
    gmsh->physicals.name = name;
}

static Point *read_points(long num, Parse *file, int mode)
{
    Point *point = teal_alloc(num, sizeof(*point));
    for (long i = 0; i < num; i++) {
        parse(file, &point[i].tag, 1, MPI_INT32_T, mode);
        parse(file, &point[i].x, 1, MPI_DOUBLE, mode);
        parse(file, &point[i].y, 1, MPI_DOUBLE, mode);
        parse(file, &point[i].z, 1, MPI_DOUBLE, mode);
        parse(file, &point[i].num_physical_tags, 1, MPI_UINT64_T, mode);
        uint64_t num_physical_tags = point[i].num_physical_tags;
        int32_t *physical_tag = teal_alloc(num_physical_tags, sizeof(*physical_tag));
        parse(file, physical_tag, num_physical_tags, MPI_INT32_T, mode);
        point[i].physical_tag = physical_tag;
    }
    return point;
}

static Curve *read_curves(long num, Parse *file, int mode)
{
    Curve *curve = teal_alloc(num, sizeof(*curve));
    for (long i = 0; i < num; i++) {
        parse(file, &curve[i].tag, 1, MPI_INT32_T, mode);
        parse(file, &curve[i].min_x, 1, MPI_DOUBLE, mode);
        parse(file, &curve[i].min_y, 1, MPI_DOUBLE, mode);
        parse(file, &curve[i].min_z, 1, MPI_DOUBLE, mode);
        parse(file, &curve[i].max_x, 1, MPI_DOUBLE, mode);
        parse(file, &curve[i].max_y, 1, MPI_DOUBLE, mode);
        parse(file, &curve[i].max_z, 1, MPI_DOUBLE, mode);
        parse(file, &curve[i].num_physical_tags, 1, MPI_UINT64_T, mode);
        uint64_t num_physical_tags = curve[i].num_physical_tags;
        int32_t *physical_tag = teal_alloc(num_physical_tags, sizeof(*physical_tag));
        parse(file, physical_tag, num_physical_tags, MPI_INT32_T, mode);
        curve[i].physical_tag = physical_tag;
        parse(file, &curve[i].num_point_tags, 1, MPI_UINT64_T, mode);
        uint64_t num_point_tags = curve[i].num_point_tags;
        int32_t *point_tag = teal_alloc(num_point_tags, sizeof(*point_tag));
        parse(file, point_tag, num_point_tags, MPI_INT32_T, mode);
        curve[i].point_tag = point_tag;
    }
    return curve;
}

static Surface *read_surfaces(long num, Parse *file, int mode)
{
    Surface *surface = teal_alloc(num, sizeof(*surface));
    for (long i = 0; i < num; i++) {
        parse(file, &surface[i].tag, 1, MPI_INT32_T, mode);
        parse(file, &surface[i].min_x, 1, MPI_DOUBLE, mode);
        parse(file, &surface[i].min_y, 1, MPI_DOUBLE, mode);
        parse(file, &surface[i].min_z, 1, MPI_DOUBLE, mode);
        parse(file, &surface[i].max_x, 1, MPI_DOUBLE, mode);
        parse(file, &surface[i].max_y, 1, MPI_DOUBLE, mode);
        parse(file, &surface[i].max_z, 1, MPI_DOUBLE, mode);
        parse(file, &surface[i].num_physical_tags, 1, MPI_UINT64_T, mode);
        uint64_t num_physical_tags = surface[i].num_physical_tags;
        int32_t *physical_tag = teal_alloc(num_physical_tags, sizeof(*physical_tag));
        parse(file, physical_tag, num_physical_tags, MPI_INT32_T, mode);
        surface[i].physical_tag = physical_tag;
        parse(file, &surface[i].num_curve_tags, 1, MPI_UINT64_T, mode);
        uint64_t num_curve_tags = surface[i].num_curve_tags;
        int32_t *curve_tag = teal_alloc(num_curve_tags, sizeof(*curve_tag));
        parse(file, curve_tag, num_curve_tags, MPI_INT32_T, mode);
        surface[i].curve_tag = curve_tag;
    }
    return surface;
}

static Volume *read_volumes(long num, Parse *file, int mode)
{
    Volume *volume = teal_alloc(num, sizeof(*volume));
    for (long i = 0; i < num; i++) {
        parse(file, &volume[i].tag, 1, MPI_INT32_T, mode);
        parse(file, &volume[i].min_x, 1, MPI_DOUBLE, mode);
        parse(file, &volume[i].min_y, 1, MPI_DOUBLE, mode);
        parse(file, &volume[i].min_z, 1, MPI_DOUBLE, mode);
        parse(file, &volume[i].max_x, 1, MPI_DOUBLE, mode);
        parse(file, &volume[i].max_y, 1, MPI_DOUBLE, mode);
        parse(file, &volume[i].max_z, 1, MPI_DOUBLE, mode);
        parse(file, &volume[i].num_physical_tags, 1, MPI_UINT64_T, mode);
        uint64_t num_physical_tags = volume[i].num_physical_tags;
        int32_t *physical_tag = teal_alloc(num_physical_tags, sizeof(*physical_tag));
        parse(file, physical_tag, num_physical_tags, MPI_INT32_T, mode);
        volume[i].physical_tag = physical_tag;
        parse(file, &volume[i].num_surface_tags, 1, MPI_UINT64_T, mode);
        uint64_t num_surface_tags = volume[i].num_surface_tags;
        int32_t *surface_tag = teal_alloc(num_surface_tags, sizeof(*surface_tag));
        parse(file, surface_tag, num_surface_tags, MPI_INT32_T, mode);
        volume[i].surface_tag = surface_tag;
    }
    return volume;
}

static void read_entities(Gmsh *gmsh, Parse *file, int mode)
{
    parse(file, &gmsh->entities.num_points, 1, MPI_UINT64_T, mode);
    parse(file, &gmsh->entities.num_curves, 1, MPI_UINT64_T, mode);
    parse(file, &gmsh->entities.num_surfaces, 1, MPI_UINT64_T, mode);
    parse(file, &gmsh->entities.num_volumes, 1, MPI_UINT64_T, mode);
    gmsh->entities.point = read_points(gmsh->entities.num_points, file, mode);
    gmsh->entities.curve = read_curves(gmsh->entities.num_curves, file, mode);
    gmsh->entities.surface = read_surfaces(gmsh->entities.num_surfaces, file, mode);
    gmsh->entities.volume = read_volumes(gmsh->entities.num_volumes, file, mode);
}

static long read_node_block(NodeBlock *block, long beg, long end, long off, Parse *file, int mode)
{
    parse(file, &block->entity_dim, 1, MPI_INT32_T, mode);
    parse(file, &block->entity_tag, 1, MPI_INT32_T, mode);
    parse(file, &block->parametric, 1, MPI_INT32_T, mode);
    parse(file, &block->num_nodes, 1, MPI_UINT64_T, mode);

    long tot_nodes = block->num_nodes;
    long num_nodes = max(0, min(end, off + tot_nodes) - max(beg, off));
    assert(tot_nodes == sync_lsum(num_nodes));

    block->num_nodes = num_nodes;
    block->tag = teal_alloc(num_nodes, sizeof(*block->tag));
    parse(file, block->tag, num_nodes, MPI_UINT64_T, mode | SPLIT);

    long len = 3 + (block->parametric ? block->entity_dim : 0);
    block->coord = teal_alloc(num_nodes * len, sizeof(*block->coord));
    parse(file, block->coord, num_nodes * len, MPI_DOUBLE, mode | SPLIT);

    return tot_nodes;
}

static void read_nodes(Gmsh *gmsh, Parse *file, int mode)
{
    parse(file, &gmsh->nodes.num_blocks, 1, MPI_UINT64_T, mode);
    parse(file, &gmsh->nodes.tot_nodes, 1, MPI_UINT64_T, mode);
    parse(file, &gmsh->nodes.min_tag, 1, MPI_UINT64_T, mode);
    parse(file, &gmsh->nodes.max_tag, 1, MPI_UINT64_T, mode);

    long num_blocks = gmsh->nodes.num_blocks;
    gmsh->nodes.block = teal_alloc(num_blocks, sizeof(*gmsh->nodes.block));

    long tot_nodes = gmsh->nodes.tot_nodes;
    long base = tot_nodes / sync.size;
    long extra = tot_nodes % sync.size;
    long beg = (sync.rank * base) + ((sync.rank < extra) ? sync.rank : extra);
    long end = beg + base + (sync.rank < extra);

    long off = 0;
    for (long i = 0; i < num_blocks; i++) {
        off += read_node_block(&gmsh->nodes.block[i], beg, end, off, file, mode);
    }
    assert(off == tot_nodes);
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
    parse(file, &block->entity_tag, 1, MPI_INT32_T, mode);
    parse(file, &block->element_type, 1, MPI_INT32_T, mode);
    parse(file, &block->num_elements, 1, MPI_UINT64_T, mode);

    long tot_elements = block->num_elements;
    long num_elements = max(0, min(end, off + tot_elements) - max(beg, off));
    assert(tot_elements == sync_lsum(num_elements));

    block->num_elements = num_elements;
    block->tag = teal_alloc(num_elements, sizeof(*block->tag));

    long len = num_node_tags(block->element_type);
    block->node_tag = teal_alloc(num_elements * len, sizeof(*block->node_tag));

    long stride = 1 + len;
    uint64_t (*tmp)[stride] = teal_alloc(num_elements, sizeof(*tmp));
    parse(file, tmp, num_elements * stride, MPI_UINT64_T, mode | SPLIT);
    for (long i = 0; i < num_elements; i++) {
        block->tag[i] = tmp[i][0];
        for (long j = 0; j < len; j++) {
            block->node_tag[(i * len) + j] = tmp[i][1 + j];
        }
    }
    teal_free(tmp);

    return tot_elements;
}

static void read_elements(Gmsh *gmsh, Parse *file, int mode)
{
    parse(file, &gmsh->elements.num_blocks, 1, MPI_UINT64_T, mode);
    parse(file, &gmsh->elements.tot_elements, 1, MPI_UINT64_T, mode);
    parse(file, &gmsh->elements.min_tag, 1, MPI_UINT64_T, mode);
    parse(file, &gmsh->elements.max_tag, 1, MPI_UINT64_T, mode);

    long num_blocks = gmsh->elements.num_blocks;
    gmsh->elements.block = teal_alloc(num_blocks, sizeof(*gmsh->elements.block));

    long tot_elements = gmsh->elements.tot_elements;
    long base = tot_elements / sync.size;
    long extra = tot_elements % sync.size;
    long beg = (sync.rank * base) + ((sync.rank < extra) ? sync.rank : extra);
    long end = beg + base + (sync.rank < extra);

    long off = 0;
    for (long i = 0; i < num_blocks; i++) {
        off += read_element_block(&gmsh->elements.block[i], beg, end, off, file, mode);
    }
    assert(off == tot_elements);
}

static void read_link(Link *link, Parse *file, int mode)
{
    parse(file, &link->entity_dim, 1, MPI_INT32_T, mode);
    parse(file, &link->entity_tag, 1, MPI_INT32_T, mode);
    parse(file, &link->entity_tag_master, 1, MPI_INT32_T, mode);
    parse(file, &link->num_affine, 1, MPI_UINT64_T, mode);

    long num_affine = link->num_affine;
    link->affine = teal_alloc(num_affine, sizeof(*link->affine));
    parse(file, link->affine, num_affine, MPI_DOUBLE, mode);

    parse(file, &link->num_nodes, 1, MPI_UINT64_T, mode);
    long num_nodes = link->num_nodes;
    link->node_tag = teal_alloc(num_nodes, sizeof(*link->node_tag));
    link->node_tag_master = teal_alloc(num_nodes, sizeof(*link->node_tag_master));

    uint64_t (*tmp)[2] = teal_alloc(num_nodes, sizeof(*tmp));
    parse(file, tmp, num_nodes * 2, MPI_UINT64_T, mode);
    for (long i = 0; i < num_nodes; i++) {
        link->node_tag[i] = tmp[i][0];
        link->node_tag_master[i] = tmp[i][1];
    }
    teal_free(tmp);
}

static void read_periodics(Gmsh *gmsh, Parse *file, int mode)
{
    parse(file, &gmsh->periodics.num, 1, MPI_UINT64_T, mode);

    long num = gmsh->periodics.num;
    gmsh->periodics.link = teal_alloc(num, sizeof(*gmsh->periodics.link));
    for (long i = 0; i < num; i++) {
        read_link(&gmsh->periodics.link[i], file, mode);
    }
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

static void gmsh_deinit(Gmsh *gmsh)
{
    if (!gmsh) {
        return;
    }

    teal_free(gmsh->physicals.name);

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

    if (gmsh->periodics.link) {
        for (uint64_t i = 0; i < gmsh->periodics.num; i++) {
            teal_free(gmsh->periodics.link[i].affine);
            teal_free(gmsh->periodics.link[i].node_tag);
            teal_free(gmsh->periodics.link[i].node_tag_master);
        }
        teal_free(gmsh->periodics.link);
    }

    teal_free(gmsh);
}

static void read_gmsh(Mesh *mesh, const char *fname)
{
    Gmsh *gmsh = gmsh_init(fname);
    gmsh_deinit(gmsh);
}
