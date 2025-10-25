#include <inttypes.h>
#include <stdlib.h>
#include <string.h>

#include "mesh.h"
#include "reorder.h"
#include "teal/arena.h"
#include "teal/assert.h"
#include "teal/parse.h"
#include "teal/sync.h"
#include "teal/utils.h"

static void read_format(double *version, Mode *mode, Type *SIZE, bool *swap, File file)
{
    char line[128];
    parse_ascii(STR, line, sizeof(line), file);
    assert(!strcmp(line, "$MeshFormat"));

    parse_ascii(F64, version, 1, file);

    int32_t file_type;
    parse_ascii(I32, &file_type, 1, file);
    *mode = (file_type == 1) ? BINARY : ASCII;

    int32_t data_size;
    parse_ascii(I32, &data_size, 1, file);
    *SIZE = (data_size == sizeof(uint32_t)) ? U32 : U64;

    if (*mode == BINARY) {
        int32_t one;
        parse_binary(I32, &one, 1, false, file);
        *swap = (one != 1);
    }

    parse_ascii(STR, line, sizeof(line), file);
    assert(!strcmp(line, "$EndMeshFormat"));
}

typedef struct {
    int32_t dim;
    int32_t tag;
    char name[128];
} Physical;

typedef struct {
    int32_t num;
    Physical *physical;
} Physicals;

static void read_physicals(Physicals *physicals, File file)
{
    char line[128];
    parse_ascii(STR, line, sizeof(line), file);
    assert(!strcmp(line, "$PhysicalNames"));

    parse_ascii(I32, &physicals->num, 1, file);

    physicals->physical = arena_malloc(physicals->num, sizeof(*physicals->physical));
    for (long i = 0; i < physicals->num; i++) {
        parse_ascii(I32, &physicals->physical[i].dim, 1, file);
        parse_ascii(I32, &physicals->physical[i].tag, 1, file);
        parse_ascii(STR, physicals->physical[i].name, sizeof(physicals->physical[i].name), file);
    }

    parse_ascii(STR, line, sizeof(line), file);
    assert(!strcmp(line, "$EndPhysicalNames"));
}

typedef union {
    uint32_t u32;
    uint64_t u64;
} Size;

typedef struct {
    int32_t tag;
    double coords[3];
    Size num_physical_tags;
    int32_t *physical_tag;
} Point;

static Point *read_points(Mode mode, Type SIZE, long num, bool swap, File file)
{
    Point *point = arena_malloc(num, sizeof(*point));
    for (long i = 0; i < num; i++) {
        parse(mode, I32, &point[i].tag, 1, swap, file);
        parse(mode, F64, point[i].coords, 3, swap, file);
        parse(mode, SIZE, &point[i].num_physical_tags, 1, swap, file);
        long num_physical_tags = parse_data_to_long(SIZE, &point[i].num_physical_tags, 0);
        point[i].physical_tag = arena_malloc(num_physical_tags, sizeof(*point[i].physical_tag));
        parse(mode, I32, point[i].physical_tag, num_physical_tags, swap, file);
    }
    return point;
}

typedef struct {
    int32_t tag;
    double bounds[6];
    Size num_physical_tags;
    int32_t *physical_tag;
    Size num_bounding_points;
    int32_t *point_tag;
} Curve;

static Curve *read_curves(Mode mode, Type SIZE, long num, bool swap, File file)
{
    Curve *curve = arena_malloc(num, sizeof(*curve));
    for (long i = 0; i < num; i++) {
        parse(mode, I32, &curve[i].tag, 1, swap, file);
        parse(mode, F64, curve[i].bounds, 6, swap, file);
        parse(mode, SIZE, &curve[i].num_physical_tags, 1, swap, file);
        long num_physical_tags = parse_data_to_long(SIZE, &curve[i].num_physical_tags, 0);
        curve[i].physical_tag = arena_malloc(num_physical_tags, sizeof(*curve[i].physical_tag));
        parse(mode, I32, curve[i].physical_tag, num_physical_tags, swap, file);
        parse(mode, SIZE, &curve[i].num_bounding_points, 1, swap, file);
        long num_bounding_points = parse_data_to_long(SIZE, &curve[i].num_bounding_points, 0);
        curve[i].point_tag = arena_malloc(num_bounding_points, sizeof(*curve[i].point_tag));
        parse(mode, I32, curve[i].point_tag, num_bounding_points, swap, file);
    }
    return curve;
}

typedef struct {
    int32_t tag;
    double bounds[6];
    Size num_physical_tags;
    int32_t *physical_tag;
    Size num_bounding_curves;
    int32_t *curve_tag;
} Surface;

static Surface *read_surfaces(Mode mode, Type SIZE, long num, bool swap, File file)
{
    Surface *surface = arena_malloc(num, sizeof(*surface));
    for (long i = 0; i < num; i++) {
        parse(mode, I32, &surface[i].tag, 1, swap, file);
        parse(mode, F64, surface[i].bounds, 6, swap, file);
        parse(mode, SIZE, &surface[i].num_physical_tags, 1, swap, file);
        long num_physical_tags = parse_data_to_long(SIZE, &surface[i].num_physical_tags, 0);
        surface[i].physical_tag = arena_malloc(num_physical_tags, sizeof(*surface[i].physical_tag));
        parse(mode, I32, surface[i].physical_tag, num_physical_tags, swap, file);
        parse(mode, SIZE, &surface[i].num_bounding_curves, 1, swap, file);
        long num_bounding_curves = parse_data_to_long(SIZE, &surface[i].num_bounding_curves, 0);
        surface[i].curve_tag = arena_malloc(num_bounding_curves, sizeof(*surface[i].curve_tag));
        parse(mode, I32, surface[i].curve_tag, num_bounding_curves, swap, file);
    }
    return surface;
}

typedef struct {
    int32_t tag;
    double bounds[6];
    Size num_physical_tags;
    int32_t *physical_tag;
    Size num_bounding_surfaces;
    int32_t *surface_tag;
} Volume;

static Volume *read_volumes(Mode mode, Type SIZE, long num, bool swap, File file)
{
    Volume *volume = arena_malloc(num, sizeof(*volume));
    for (long i = 0; i < num; i++) {
        parse(mode, I32, &volume[i].tag, 1, swap, file);
        parse(mode, F64, volume[i].bounds, 6, swap, file);
        parse(mode, SIZE, &volume[i].num_physical_tags, 1, swap, file);
        long num_physical_tags = parse_data_to_long(SIZE, &volume[i].num_physical_tags, 0);
        volume[i].physical_tag = arena_malloc(num_physical_tags, sizeof(*volume[i].physical_tag));
        parse(mode, I32, volume[i].physical_tag, num_physical_tags, swap, file);
        parse(mode, SIZE, &volume[i].num_bounding_surfaces, 1, swap, file);
        long num_bounding_surfaces = parse_data_to_long(SIZE, &volume[i].num_bounding_surfaces, 0);
        volume[i].surface_tag = arena_malloc(num_bounding_surfaces, sizeof(*volume[i].surface_tag));
        parse(mode, I32, volume[i].surface_tag, num_bounding_surfaces, swap, file);
    }
    return volume;
}

typedef struct {
    Size num_points;
    Size num_curves;
    Size num_surfaces;
    Size num_volumes;
    Point *point;
    Curve *curve;
    Surface *surface;
    Volume *volume;
} Entities;

static void read_entities(Entities *entities, Mode mode, Type SIZE, bool swap, File file)
{
    char line[128];
    parse_ascii(STR, line, sizeof(line), file);
    assert(!strcmp(line, "$Entities"));

    parse(mode, SIZE, &entities->num_points, 1, swap, file);
    parse(mode, SIZE, &entities->num_curves, 1, swap, file);
    parse(mode, SIZE, &entities->num_surfaces, 1, swap, file);
    parse(mode, SIZE, &entities->num_volumes, 1, swap, file);

    long num_points = parse_data_to_long(SIZE, &entities->num_points, 0);
    entities->point = read_points(mode, SIZE, num_points, swap, file);

    long num_curves = parse_data_to_long(SIZE, &entities->num_curves, 0);
    entities->curve = read_curves(mode, SIZE, num_curves, swap, file);

    long num_surfaces = parse_data_to_long(SIZE, &entities->num_surfaces, 0);
    entities->surface = read_surfaces(mode, SIZE, num_surfaces, swap, file);

    long num_volumes = parse_data_to_long(SIZE, &entities->num_volumes, 0);
    entities->volume = read_volumes(mode, SIZE, num_volumes, swap, file);

    parse_ascii(STR, line, sizeof(line), file);
    assert(!strcmp(line, "$EndEntities"));
}

typedef union {
    uint32_t *u32;
    uint64_t *u64;
} SizePtr;

typedef struct {
    int32_t entity_dim;
    int32_t entity_tag;
    int32_t parametric;
    Size num_nodes;
    SizePtr tag;
    double *coord;
} NodeBlock;

static long read_node_block(NodeBlock *block, long beg, long end, long off, Mode mode, Type SIZE,
                            bool swap, File file)
{
    parse(mode, I32, &block->entity_dim, 1, swap, file);
    parse(mode, I32, &block->entity_tag, 1, swap, file);
    parse(mode, I32, &block->parametric, 1, swap, file);
    parse(mode, SIZE, &block->num_nodes, 1, swap, file);

    long tot_nodes = parse_data_to_long(SIZE, &block->num_nodes, 0);
    long beg_nodes = lmax(beg, off);
    long end_nodes = lmin(end, off + tot_nodes);
    long num_nodes = lmax(0, end_nodes - beg_nodes);
    assert(tot_nodes == sync_lsum(num_nodes));

    switch (SIZE) {
        case U32: {
            block->num_nodes.u32 = num_nodes;
            block->tag.u32 = arena_malloc(num_nodes, sizeof(*block->tag.u32));
            long stride = (mode == ASCII) ? 1 : sizeof(*block->tag.u32);
            parse_split(mode, U32, block->tag.u32, num_nodes, 1, stride, swap, file);
            break;
        }
        case U64: {
            block->num_nodes.u64 = num_nodes;
            block->tag.u64 = arena_malloc(num_nodes, sizeof(*block->tag.u64));
            long stride = (mode == ASCII) ? 1 : sizeof(*block->tag.u64);
            parse_split(mode, U64, block->tag.u64, num_nodes, 1, stride, swap, file);
            break;
        }
        default: assert(false);
    }

    long len = 3 + (block->parametric ? block->entity_dim : 0);
    block->coord = arena_malloc(num_nodes * len, sizeof(*block->coord));
    long stride = (mode == ASCII) ? len : len * sizeof(*block->coord);
    parse_split(mode, F64, block->coord, num_nodes, len, stride, swap, file);

    return off + tot_nodes;
}

typedef struct {
    Size num_blocks;
    Size tot_nodes;
    Size min_tag;
    Size max_tag;
    NodeBlock *block;
} Nodes;

static void read_nodes(Nodes *nodes, Mode mode, Type SIZE, bool swap, File file)
{
    char line[128];
    parse_ascii(STR, line, sizeof(line), file);
    assert(!strcmp(line, "$Nodes"));

    parse(mode, SIZE, &nodes->num_blocks, 1, swap, file);
    parse(mode, SIZE, &nodes->tot_nodes, 1, swap, file);
    parse(mode, SIZE, &nodes->min_tag, 1, swap, file);
    parse(mode, SIZE, &nodes->max_tag, 1, swap, file);

    long num_blocks = parse_data_to_long(SIZE, &nodes->num_blocks, 0);
    nodes->block = arena_malloc(num_blocks, sizeof(*nodes->block));

    long tot_nodes = parse_data_to_long(SIZE, &nodes->tot_nodes, 0);
    long base = tot_nodes / sync.size;
    long extra = tot_nodes % sync.size;
    long beg = (sync.rank * base) + ((sync.rank < extra) ? sync.rank : extra);
    long end = beg + base + (sync.rank < extra);

    long off = 0;
    for (long i = 0; i < num_blocks; i++) {
        off = read_node_block(&nodes->block[i], beg, end, off, mode, SIZE, swap, file);
    }
    assert(off == tot_nodes);

    parse_ascii(STR, line, sizeof(line), file);
    assert(!strcmp(line, "$EndNodes"));
}

typedef struct {
    int32_t entity_dim;
    int32_t entity_tag;
    int32_t element_type;
    Size num_elements;
    SizePtr tag;
    SizePtr node_tag;
} ElementBlock;

static long num_node_tags(long element_type)
{
    switch (element_type) {
        case 2: return 3;  // triangle
        case 3:            // quadrangle
        case 4: return 4;  // tetrahedron
        case 5: return 8;  // hexahedron
        case 6: return 6;  // prism
        case 7: return 5;  // pyramid
        default: assert(false);
    }
}

static long read_element_block(ElementBlock *block, long beg, long end, long off, Mode mode,
                               Type SIZE, bool swap, File file)
{
    parse(mode, I32, &block->entity_dim, 1, swap, file);
    parse(mode, I32, &block->entity_tag, 1, swap, file);
    parse(mode, I32, &block->element_type, 1, swap, file);
    parse(mode, SIZE, &block->num_elements, 1, swap, file);

    long tot_elements = parse_data_to_long(SIZE, &block->num_elements, 0);
    long beg_elements = lmax(beg, off);
    long end_elements = lmin(end, off + tot_elements);
    long num_elements = lmax(0, end_elements - beg_elements);
    assert(tot_elements == sync_lsum(num_elements));

    long len = num_node_tags(block->element_type);
    long offset;
    switch (SIZE) {
        case U32: {
            block->num_elements.u32 = num_elements;
            block->tag.u32 = arena_malloc(num_elements, sizeof(*block->tag.u32));
            block->node_tag.u32 = arena_malloc(num_elements * len, sizeof(*block->node_tag.u32));
            long stride = (mode == ASCII)
                              ? 1 + len
                              : sizeof(*block->tag.u32) + (len * sizeof(*block->node_tag.u32));
            parse_split(mode, U32, block->tag.u32, num_elements, 1, stride, swap, file);
            offset =
                parse_split(mode, U32, block->node_tag.u32, num_elements, len, stride, swap, file);
            break;
        }
        case U64: {
            block->num_elements.u64 = num_elements;
            block->tag.u64 = arena_malloc(num_elements, sizeof(*block->tag.u64));
            block->node_tag.u64 = arena_malloc(num_elements * len, sizeof(*block->node_tag.u64));
            long stride = (mode == ASCII)
                              ? 1 + len
                              : sizeof(*block->tag.u64) + (len * sizeof(*block->node_tag.u64));
            parse_split(mode, U64, block->tag.u64, num_elements, 1, stride, swap, file);
            offset =
                parse_split(mode, U64, block->node_tag.u64, num_elements, len, stride, swap, file);
            break;
        }
        default: assert(false);
    }
    parse_set_offset(file, offset);

    return off + tot_elements;
}

typedef struct {
    Size num_blocks;
    Size tot_elements;
    Size min_tag;
    Size max_tag;
    ElementBlock *block;
} Elements;

static void read_elements(Elements *elements, Mode mode, Type SIZE, bool swap, File file)
{
    char line[128];
    parse_ascii(STR, line, sizeof(line), file);
    assert(!strcmp(line, "$Elements"));

    parse(mode, SIZE, &elements->num_blocks, 1, swap, file);
    parse(mode, SIZE, &elements->tot_elements, 1, swap, file);
    parse(mode, SIZE, &elements->min_tag, 1, swap, file);
    parse(mode, SIZE, &elements->max_tag, 1, swap, file);

    long num_blocks = parse_data_to_long(SIZE, &elements->num_blocks, 0);
    elements->block = arena_malloc(num_blocks, sizeof(*elements->block));

    long tot_elements = parse_data_to_long(SIZE, &elements->tot_elements, 0);
    long base = tot_elements / sync.size;
    long extra = tot_elements % sync.size;
    long beg = (sync.rank * base) + ((sync.rank < extra) ? sync.rank : extra);
    long end = beg + base + (sync.rank < extra);

    long off = 0;
    for (long i = 0; i < num_blocks; i++) {
        off = read_element_block(&elements->block[i], beg, end, off, mode, SIZE, swap, file);
    }
    assert(off == tot_elements);

    parse_ascii(STR, line, sizeof(line), file);
    assert(!strcmp(line, "$EndElements"));
}

typedef struct {
    Physicals physicals;
    Entities entities;
    Nodes nodes;
    Elements elements;
} Gmsh;

__attribute((unused)) static void gmsh_dump(FILE *stream, double version, Mode mode, Type SIZE,
                                            bool swap, const Gmsh *gmsh)
{
    assert(stream && gmsh);
    long off = 0;
    for (long rank = 0; rank < sync.size; rank++) {
        if (rank == sync.rank) {
            fseek(stream, off, SEEK_SET);
            fprintf(stream, "rank %d:\n", sync.rank);

            fprintf(stream, "\t $MeshFormat\n");
            fprintf(stream, "\t %g\n", version);
            fprintf(stream, "\t %d\n", (mode == ASCII) ? 0 : 1);
            fprintf(stream, "\t %d\n", (SIZE == U32) ? 4 : 8);
            fprintf(stream, "\t %s\n", swap ? "true" : "false");
            fprintf(stream, "\t $EndMeshFormat\n");

            Physicals physicals = gmsh->physicals;
            fprintf(stream, "\t $PhysicalNames\n");
            fprintf(stream, "\t %d\n", physicals.num);
            for (long i = 0; i < physicals.num; i++) {
                fprintf(stream, "\t %d %d %s\n", physicals.physical[i].dim,
                        physicals.physical[i].tag, physicals.physical[i].name);
            }
            fprintf(stream, "\t $EndPhysicalNames\n");

            Entities entities = gmsh->entities;
            fprintf(stream, "\t $Entities\n");
            switch (SIZE) {
                case U32:
                    fprintf(stream, "\t %" PRIu32 " %" PRIu32 " %" PRIu32 " %" PRIu32 "\n",
                            entities.num_points.u32, entities.num_curves.u32,
                            entities.num_surfaces.u32, entities.num_volumes.u32);
                    break;
                case U64:
                    fprintf(stream, "\t %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 "\n",
                            entities.num_points.u64, entities.num_curves.u64,
                            entities.num_surfaces.u64, entities.num_volumes.u64);
                    break;
                default: assert(false);
            }
            for (long i = 0; i < parse_data_to_long(SIZE, &entities.num_points, 0); i++) {
                Point point = entities.point[i];
                switch (SIZE) {
                    case U32:
                        fprintf(stream, "\t %d %g %g %g %" PRIu32, point.tag, point.coords[0],
                                point.coords[1], point.coords[2], point.num_physical_tags.u32);
                        for (long j = 0; j < (long)point.num_physical_tags.u32; j++) {
                            fprintf(stream, " %d", point.physical_tag[j]);
                        }
                        break;
                    case U64:
                        fprintf(stream, "\t %d %g %g %g %" PRIu64, point.tag, point.coords[0],
                                point.coords[1], point.coords[2], point.num_physical_tags.u64);
                        for (long j = 0; j < (long)point.num_physical_tags.u64; j++) {
                            fprintf(stream, " %d", point.physical_tag[j]);
                        }
                        break;
                    default: assert(false);
                }
                fprintf(stream, "\n");
            }
            for (long i = 0; i < parse_data_to_long(SIZE, &entities.num_curves, 0); i++) {
                Curve curve = entities.curve[i];
                switch (SIZE) {
                    case U32:
                        fprintf(stream, "\t %d %g %g %g %g %g %g %" PRIu32, curve.tag,
                                curve.bounds[0], curve.bounds[1], curve.bounds[2], curve.bounds[3],
                                curve.bounds[4], curve.bounds[5], curve.num_physical_tags.u32);
                        for (long j = 0; j < (long)curve.num_physical_tags.u32; j++) {
                            fprintf(stream, " %d", curve.physical_tag[j]);
                        }
                        fprintf(stream, " %" PRIu32, curve.num_bounding_points.u32);
                        for (long j = 0; j < (long)curve.num_bounding_points.u32; j++) {
                            fprintf(stream, " %d", curve.point_tag[j]);
                        }
                        break;
                    case U64:
                        fprintf(stream, "\t %d %g %g %g %g %g %g %" PRIu64, curve.tag,
                                curve.bounds[0], curve.bounds[1], curve.bounds[2], curve.bounds[3],
                                curve.bounds[4], curve.bounds[5], curve.num_physical_tags.u64);
                        for (long j = 0; j < (long)curve.num_physical_tags.u64; j++) {
                            fprintf(stream, " %d", curve.physical_tag[j]);
                        }
                        fprintf(stream, " %" PRIu64, curve.num_bounding_points.u64);
                        for (long j = 0; j < (long)curve.num_bounding_points.u64; j++) {
                            fprintf(stream, " %d", curve.point_tag[j]);
                        }
                        break;
                    default: assert(false);
                }
                fprintf(stream, "\n");
            }
            for (long i = 0; i < parse_data_to_long(SIZE, &entities.num_surfaces, 0); i++) {
                Surface surface = entities.surface[i];
                switch (SIZE) {
                    case U32:
                        fprintf(stream, "\t %d %g %g %g %g %g %g %" PRIu32, surface.tag,
                                surface.bounds[0], surface.bounds[1], surface.bounds[2],
                                surface.bounds[3], surface.bounds[4], surface.bounds[5],
                                surface.num_physical_tags.u32);
                        for (long j = 0; j < (long)surface.num_physical_tags.u32; j++) {
                            fprintf(stream, " %d", surface.physical_tag[j]);
                        }
                        fprintf(stream, " %" PRIu32, surface.num_bounding_curves.u32);
                        for (long j = 0; j < (long)surface.num_bounding_curves.u32; j++) {
                            fprintf(stream, " %d", surface.curve_tag[j]);
                        }
                        break;
                    case U64:
                        fprintf(stream, "\t %d %g %g %g %g %g %g %" PRIu64, surface.tag,
                                surface.bounds[0], surface.bounds[1], surface.bounds[2],
                                surface.bounds[3], surface.bounds[4], surface.bounds[5],
                                surface.num_physical_tags.u64);
                        for (long j = 0; j < (long)surface.num_physical_tags.u64; j++) {
                            fprintf(stream, " %d", surface.physical_tag[j]);
                        }
                        fprintf(stream, " %" PRIu64, surface.num_bounding_curves.u64);
                        for (long j = 0; j < (long)surface.num_bounding_curves.u64; j++) {
                            fprintf(stream, " %d", surface.curve_tag[j]);
                        }
                        break;
                    default: assert(false);
                }
                fprintf(stream, "\n");
            }
            for (long i = 0; i < parse_data_to_long(SIZE, &entities.num_volumes, 0); i++) {
                Volume volume = entities.volume[i];
                switch (SIZE) {
                    case U32:
                        fprintf(stream, "\t %d %g %g %g %g %g %g %" PRIu32, volume.tag,
                                volume.bounds[0], volume.bounds[1], volume.bounds[2],
                                volume.bounds[3], volume.bounds[4], volume.bounds[5],
                                volume.num_physical_tags.u32);
                        for (long j = 0; j < (long)volume.num_physical_tags.u32; j++) {
                            fprintf(stream, " %d", volume.physical_tag[j]);
                        }
                        fprintf(stream, " %" PRIu32, volume.num_bounding_surfaces.u32);
                        for (long j = 0; j < (long)volume.num_bounding_surfaces.u32; j++) {
                            fprintf(stream, " %d", volume.surface_tag[j]);
                        }
                        break;
                    case U64:
                        fprintf(stream, "\t %d %g %g %g %g %g %g %" PRIu64, volume.tag,
                                volume.bounds[0], volume.bounds[1], volume.bounds[2],
                                volume.bounds[3], volume.bounds[4], volume.bounds[5],
                                volume.num_physical_tags.u64);
                        for (long j = 0; j < (long)volume.num_physical_tags.u64; j++) {
                            fprintf(stream, " %d", volume.physical_tag[j]);
                        }
                        fprintf(stream, " %" PRIu64, volume.num_bounding_surfaces.u64);
                        for (long j = 0; j < (long)volume.num_bounding_surfaces.u64; j++) {
                            fprintf(stream, " %d", volume.surface_tag[j]);
                        }
                        break;
                    default: assert(false);
                }
                fprintf(stream, "\n");
            }
            fprintf(stream, "\t $EndEntities\n");

            Nodes nodes = gmsh->nodes;
            fprintf(stream, "\t $Nodes\n");
            switch (SIZE) {
                case U32:
                    fprintf(stream, "\t %" PRIu32 " %" PRIu32 " %" PRIu32 " %" PRIu32 "\n",
                            nodes.num_blocks.u32, nodes.tot_nodes.u32, nodes.min_tag.u32,
                            nodes.max_tag.u32);
                    break;
                case U64:
                    fprintf(stream, "\t %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 "\n",
                            nodes.num_blocks.u64, nodes.tot_nodes.u64, nodes.min_tag.u64,
                            nodes.max_tag.u64);
                    break;
                default: assert(false);
            }
            for (long i = 0; i < parse_data_to_long(SIZE, &nodes.num_blocks, 0); i++) {
                NodeBlock block = nodes.block[i];
                switch (SIZE) {
                    case U32:
                        fprintf(stream, "\t %d %d %d %" PRIu32 "\n", block.entity_dim,
                                block.entity_tag, block.parametric, block.num_nodes.u32);
                        for (long j = 0; j < (long)block.num_nodes.u32; j++) {
                            fprintf(stream, "\t %" PRIu32 "\n", block.tag.u32[j]);
                        }
                        break;
                    case U64:
                        fprintf(stream, "\t %d %d %d %" PRIu64 "\n", block.entity_dim,
                                block.entity_tag, block.parametric, block.num_nodes.u64);
                        for (long j = 0; j < (long)block.num_nodes.u64; j++) {
                            fprintf(stream, "\t %" PRIu64 "\n", block.tag.u64[j]);
                        }
                        break;
                    default: assert(false);
                }
                long len = 3 + (block.parametric ? block.entity_dim : 0);
                for (long j = 0; j < (long)block.num_nodes.u64; j++) {
                    fprintf(stream, "\t");
                    for (long k = 0; k < len; k++) {
                        fprintf(stream, " %g", block.coord[(j * len) + k]);
                    }
                    fprintf(stream, "\n");
                }
            }
            fprintf(stream, "\t $EndNodes\n");

            Elements elements = gmsh->elements;
            fprintf(stream, "\t $Elements\n");
            switch (SIZE) {
                case U32:
                    fprintf(stream, "\t %" PRIu32 " %" PRIu32 " %" PRIu32 " %" PRIu32 "\n",
                            elements.num_blocks.u32, elements.tot_elements.u32,
                            elements.min_tag.u32, elements.max_tag.u32);
                    break;
                case U64:
                    fprintf(stream, "\t %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 "\n",
                            elements.num_blocks.u64, elements.tot_elements.u64,
                            elements.min_tag.u64, elements.max_tag.u64);
                    break;
                default: assert(false);
            }
            for (long i = 0; i < parse_data_to_long(SIZE, &elements.num_blocks, 0); i++) {
                ElementBlock block = elements.block[i];
                long len = num_node_tags(block.element_type);
                switch (SIZE) {
                    case U32:
                        fprintf(stream, "\t %d %d %d %" PRIu32 "\n", block.entity_dim,
                                block.entity_tag, block.element_type, block.num_elements.u32);
                        for (long j = 0; j < (long)block.num_elements.u32; j++) {
                            fprintf(stream, "\t %" PRIu32, block.tag.u32[j]);
                            for (long k = 0; k < len; k++) {
                                fprintf(stream, " %" PRIu32, block.node_tag.u32[(j * len) + k]);
                            }
                            fprintf(stream, "\n");
                        }
                        break;
                    case U64:
                        fprintf(stream, "\t %d %d %d %" PRIu64 "\n", block.entity_dim,
                                block.entity_tag, block.element_type, block.num_elements.u64);
                        for (long j = 0; j < (long)block.num_elements.u64; j++) {
                            fprintf(stream, "\t %" PRIu64, block.tag.u64[j]);
                            for (long k = 0; k < len; k++) {
                                fprintf(stream, " %" PRIu64, block.node_tag.u64[(j * len) + k]);
                            }
                            fprintf(stream, "\n");
                        }
                        break;
                    default: assert(false);
                }
            }
            fprintf(stream, "\t $EndElements\n");

            fflush(stream);
            off = ftell(stream);
        }
        MPI_Barrier(sync.comm);
        MPI_Bcast(&off, 1, MPI_LONG, rank, sync.comm);
    }
}

static void create_nodes(MeshNodes *nodes, Type SIZE, const Gmsh *gmsh)
{
    long num_blocks = parse_data_to_long(SIZE, &gmsh->nodes.num_blocks, 0);
    long num = 0;
    for (long i = 0; i < num_blocks; i++) {
        num += parse_data_to_long(SIZE, &gmsh->nodes.block[i].num_nodes, 0);
    }
    nodes->num = num;
    assert(nodes->num > 0);

    nodes->coord = malloc(nodes->num * sizeof(*nodes->coord));
    assert(nodes->coord);

    num = 0;
    for (long i = 0; i < num_blocks; i++) {
        NodeBlock block = gmsh->nodes.block[i];
        long num_nodes = parse_data_to_long(SIZE, &block.num_nodes, 0);
        long len = 3 + (block.parametric ? block.entity_dim : 0);
        for (long j = 0; j < num_nodes; j++) {
            nodes->coord[num].x = block.coord[(j * len) + 0];
            nodes->coord[num].y = block.coord[(j * len) + 1];
            nodes->coord[num].z = block.coord[(j * len) + 2];
            num += 1;
        }
    }
    assert(num == nodes->num);
}

typedef struct {
    long tag;
    long idx;
} Map;

static int cmp_map(const void *lhs, const void *rhs)
{
    const Map *_lhs = lhs;
    const Map *_rhs = rhs;
    return cmp_asc(_lhs->tag, _rhs->tag);
}

static void convert_node_tags_to_indices(const MeshNodes *nodes, MeshCells *cells, Type SIZE,
                                         const Gmsh *gmsh)
{
    Arena save = arena_save();

    long cap = sync_lmax(nodes->num);
    Map *map = arena_calloc(cap, sizeof(*map));

    long off_nodes = sync_lexsum(nodes->num);
    long num = 0;
    for (long i = 0; i < parse_data_to_long(SIZE, &gmsh->nodes.num_blocks, 0); i++) {
        NodeBlock block = gmsh->nodes.block[i];
        for (long j = 0; j < parse_data_to_long(SIZE, &block.num_nodes, 0); j++) {
            switch (SIZE) {
                case U32: map[num].tag = block.tag.u32[j]; break;
                case U64: map[num].tag = block.tag.u64[j]; break;
                default: assert(false);
            }
            map[num].idx = num + off_nodes;
            num += 1;
        }
    }
    assert(num == nodes->num);
    qsort(map, num, sizeof(*map), cmp_map);

    long num_tags = cells->node.off[cells->num];
    bool *converted = arena_calloc(num_tags, sizeof(*converted));

    MPI_Datatype type;
    MPI_Type_contiguous(sizeof(Map), MPI_BYTE, &type);
    MPI_Type_commit(&type);

    int dst = (sync.rank + 1) % sync.size;
    int src = (sync.rank - 1 + sync.size) % sync.size;
    int tag = sync_tag();
    for (long step = 0; step < sync.size; step++) {
        for (long i = 0; i < num_tags; i++) {
            if (!converted[i]) {
                Map key = {.tag = cells->node.idx[i]};
                Map *val = bsearch(&key, map, num, sizeof(*map), cmp_map);
                if (val) {
                    cells->node.idx[i] = val->idx;
                    converted[i] = true;
                }
            }
        }
        MPI_Sendrecv_replace(&num, 1, MPI_LONG, dst, tag, src, tag, sync.comm, MPI_STATUS_IGNORE);
        MPI_Sendrecv_replace(map, cap, type, dst, tag, src, tag, sync.comm, MPI_STATUS_IGNORE);
    }
    MPI_Type_free(&type);

    for (long i = 0; i < num_tags; i++) {
        assert(converted[i]);
    }

    arena_load(save);
}

static void create_cells(const MeshNodes *nodes, MeshCells *cells, Type SIZE, const Gmsh *gmsh)
{
    long num_blocks = parse_data_to_long(SIZE, &gmsh->elements.num_blocks, 0);
    long num = 0;
    for (long i = 0; i < num_blocks; i++) {
        num += parse_data_to_long(SIZE, &gmsh->elements.block[i].num_elements, 0);
    }
    cells->num = num;
    assert(cells->num > 0);

    long *off = malloc((cells->num + 1) * sizeof(*off));
    assert(off);

    long *idx = malloc((cells->num * MAX_CELL_NODES) * sizeof(*idx));
    assert(idx);

    off[0] = 0;

    num = 0;
    for (long i = 0; i < num_blocks; i++) {
        ElementBlock block = gmsh->elements.block[i];
        long num_elements = parse_data_to_long(SIZE, &block.num_elements, 0);
        long len = num_node_tags(block.element_type);
        for (long j = 0; j < num_elements; j++) {
            off[num + 1] = off[num] + len;
            for (long k = 0; k < len; k++) {
                switch (SIZE) {
                    case U32: idx[off[num] + k] = block.node_tag.u32[(j * len) + k]; break;
                    case U64: idx[off[num] + k] = block.node_tag.u64[(j * len) + k]; break;
                    default: assert(false);
                }
            }
            num += 1;
        }
    }
    assert(num == cells->num);

    cells->node.off = off;
    cells->node.idx = realloc(idx, off[cells->num] * sizeof(*idx));
    assert(cells->node.idx);

    convert_node_tags_to_indices(nodes, cells, SIZE, gmsh);
}

static int cmp_physical(const void *lhs_, const void *rhs_)
{
    const Physical *lhs = lhs_;
    const Physical *rhs = rhs_;
    if (lhs->dim != rhs->dim) {
        return cmp_dsc(lhs->dim, rhs->dim);
    }
    bool lhs_periodic = strchr(lhs->name, ':');
    bool rhs_periodic = strchr(rhs->name, ':');
    if (lhs_periodic != rhs_periodic) {
        return cmp_asc(lhs_periodic, rhs_periodic);
    }
    return cmp_asc(lhs->tag, rhs->tag);
}

static void create_entities(MeshEntities *entities, const Gmsh *gmsh)
{
    Arena save = arena_save();

    long num = gmsh->physicals.num;
    Physical *physical = arena_memdup(gmsh->physicals.physical, num, sizeof(*physical));
    qsort(physical, num, sizeof(*physical), cmp_physical);

    Name *name = calloc(num, sizeof(*name));
    assert(name);

    long num_inner = 0;
    long num_ghost = 0;
    for (long i = 0; i < num; i++) {
        strcpy(name[i], physical[i].name);
        if (physical[i].dim == 3) {
            num_inner += 1;
        }
        else {
            assert(physical[i].dim == 2);
            if (!strchr(name[i], ':')) {
                num_ghost += 1;
            }
        }
    }

    arena_load(save);

    entities->num = num;
    entities->num_inner = num_inner;
    entities->num_ghost = num_ghost;
    entities->name = name;
}

static int cmp_surface(const void *lhs_, const void *rhs_)
{
    const Surface *lhs = lhs_;
    const Surface *rhs = rhs_;
    return cmp_asc(lhs->tag, rhs->tag);
}

static int cmp_volume(const void *lhs_, const void *rhs_)
{
    const Volume *lhs = lhs_;
    const Volume *rhs = rhs_;
    return cmp_asc(lhs->tag, rhs->tag);
}

static void compute_entity_map(long *entity, const MeshCells *cells, Type SIZE, const Gmsh *gmsh)
{
    Arena save = arena_save();

    long num_physicals = gmsh->physicals.num;
    Physical *physical = arena_memdup(gmsh->physicals.physical, num_physicals, sizeof(*physical));
    qsort(physical, num_physicals, sizeof(*physical), cmp_physical);

    long num_surfaces = parse_data_to_long(SIZE, &gmsh->entities.num_surfaces, 0);
    Surface *surface = arena_memdup(gmsh->entities.surface, num_surfaces, sizeof(*surface));
    qsort(surface, num_surfaces, sizeof(*surface), cmp_surface);

    long num_volumes = parse_data_to_long(SIZE, &gmsh->entities.num_volumes, 0);
    Volume *volume = arena_memdup(gmsh->entities.volume, num_volumes, sizeof(*volume));
    qsort(volume, num_volumes, sizeof(*volume), cmp_volume);

    long num = 0;
    for (long i = 0; i < parse_data_to_long(SIZE, &gmsh->elements.num_blocks, 0); i++) {
        ElementBlock block = gmsh->elements.block[i];
        long physical_tag;
        switch (block.entity_dim) {
            case 2: {
                Surface key = {.tag = block.entity_tag};
                Surface *val = bsearch(&key, surface, num_surfaces, sizeof(*surface), cmp_surface);
                assert(val);
                long num_physical_tags = parse_data_to_long(SIZE, &val->num_physical_tags, 0);
                assert(num_physical_tags == 1);
                physical_tag = val->physical_tag[0];
                break;
            }
            case 3: {
                Volume key = {.tag = block.entity_tag};
                Volume *val = bsearch(&key, volume, num_volumes, sizeof(*volume), cmp_volume);
                assert(val);
                long num_physical_tags = parse_data_to_long(SIZE, &val->num_physical_tags, 0);
                assert(num_physical_tags == 1);
                physical_tag = val->physical_tag[0];
                break;
            }
            default: assert(false);
        }
        long idx = -1;
        for (long j = 0; j < num_physicals; j++) {
            if (physical[j].dim == block.entity_dim && physical[j].tag == physical_tag) {
                idx = j;
                break;
            }
        }
        assert(idx != -1);
        for (long j = 0; j < parse_data_to_long(SIZE, &block.num_elements, 0); j++) {
            entity[num++] = idx;
        }
    }
    assert(num == cells->num);

    arena_load(save);
}

static void reorder(MeshCells *cells, MeshEntities *entities, Type SIZE, const Gmsh *gmsh)
{
    Arena save = arena_save();

    long *entity = arena_malloc(cells->num, sizeof(*entity));
    compute_entity_map(entity, cells, SIZE, gmsh);

    long *num_cells = arena_calloc(entities->num, sizeof(*num_cells));
    for (long i = 0; i < cells->num; i++) {
        num_cells[entity[i]] += 1;
    }

    long *cell_off = malloc((entities->num + 1) * sizeof(*cell_off));
    assert(cell_off);

    cell_off[0] = 0;
    for (long i = 0; i < entities->num; i++) {
        cell_off[i + 1] = cell_off[i] + num_cells[i];
    }

    long *map = arena_malloc(cells->num, sizeof(*map));
    for (long i = 0; i < entities->num; i++) {
        cell_off[i + 1] -= num_cells[i];
    }
    for (long i = 0; i < cells->num; i++) {
        map[i] = cell_off[entity[i] + 1]++;
    }
    mesh_reorder_cells(cells, 0, 0, cells->num, map);

    arena_load(save);

    entities->cell_off = cell_off;
}

void mesh_read_gmsh(Mesh *mesh, const char *fname)
{
    Arena save = arena_save();

    Gmsh *gmsh = arena_calloc(1, sizeof(*gmsh));

    File file = parse_open(fname);

    double version;
    Mode mode;
    Type SIZE;
    bool swap = false;
    read_format(&version, &mode, &SIZE, &swap, file);
    assert(isclose(version, 4.1));  // NOLINT(readability-magic-numbers)

    read_physicals(&gmsh->physicals, file);
    read_entities(&gmsh->entities, mode, SIZE, swap, file);
    read_nodes(&gmsh->nodes, mode, SIZE, swap, file);
    read_elements(&gmsh->elements, mode, SIZE, swap, file);

    parse_close(file);

    // gmsh_dump(stdout, version, mode, SIZE, swap, gmsh);

    create_nodes(&mesh->nodes, SIZE, gmsh);
    create_cells(&mesh->nodes, &mesh->cells, SIZE, gmsh);
    create_entities(&mesh->entities, gmsh);

    reorder(&mesh->cells, &mesh->entities, SIZE, gmsh);

    arena_load(save);
}
