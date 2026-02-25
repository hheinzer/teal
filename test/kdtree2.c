#include "kdtree2.h"

#include "teal2.h"
#include "test.h"
#include "utils2.h"

typedef struct {
    int id;
    Vector pos;
} Data;

static void test_empty(void)
{
    Kdtree2 *tree = kdtree2_init(sizeof(Data));

    Vector pos = {0, 0, 0};

    Data out[4];
    TEST(kdtree2_nearest(tree, pos, out, 4) == 0);
    TEST(kdtree2_radius(tree, pos, out, 4, 1.0) == 0);

    kdtree2_deinit(tree);
}

static double dist(Vector lhs, Vector rhs)
{
    return vector2_norm(vector2_sub(lhs, rhs));
}

static void test_single_and_copy(void)
{
    Kdtree2 *tree = kdtree2_init(sizeof(Data));

    Vector pos = {1, 2, 3};

    Data src = {42, pos};
    kdtree2_insert(tree, src.pos, &src);

    src.id = 24;
    src.pos = (Vector){3, 2, 1};

    Data out[2];
    TEST(kdtree2_nearest(tree, pos, out, 2) == 1);
    TEST(out[0].id == 42);
    TEST(isclose(dist(out[0].pos, pos), 0));

    TEST(kdtree2_radius(tree, pos, out, 2, 1e-8) == 1);
    TEST(out[0].id == 42);
    TEST(isclose(dist(out[0].pos, pos), 0));

    kdtree2_deinit(tree);
}

static void test_nearest_order_and_cap(void)
{
    Kdtree2 *tree = kdtree2_init(sizeof(Data));

    Data src[4] = {
        {1, {1, 0, 0}},
        {2, {2, 0, 0}},
        {3, {3, 0, 0}},
        {4, {4, 0, 0}},
    };
    for (int i = 0; i < 4; i++) {
        kdtree2_insert(tree, src[i].pos, &src[i]);
    }

    Vector pos = {0, 0, 0};

    Data out[3];
    TEST(kdtree2_nearest(tree, pos, out, 3) == 3);
    TEST(out[0].id == 1);
    TEST(out[1].id == 2);
    TEST(out[2].id == 3);
    for (int i = 1; i < 3; i++) {
        TEST(dist(out[i - 1].pos, pos) <= dist(out[i].pos, pos));
    }

    kdtree2_deinit(tree);
}

static int has_id(const Data *data, int num, int id)
{
    for (int i = 0; i < num; i++) {
        if (data[i].id == id) {
            return 1;
        }
    }
    return 0;
}

static void test_radius_count_and_truncation(void)
{
    Kdtree2 *tree = kdtree2_init(sizeof(Data));

    Data src[4] = {
        {1, {1, 0, 0}},
        {2, {2, 0, 0}},
        {3, {3, 0, 0}},
        {4, {4, 0, 0}},
    };
    for (int i = 0; i < 4; i++) {
        kdtree2_insert(tree, src[i].pos, &src[i]);
    }

    Vector pos = {0, 0, 0};

    Data out[2];
    TEST(kdtree2_radius(tree, pos, out, 2, 3) == 3);
    for (int i = 0; i < 2; i++) {
        TEST(has_id(src, 4, out[i].id));
        TEST(dist(out[i].pos, pos) <= 3);
    }

    kdtree2_deinit(tree);
}

static void test_duplicates(void)
{
    Kdtree2 *tree = kdtree2_init(sizeof(Data));

    Vector pos = {5, 5, 5};

    Data a = {21, pos};
    kdtree2_insert(tree, a.pos, &a);

    Data b = {22, pos};
    kdtree2_insert(tree, b.pos, &b);

    Data c = {23, {8, 8, 8}};
    kdtree2_insert(tree, c.pos, &c);

    Data out[4];
    TEST(kdtree2_radius(tree, pos, out, 4, 1e-8) == 2);
    TEST(has_id(out, 2, 21));
    TEST(has_id(out, 2, 22));

    TEST(kdtree2_nearest(tree, pos, out, 2) == 2);
    TEST(has_id(out, 2, 21));
    TEST(has_id(out, 2, 22));

    kdtree2_deinit(tree);
}

int main(int argc, char **argv)
{
    teal2_init(&argc, &argv);

    test_empty();
    test_single_and_copy();
    test_nearest_order_and_cap();
    test_radius_count_and_truncation();
    test_duplicates();

    teal2_deinit();
}
