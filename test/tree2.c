#include "tree2.h"

#include "teal2.h"
#include "test.h"

#define EPS 1e-8

static void test_empty(void)
{
    Tree2 *tree = tree2_init(0, 0, 0);

    Vector point = {0, 0, 0};

    int idx[4];
    TEST(tree2_nearest(tree, point, idx, 4) == 0);
    TEST(tree2_radius(tree, point, 1.0, idx, 4) == 0);

    tree2_deinit(tree);
}

static void test_single(void)
{
    Vector point[] = {
        {1, 2, 3},
    };
    Tree2 *tree = tree2_init(point, 1, 0);

    int idx[2];
    TEST(tree2_nearest(tree, point[0], idx, 2) == 1);
    TEST(idx[0] == 0);

    TEST(tree2_radius(tree, point[0], EPS, idx, 2) == 1);
    TEST(idx[0] == 0);

    tree2_deinit(tree);
}

static double dist(Vector lhs, Vector rhs)
{
    return vector2_norm(vector2_sub(lhs, rhs));
}

static void test_nearest_order_and_cap(void)
{
    Vector point[] = {
        {1, 0, 0},
        {2, 0, 0},
        {3, 0, 0},
        {4, 0, 0},
    };
    Tree2 *tree = tree2_init(point, 4, 0);

    Vector pos = {0, 0, 0};

    int idx[3];
    TEST(tree2_nearest(tree, pos, idx, 3) == 3);
    TEST(idx[0] == 0);
    TEST(idx[1] == 1);
    TEST(idx[2] == 2);

    for (int i = 1; i < 3; i++) {
        TEST(dist(point[idx[i - 1]], pos) <= dist(point[idx[i]], pos));
    }

    tree2_deinit(tree);
}

static void test_radius_count_and_truncation(void)
{
    Vector point[] = {
        {1, 0, 0},
        {2, 0, 0},
        {3, 0, 0},
        {4, 0, 0},
    };
    Tree2 *tree = tree2_init(point, 4, 0);

    Vector pos = {0, 0, 0};

    int idx[2];
    TEST(tree2_radius(tree, pos, 3, idx, 2) == 3);

    for (int i = 0; i < 2; i++) {
        TEST(0 <= idx[i] && idx[i] < 4);
        TEST(dist(point[idx[i]], pos) <= 3);
    }

    tree2_deinit(tree);
}

static int has_idx(const int *idx, int num, int val)
{
    for (int i = 0; i < num; i++) {
        if (idx[i] == val) {
            return 1;
        }
    }
    return 0;
}

static void test_duplicates(void)
{
    Vector point[] = {
        {5, 5, 5},
        {5, 5, 5},
        {8, 8, 8},
    };
    Tree2 *tree = tree2_init(point, 3, 0);

    Vector pos = {5, 5, 5};

    int idx[4];
    TEST(tree2_radius(tree, pos, EPS, idx, 4) == 2);
    TEST(has_idx(idx, 2, 0));
    TEST(has_idx(idx, 2, 1));

    TEST(tree2_nearest(tree, pos, idx, 2) == 2);
    TEST(has_idx(idx, 2, 0));
    TEST(has_idx(idx, 2, 1));

    tree2_deinit(tree);
}

int main(int argc, char **argv)
{
    teal2_init(&argc, &argv);

    test_empty();
    test_single();
    test_nearest_order_and_cap();
    test_radius_count_and_truncation();
    test_duplicates();

    teal2_deinit();
}
