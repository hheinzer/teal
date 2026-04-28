#include "kdtree.h"

#include "../test.h"
#include "teal.h"

static void test_empty(void)
{
    Kdtree *tree = kdtree_init(0, 0);

    Vector point = {0, 0, 0};

    int idx[4];
    test(kdtree_nearest(tree, point, idx, 4) == 0);

    kdtree_deinit(tree);
}

static void test_single(void)
{
    Vector point[] = {
        {1, 2, 3},
    };
    Kdtree *tree = kdtree_init(point, 1);

    int idx[2];
    test(kdtree_nearest(tree, point[0], idx, 2) == 1);
    test(idx[0] == 0);

    kdtree_deinit(tree);
}

static double dist(Vector lhs, Vector rhs)
{
    return vector_norm(vector_sub(lhs, rhs));
}

static void test_nearest_order_and_cap(void)
{
    Vector point[] = {
        {1, 0, 0},
        {2, 0, 0},
        {3, 0, 0},
        {4, 0, 0},
    };
    Kdtree *tree = kdtree_init(point, 4);

    Vector pos = {0, 0, 0};

    int idx[3];
    test(kdtree_nearest(tree, pos, idx, 3) == 3);
    test(idx[0] == 0);
    test(idx[1] == 1);
    test(idx[2] == 2);

    for (int i = 1; i < 3; i++) {
        test(dist(point[idx[i - 1]], pos) <= dist(point[idx[i]], pos));
    }

    kdtree_deinit(tree);
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
    Kdtree *tree = kdtree_init(point, 3);

    Vector pos = {5, 5, 5};

    int idx[4];
    test(kdtree_nearest(tree, pos, idx, 2) == 2);
    test(has_idx(idx, 2, 0));
    test(has_idx(idx, 2, 1));

    kdtree_deinit(tree);
}

int main(int argc, char **argv)
{
    teal_init(&argc, &argv);

    test_empty();
    test_single();
    test_nearest_order_and_cap();
    test_duplicates();

    teal_deinit();
}
