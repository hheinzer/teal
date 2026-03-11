#include <assert.h>

#include "defer2.h"
#include "teal2.h"

typedef struct item {
    void *ptr;
    void (*func)(void *);
    struct item *prev;
} Item;

static Item *stack;

void defer2_init(void)
{
    stack = 0;
}

void defer2_deinit(void)
{
    Item *head = stack;
    while (head) {
        Item *prev = head->prev;
        if (head->func) {
            head->func(head->ptr);
        }
        teal2_free(head);
        head = prev;
    }
}

void *defer2_push(void *ptr, void (*func)(void *))
{
    assert(ptr && func);
    Item *next = teal2_calloc(1, sizeof(*next));
    next->ptr = ptr;
    next->func = func;
    next->prev = stack;
    stack = next;
    return ptr;
}

void *defer2_pop(void *ptr)
{
    assert(ptr);
    for (Item *item = stack; item; item = item->prev) {
        if (item->ptr == ptr) {
            item->func = 0;
            break;
        }
    }
    return ptr;
}
