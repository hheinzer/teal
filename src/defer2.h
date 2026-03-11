#pragma once

// Initialize the defer stack.
void defer2_init(void);

// Run all deferred calls in LIFO order.
void defer2_deinit(void);

// Push a deferred call onto the stack. Return `ptr`.
void *defer2_push(void *ptr, void (*func)(void *));

// Cancel the deferred call associated with `ptr` without running it. Return `ptr`.
void *defer2_pop(void *ptr);
