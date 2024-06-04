#include "hash.h"

uint64_t hash_fnv_64a(const void *ptr, const long nmemb, const long size)
{
    // FNV-1a 64 bit: https://en.wikipedia.org/wiki/Fowler%E2%80%93Noll%E2%80%93Vo_hash_function
    const unsigned char *byte = (typeof(byte))ptr;
    const unsigned char *byte_end = byte + nmemb * size;
    uint64_t hash = 0xcbf29ce484222325;
    while (byte < byte_end) {
        hash ^= (typeof(hash))*byte++;
        hash *= 0x00000100000001b3;
    }
    return hash;
}
