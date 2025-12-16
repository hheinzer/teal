#include <string.h>

#include "mesh.h"
#include "teal/option.h"
#include "teal/utils.h"

int main(int argc, char **argv)
{
    option.quiet = true;
    teal_initialize(&argc, &argv);

    if (argc < 2) {
        error("usage: %s <input> [prefix]", argv[0]);
    }

    const char *input = argv[1];

    char prefix[128];
    if (argc > 2) {
        strcpy(prefix, argv[2]);
    }
    else {
        strcpy(prefix, input);
        char *dot = strrchr(prefix, '.');
        if (dot) {
            *dot = 0;
        }
    }

    Mesh *mesh = mesh_read(input);
    mesh_write(mesh, prefix);

    teal_finalize();
}
