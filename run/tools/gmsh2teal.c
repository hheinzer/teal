#include <stdio.h>
#include <string.h>

#include "mesh.h"
#include "teal/option.h"
#include "teal/sync.h"
#include "teal/utils.h"

int main(int argc, char **argv)
{
    option.quiet = true;
    teal_initialize(&argc, &argv);

    if (argc < 2) {
        error("usage: %s <input> [output]", argv[0]);
    }

    const char *input = argv[1];
    Mesh *mesh = mesh_read(input);
    mesh_write(mesh, input);

    char output[128];
    if (argc > 2) {
        strcpy(output, argv[2]);
    }
    else {
        strcpy(output, input);
        char *dot = strrchr(output, '.');
        if (dot) {
            strcpy(dot, ".vtkhdf");
        }
        else {
            strcat(output, ".vtkhdf");
        }
    }

    if (sync.rank == 0) {
        char fname[128];
        sprintf(fname, "%s_mesh.vtkhdf", input);
        if (rename(fname, output) != 0) {
            error("could not rename (%s) to (%s)", fname, output);
        }
        printf("%s -> %s\n", input, output);
    }

    teal_finalize();
}
