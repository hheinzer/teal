#include <stdio.h>

#include "mesh.h"
#include "teal/sync.h"

static void print_null(void)
{
    printf(" - |");
}

static void print_long(const long *arr, long idx)
{
    if (arr) {
        printf(" %ld |", arr[idx]);
    }
    else {
        print_null();
    }
}

static void print_long_pair(const long *arr, long idx)
{
    if (arr) {
        printf(" %ld %ld |", arr[idx], arr[idx + 1]);
    }
    else {
        print_null();
    }
}

static void print_double(const double *arr, long idx)
{
    if (arr) {
        printf(" %g |", arr[idx]);
    }
    else {
        print_null();
    }
}

static void print_string(const string *arr, long idx)
{
    if (arr) {
        printf(" %s |", arr[idx]);
    }
    else {
        print_null();
    }
}

static void print_vector(const vector *arr, long idx)
{
    if (arr) {
        vector vec = arr[idx];
        printf(" %g %g %g |", vec.x, vec.y, vec.z);
    }
    else {
        print_null();
    }
}

static void print_matrix(const matrix *arr, long idx)
{
    if (arr) {
        matrix mat = arr[idx];
        printf(" %g %g %g  %g %g %g  %g %g %g |", mat.x.x, mat.x.y, mat.x.z, mat.y.x, mat.y.y,
               mat.y.z, mat.z.x, mat.z.y, mat.z.z);
    }
    else {
        print_null();
    }
}

static void print_graph(const MeshGraph *graph, long idx)
{
    if (graph->off && graph->idx) {
        for (long j = graph->off[idx]; j < graph->off[idx + 1]; j++) {
            printf(" %ld", graph->idx[j]);
        }
        printf(" |");
    }
    else {
        print_null();
    }
}

static void print_face_cell(const MeshFaceCell *cell, long idx)
{
    if (cell) {
        printf(" %ld %ld |", cell[idx].left, cell[idx].right);
    }
    else {
        print_null();
    }
}

void mesh_print(const Mesh *mesh)
{
    for (long rank = 0; rank < sync.size; rank++) {
        if (rank == sync.rank) {
            printf("rank %d:\n", sync.rank);

            if (mesh->nodes.num) {
                printf("\t mesh->nodes.{num|inner} = %ld | %ld\n", mesh->nodes.num,
                       mesh->nodes.num_inner);
                printf("\t mesh->nodes.{global|coord} = [\n");
                for (long i = 0; i < mesh->nodes.num; i++) {
                    printf("\t\t [%ld]", i);
                    print_long(mesh->nodes.global, i);
                    print_vector(mesh->nodes.coord, i);
                    printf("\n");
                }
                printf("\t ]\n");
            }

            if (mesh->cells.num) {
                printf("\t mesh->cells.{num|inner|ghost|periodic} = %ld | %ld | %ld | %ld\n",
                       mesh->cells.num, mesh->cells.num_inner, mesh->cells.num_ghost,
                       mesh->cells.num_periodic);
                printf("\t mesh->cells.{node|cell|volume|center|projection} = [\n");
                for (long i = 0; i < mesh->cells.num; i++) {
                    printf("\t\t [%ld]", i);
                    print_graph(&mesh->cells.node, i);
                    print_graph(&mesh->cells.cell, i);
                    print_double(mesh->cells.volume, i);
                    print_vector(mesh->cells.center, i);
                    print_vector(mesh->cells.projection, i);
                    printf("\n");
                }
                printf("\t ]\n");
            }

            if (mesh->faces.num) {
                printf("\t mesh->faces.{num|inner|ghost} = %ld | %ld | %ld\n", mesh->faces.num,
                       mesh->faces.num_inner, mesh->faces.num_ghost);
                printf("\t mesh->faces.{node|cell|area|center|basis|weight} = [\n");
                for (long i = 0; i < mesh->faces.num; i++) {
                    printf("\t\t [%ld]", i);
                    print_graph(&mesh->faces.node, i);
                    print_face_cell(mesh->faces.cell, i);
                    print_double(mesh->faces.area, i);
                    print_vector(mesh->faces.center, i);
                    print_matrix(mesh->faces.basis, i);
                    print_vector(mesh->faces.weight, i);
                    printf("\n");
                }
                printf("\t ]\n");
            }

            if (mesh->entities.num) {
                printf("\t mesh->entities.{num|inner|ghost} = %ld | %ld | %ld\n",
                       mesh->entities.num, mesh->entities.num_inner, mesh->entities.num_ghost);
                printf("\t mesh->entities.{name|cell|face|offset} = [\n");
                for (long i = 0; i < mesh->entities.num; i++) {
                    printf("\t\t [%ld]", i);
                    print_string(mesh->entities.name, i);
                    print_long_pair(mesh->entities.cell_off, i);
                    print_long_pair(mesh->entities.face_off, i);
                    print_vector(mesh->entities.offset, i);
                    printf("\n");
                }
                printf("\t ]\n");
            }

            if (mesh->neighbors.num) {
                printf("\t mesh->neighbors.num = %ld\n", mesh->neighbors.num);
                printf("\t mesh->neighbors.{rank|recv|send} = [\n");
                for (long i = 0; i < mesh->neighbors.num; i++) {
                    printf("\t\t [%ld]", i);
                    print_long(mesh->neighbors.rank, i);
                    print_long_pair(mesh->neighbors.recv_off, i);
                    print_graph(&mesh->neighbors.send, i);
                    printf("\n");
                }
                printf("\t ]\n");
            }
        }
        MPI_Barrier(sync.comm);
    }
    MPI_Barrier(sync.comm);
}
