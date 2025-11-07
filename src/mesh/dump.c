#include <stdio.h>

#include "mesh.h"
#include "teal/assert.h"
#include "teal/sync.h"

static void print_null(FILE *stream)
{
    fprintf(stream, " - |");
}

static void print_number(FILE *stream, const number *arr, number idx)
{
    if (arr) {
        fprintf(stream, " %td |", arr[idx]);
    }
    else {
        print_null(stream);
    }
}

static void print_number_pair(FILE *stream, const number *arr, number idx)
{
    if (arr) {
        fprintf(stream, " %td %td |", arr[idx], arr[idx + 1]);
    }
    else {
        print_null(stream);
    }
}

static void print_scalar(FILE *stream, const scalar *arr, number idx)
{
    if (arr) {
        fprintf(stream, " %g |", arr[idx]);
    }
    else {
        print_null(stream);
    }
}

static void print_name(FILE *stream, const Name *arr, number idx)
{
    if (arr) {
        fprintf(stream, " %s |", arr[idx]);
    }
    else {
        print_null(stream);
    }
}

static void print_vector(FILE *stream, const vector *arr, number idx)
{
    if (arr) {
        vector vec = arr[idx];
        fprintf(stream, " %g %g %g |", vec.x, vec.y, vec.z);
    }
    else {
        print_null(stream);
    }
}

static void print_matrix(FILE *stream, const matrix *arr, number idx)
{
    if (arr) {
        matrix mat = arr[idx];
        fprintf(stream, " %g %g %g  %g %g %g  %g %g %g |", mat.x.x, mat.x.y, mat.x.z, mat.y.x,
                mat.y.y, mat.y.z, mat.z.x, mat.z.y, mat.z.z);
    }
    else {
        print_null(stream);
    }
}

static void print_graph(FILE *stream, Graph graph, number idx)
{
    if (graph.off && graph.idx) {
        for (number j = graph.off[idx]; j < graph.off[idx + 1]; j++) {
            fprintf(stream, " %td", graph.idx[j]);
        }
        fprintf(stream, " |");
    }
    else {
        print_null(stream);
    }
}

static void print_adjacent(FILE *stream, const Adjacent *cell, number idx)
{
    if (cell) {
        fprintf(stream, " %td %td |", cell[idx].left, cell[idx].right);
    }
    else {
        print_null(stream);
    }
}

void mesh_dump(FILE *stream, const Mesh *mesh)
{
    assert(stream && mesh);
    for (number rank = 0; rank < sync.size; rank++) {
        if (rank == sync.rank) {
            fprintf(stream, "rank %d:\n", sync.rank);

            if (mesh->nodes.num) {
                fprintf(stream, "\t mesh->nodes.{num|inner} = %td | %td\n", mesh->nodes.num,
                        mesh->nodes.num_inner);
                fprintf(stream, "\t mesh->nodes.{global|coord} = [\n");
                for (number i = 0; i < mesh->nodes.num; i++) {
                    fprintf(stream, "\t\t [%td]", i);
                    print_number(stream, mesh->nodes.global, i);
                    print_vector(stream, mesh->nodes.coord, i);
                    fprintf(stream, "\n");
                }
                fprintf(stream, "\t ]\n");
            }

            if (mesh->cells.num) {
                fprintf(stream,
                        "\t mesh->cells.{num|num_inner|off_ghost|off_periodic} ="
                        " %td | %td | %td | %td\n",
                        mesh->cells.num, mesh->cells.num_inner, mesh->cells.off_ghost,
                        mesh->cells.off_periodic);
                fprintf(stream, "\t mesh->cells.{node|cell|volume|center|projection} = [\n");
                for (number i = 0; i < mesh->cells.num; i++) {
                    fprintf(stream, "\t\t [%td]", i);
                    print_graph(stream, mesh->cells.node, i);
                    print_graph(stream, mesh->cells.cell, i);
                    print_scalar(stream, mesh->cells.volume, i);
                    print_vector(stream, mesh->cells.center, i);
                    print_vector(stream, mesh->cells.projection, i);
                    fprintf(stream, "\n");
                }
                fprintf(stream, "\t ]\n");
            }

            if (mesh->faces.num) {
                fprintf(stream, "\t mesh->faces.{num|num_inner|off_ghost} = %td | %td | %td\n",
                        mesh->faces.num, mesh->faces.num_inner, mesh->faces.off_ghost);
                fprintf(stream, "\t mesh->faces.{node|cell|area|center|basis|weight} = [\n");
                for (number i = 0; i < mesh->faces.num; i++) {
                    fprintf(stream, "\t\t [%td]", i);
                    print_graph(stream, mesh->faces.node, i);
                    print_adjacent(stream, mesh->faces.cell, i);
                    print_scalar(stream, mesh->faces.area, i);
                    print_vector(stream, mesh->faces.center, i);
                    print_matrix(stream, mesh->faces.basis, i);
                    print_vector(stream, mesh->faces.weight, i);
                    fprintf(stream, "\n");
                }
                fprintf(stream, "\t ]\n");
            }

            if (mesh->entities.num) {
                fprintf(stream, "\t mesh->entities.{num|num_inner|off_ghost} = %td | %td | %td\n",
                        mesh->entities.num, mesh->entities.num_inner, mesh->entities.off_ghost);
                fprintf(stream, "\t mesh->entities.{name|cell|face|translation} = [\n");
                for (number i = 0; i < mesh->entities.num; i++) {
                    fprintf(stream, "\t\t [%td]", i);
                    print_name(stream, mesh->entities.name, i);
                    print_number_pair(stream, mesh->entities.cell_off, i);
                    print_number_pair(stream, mesh->entities.face_off, i);
                    print_vector(stream, mesh->entities.translation, i);
                    fprintf(stream, "\n");
                }
                fprintf(stream, "\t ]\n");
            }

            if (mesh->neighbors.num) {
                fprintf(stream, "\t mesh->neighbors.num = %td\n", mesh->neighbors.num);
                fprintf(stream, "\t mesh->neighbors.{rank|recv|send} = [\n");
                for (number i = 0; i < mesh->neighbors.num; i++) {
                    fprintf(stream, "\t\t [%td]", i);
                    print_number(stream, mesh->neighbors.rank, i);
                    print_number_pair(stream, mesh->neighbors.recv_off, i);
                    print_graph(stream, mesh->neighbors.send, i);
                    fprintf(stream, "\n");
                }
                fprintf(stream, "\t ]\n");
            }

            fflush(stream);
        }
        MPI_Barrier(sync.comm);
    }
}
