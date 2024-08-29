#include "sync.h"

#include <math.h>

#include "array.h"
#include "memory.h"

struct Sync sync = {0};

void sync_init(int *argc, char ***argv)
{
    MPI_Init(argc, argv);
    sync.comm = MPI_COMM_WORLD;
    MPI_Comm_rank(sync.comm, &sync.rank);
    MPI_Comm_size(sync.comm, &sync.size);
}

long x__sync_min_long(long v)
{
    long min = 0;
    MPI_Allreduce(&v, &min, 1, MPI_LONG, MPI_MIN, sync.comm);
    return min;
}
double x__sync_min_double(double v)
{
    double min = 0;
    MPI_Allreduce(&v, &min, 1, MPI_DOUBLE, MPI_MIN, sync.comm);
    return min;
}

long x__sync_max_long(long v)
{
    long max = 0;
    MPI_Allreduce(&v, &max, 1, MPI_LONG, MPI_MAX, sync.comm);
    return max;
}
double x__sync_max_double(double v)
{
    double max = 0;
    MPI_Allreduce(&v, &max, 1, MPI_DOUBLE, MPI_MAX, sync.comm);
    return max;
}

long x__sync_sum_long(long v)
{
    long sum = 0;
    MPI_Allreduce(&v, &sum, 1, MPI_LONG, MPI_SUM, sync.comm);
    return sum;
}
double x__sync_sum_double(double v)
{
    double sum = 0;
    MPI_Allreduce(&v, &sum, 1, MPI_DOUBLE, MPI_SUM, sync.comm);
    return sum;
}

long x__sync_exsum_long(long v)
{
    long exsum = 0;
    MPI_Exscan(&v, &exsum, 1, MPI_LONG, MPI_SUM, sync.comm);
    return exsum;
}
double x__sync_exsum_double(double v)
{
    double exsum = 0;
    MPI_Exscan(&v, &exsum, 1, MPI_DOUBLE, MPI_SUM, sync.comm);
    return exsum;
}

double sync_dot(const double *a, const double *b, long n)
{
    return x__sync_sum_double(array_dot(a, b, n));
}

double sync_norm(const double *a, long n)
{
    return sqrt(sync_dot(a, a, n));
}

void x__sync_gatherv_long(long *a, const long *b, long n, int root)
{
    smart int *counts = (sync.rank == root ? memory_calloc(sync.size, sizeof(*counts)) : 0);
    smart int *displs = (sync.rank == root ? memory_calloc(sync.size, sizeof(*displs)) : 0);
    MPI_Gather(&(int[]){n}, 1, MPI_INT, counts, 1, MPI_INT, root, sync.comm);
    if (sync.rank == root)
        for (int rank = 1; rank < sync.size; ++rank)
            displs[rank] = displs[rank - 1] + counts[rank - 1];
    MPI_Gatherv(b, n, MPI_LONG, a, counts, displs, MPI_LONG, root, sync.comm);
}
void x__sync_gatherv_double(double *a, const double *b, long n, int root)
{
    smart int *counts = (sync.rank == root ? memory_calloc(sync.size, sizeof(*counts)) : 0);
    smart int *displs = (sync.rank == root ? memory_calloc(sync.size, sizeof(*displs)) : 0);
    MPI_Gather(&(int[]){n}, 1, MPI_INT, counts, 1, MPI_INT, root, sync.comm);
    if (sync.rank == root)
        for (int rank = 1; rank < sync.size; ++rank)
            displs[rank] = displs[rank - 1] + counts[rank - 1];
    MPI_Gatherv(b, n, MPI_DOUBLE, a, counts, displs, MPI_DOUBLE, root, sync.comm);
}

void x__sync_scatterv_long(long *a, const long *b, long n, int root)
{
    smart int *counts = (sync.rank == root ? memory_calloc(sync.size, sizeof(*counts)) : 0);
    smart int *displs = (sync.rank == root ? memory_calloc(sync.size, sizeof(*displs)) : 0);
    MPI_Gather(&(int[]){n}, 1, MPI_INT, counts, 1, MPI_INT, root, sync.comm);
    if (sync.rank == root)
        for (int rank = 1; rank < sync.size; ++rank)
            displs[rank] = displs[rank - 1] + counts[rank - 1];
    MPI_Scatterv(b, counts, displs, MPI_LONG, a, n, MPI_LONG, root, sync.comm);
}
void x__sync_scatterv_double(double *a, const double *b, long n, int root)
{
    smart int *counts = (sync.rank == root ? memory_calloc(sync.size, sizeof(*counts)) : 0);
    smart int *displs = (sync.rank == root ? memory_calloc(sync.size, sizeof(*displs)) : 0);
    MPI_Gather(&(int[]){n}, 1, MPI_INT, counts, 1, MPI_INT, root, sync.comm);
    if (sync.rank == root)
        for (int rank = 1; rank < sync.size; ++rank)
            displs[rank] = displs[rank - 1] + counts[rank - 1];
    MPI_Scatterv(b, counts, displs, MPI_DOUBLE, a, n, MPI_DOUBLE, root, sync.comm);
}

void sync_irecv(const long *j_recv, MPI_Request *req, double *u, long ldu)
{
    for (int rank = 0; rank < sync.size; ++rank) {
        const int count = (j_recv[rank + 1] - j_recv[rank]) * ldu;
        MPI_Irecv(&u[j_recv[rank] * ldu], count, MPI_DOUBLE, rank, 0, sync.comm, &req[rank]);
    }
}

void sync_isend(const long *i_send, const long *send, const double *u, MPI_Request *req,
                double *buf, long ldu)
{
    for (int rank = 0; rank < sync.size; ++rank)
        for (long i = i_send[rank]; i < i_send[rank + 1]; ++i)
            for (long j = 0; j < ldu; ++j) buf[i * ldu + j] = u[send[i] * ldu + j];

    for (int rank = 0; rank < sync.size; ++rank) {
        const int count = (i_send[rank + 1] - i_send[rank]) * ldu;
        MPI_Isend(&buf[i_send[rank] * ldu], count, MPI_DOUBLE, rank, 0, sync.comm, &req[rank]);
    }
}

void sync_waitall(MPI_Request *req)
{
    MPI_Waitall(sync.size, req, MPI_STATUSES_IGNORE);
}

void sync_finalize(void)
{
    MPI_Finalize();
}
