#include "sync.h"

int sync_min(int val)
{
    int min = val;
    MPI_Allreduce(&val, &min, 1, MPI_INT, MPI_MIN, sync.comm);
    return min;
}

int sync_max(int val)
{
    int max = val;
    MPI_Allreduce(&val, &max, 1, MPI_INT, MPI_MAX, sync.comm);
    return max;
}

long sync_lsum(long val)
{
    long sum = val;
    MPI_Allreduce(&val, &sum, 1, MPI_LONG, MPI_SUM, sync.comm);
    return sum;
}

long sync_lexsum(long val)
{
    long exsum = 0;
    MPI_Exscan(&val, &exsum, 1, MPI_LONG, MPI_SUM, sync.comm);
    return (sync.rank == 0) ? 0 : exsum;
}

double sync_fmin(double val)
{
    double min = val;
    MPI_Allreduce(&val, &min, 1, MPI_DOUBLE, MPI_MIN, sync.comm);
    return min;
}

double sync_fmax(double val)
{
    double max = val;
    MPI_Allreduce(&val, &max, 1, MPI_DOUBLE, MPI_MAX, sync.comm);
    return max;
}

double sync_fsum(double val)
{
    double sum = val;
    MPI_Allreduce(&val, &sum, 1, MPI_DOUBLE, MPI_SUM, sync.comm);
    return sum;
}
