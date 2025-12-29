#include "sync.h"

long sync_all_lmin(long val)
{
    long min = val;
    MPI_Allreduce(&val, &min, 1, MPI_LONG, MPI_MIN, sync.comm);
    return min;
}

long sync_all_lmax(long val)
{
    long max = val;
    MPI_Allreduce(&val, &max, 1, MPI_LONG, MPI_MAX, sync.comm);
    return max;
}

long sync_all_lsum(long val)
{
    long sum = val;
    MPI_Allreduce(&val, &sum, 1, MPI_LONG, MPI_SUM, sync.comm);
    return sum;
}

double sync_all_fmin(double val)
{
    double min = val;
    MPI_Allreduce(&val, &min, 1, MPI_DOUBLE, MPI_MIN, sync.comm);
    return min;
}

double sync_all_fmax(double val)
{
    double max = val;
    MPI_Allreduce(&val, &max, 1, MPI_DOUBLE, MPI_MAX, sync.comm);
    return max;
}

double sync_all_fsum(double val)
{
    double sum = val;
    MPI_Allreduce(&val, &sum, 1, MPI_DOUBLE, MPI_SUM, sync.comm);
    return sum;
}
