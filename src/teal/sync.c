#include "sync.h"

Sync sync = {0};

MPI_Datatype vector_type;

static int tag_ub;

void sync_init(int *argc, char ***argv)
{
    MPI_Init(argc, argv);

    MPI_Comm_dup(MPI_COMM_WORLD, &sync.comm);
    MPI_Comm_rank(sync.comm, &sync.rank);
    MPI_Comm_size(sync.comm, &sync.size);

    MPI_Type_contiguous(3, MPI_DOUBLE, &vector_type);
    MPI_Type_commit(&vector_type);

    enum { MPI_TAG_UB_MIN = 32767 };
    int *attr = 0;
    int flag = 0;
    MPI_Comm_get_attr(sync.comm, MPI_TAG_UB, &attr, &flag);
    tag_ub = (flag && attr) ? *attr : MPI_TAG_UB_MIN;
}

void sync_reinit(MPI_Comm comm)
{
    MPI_Comm_free(&sync.comm);
    sync.comm = comm;
    MPI_Comm_rank(sync.comm, &sync.rank);
    MPI_Comm_size(sync.comm, &sync.size);
}

int sync_tag(void)
{
    static int tag = 0;
    return tag = (tag % tag_ub) + 1;
}

long sync_lmin(long val)
{
    long min = val;
    MPI_Allreduce(&val, &min, 1, MPI_LONG, MPI_MIN, sync.comm);
    return min;
}

long sync_lmax(long val)
{
    long max = val;
    MPI_Allreduce(&val, &max, 1, MPI_LONG, MPI_MAX, sync.comm);
    return max;
}

long sync_lsum(long val)
{
    long sum = 0;
    MPI_Allreduce(&val, &sum, 1, MPI_LONG, MPI_SUM, sync.comm);
    return sum;
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
    double sum = 0;
    MPI_Allreduce(&val, &sum, 1, MPI_DOUBLE, MPI_SUM, sync.comm);
    return sum;
}

vector sync_vmin(vector val)
{
    vector min = val;
    MPI_Allreduce(&val, &min, 3, MPI_DOUBLE, MPI_MIN, sync.comm);
    return min;
}

vector sync_vmax(vector val)
{
    vector max = val;
    MPI_Allreduce(&val, &max, 3, MPI_DOUBLE, MPI_MAX, sync.comm);
    return max;
}

vector sync_vsum(vector val)
{
    vector sum = {0};
    MPI_Allreduce(&val, &sum, 3, MPI_DOUBLE, MPI_SUM, sync.comm);
    return sum;
}

long sync_exsum(long val)
{
    long exsum = 0;
    MPI_Exscan(&val, &exsum, 1, MPI_LONG, MPI_SUM, sync.comm);
    return exsum;
}

void sync_finalize(void)
{
    MPI_Type_free(&vector_type);
    MPI_Comm_free(&sync.comm);
    MPI_Finalize();
}
