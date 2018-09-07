#include <mpi.h>
#include <iostream>
#include <assert.h>

int main (int argc, char* argv[])
{
    int status, rank, size;
    status = MPI_Init (&argc, &argv);
    assert(status == MPI_SUCCESS);

    status = MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    assert(status == MPI_SUCCESS);

    std::cout << "Hello process #" << rank << '\n';

    status = MPI_Finalize();
    assert(status == MPI_SUCCESS);

    return 0;
}
