#include <mpi.h>
#include <iostream>
#include <assert.h>

#define MainProc 0
int main(int argc, char* argv[])
{   
	// mpi variables
	int status, rank, size;
	double t1, t2;
	//usual variables
	int  column,
		rows,
		arrSize,
		partialSum = 0,
		arrSum = 0,
		nProc = 0,
		N; // number of elem that will be given to one process
	int *arr, //matrix
		*buff; // buffer for message exchanging
	//init
	rows = (atoi(argv[1]) > 0) ? atoi(argv[1]) : 1;
	column = (atoi(argv[2]) > 0) ? atoi(argv[2]) : 1;
	arrSize = rows * column;

	//mpi part
	status = MPI_Init(&argc, &argv);
	assert(status == MPI_SUCCESS);

	status = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	assert(status == MPI_SUCCESS);

	status = MPI_Comm_size(MPI_COMM_WORLD, &size);
	assert(status == MPI_SUCCESS);

	if (size > 64) return -1; // limit 
	if (rank == MainProc) t1 = MPI_Wtime();
	
	nProc = (arrSize >= size) ? size : arrSize;
	N = arrSize / nProc;
	
	if (rank == MainProc) // sending
	{
		arr = new int[arrSize];

		for (int i = 0; i < arrSize; i++)
			arr[i] = 1;

		buff = arr;

		for (int i = 1; i < nProc - 1; i++)
			status = MPI_Send(&arr[i*N], N, MPI_INT, i, i, MPI_COMM_WORLD);

		if (size != 1) // sending last part separately in case if arrSize % size != 0
		    MPI_Send(&arr[N*(nProc - 1)], arrSize - N * (nProc - 1), MPI_INT, nProc - 1, nProc - 1, MPI_COMM_WORLD);
	

	}
	else //receiving
	{
		if (rank == nProc - 1)
		{
			N = arrSize - N * (nProc - 1);
			buff = new int[N];
		    MPI_Recv(buff, N, MPI_INT, MainProc, rank, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		}
		else
		{
			buff = new int[N];
			MPI_Recv(buff, N, MPI_INT, MainProc, rank, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		}
	}

	//common part 
	for (int i = 0; i < N; i++)
		partialSum += buff[i];

	if (rank == MainProc)
	{
		arrSum += partialSum;
		
		for (int i = 1; i < size; i++)
		{
			MPI_Recv(&partialSum, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			arrSum += partialSum;
		}

		std::cout << arrSum << std::endl;
	}
	else
	{
		MPI_Send(&partialSum, 1, MPI_INT, MainProc, 0, MPI_COMM_WORLD);
	}



	if (rank == 0)
		delete[] arr;
	else
		delete[] buff;
		
	if (rank == MainProc)
	{
		t2 = MPI_Wtime();
		std::cout << "Wtime: " << t2 - t1 << std::endl;
	}
	status = MPI_Finalize();
	assert(status == MPI_SUCCESS);
	return 0;
}
