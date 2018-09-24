#include <mpi.h>
#include <iostream>
#include <assert.h>
#define arrSize 2
//Q?
//Запуск с командной строки
//Запуск из студии
//Отладка
int main(int argc, char* argv[])
{
	int status, rank, size;
	int** arr = new int*[arrSize];
	
		for (int i = 0; i < arrSize; i++)
		{
			arr[i] = new int[arrSize];
			// init
			for (int j = 0; j < arrSize; j++)
			{
				arr[i][j] = 1;
			//	std::cout << arr[i][j] << " ";
			}
		//	std::cout << std::endl;
		}

	status = MPI_Init(&argc, &argv);
	assert(status == MPI_SUCCESS);

	status = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	assert(status == MPI_SUCCESS);
	
	status = MPI_Comm_size(MPI_COMM_WORLD, &size);
	assert(status == MPI_SUCCESS);

	/*std::cout << "Process #" << rank << '\n';
	std::cout << "Count process: " << size << '\n';*/
	int sum[2] = { 0,0 };
	int arrSum = 0;
	if (rank == 0)
	{
		/*for (int i = 0; i < arrSize; i++)
			sum += arr[0][i];*/
		MPI_Send(arr[0],arrSize, MPI_INT, 1, 0, MPI_COMM_WORLD);
		MPI_Send(arr[1], arrSize, MPI_INT, 2, 1, MPI_COMM_WORLD);
		MPI_Recv(&sum[0], 1, MPI_INT, 1, 2, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		MPI_Recv(&sum[1], 1, MPI_INT, 2, 3, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
	}
	else 
		if (rank == 1)
		{
			MPI_Recv(arr[0], arrSize, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			for (int i = 0; i < arrSize; i++)
				sum[0] += arr[0][i];
			MPI_Send(&sum[0], 1, MPI_INT, 0, 2, MPI_COMM_WORLD);

		}
		else 
			if (rank == 2)
			{
				MPI_Recv(arr[0], arrSize, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
				for (int i = 0; i < arrSize; i++)
					sum[1] += arr[0][i];
				MPI_Send(&sum[1], 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
			}
		
	status = MPI_Finalize();
	assert(status == MPI_SUCCESS);
	arrSum = sum[0] + sum[1];
	std::cout << arrSum << std::endl;
	return 0;
}