#include <mpi.h>
#include <iostream>
#include <assert.h>
#define arrSize 100
#define MainProc 0
//Q?
//Запуск с командной строки
//Запуск из студии
//Отладка
// -np при запуску или -n
// когда собирается exe?
int main(int argc, char* argv[])
{
	int status, rank, size;
	int* arr = new int[arrSize];
	int* buff;
	int nProc;

	status = MPI_Init(&argc, &argv);
	assert(status == MPI_SUCCESS);

	status = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	assert(status == MPI_SUCCESS);
	
	status = MPI_Comm_size(MPI_COMM_WORLD, &size);
	assert(status == MPI_SUCCESS);

	/*std::cout << "Process #" << rank << '\n';
	std::cout << "Count process: " << size << '\n';*/
	int partialSum = 0;
	int arrSum = 0;
	nProc = arrSize / size;
	buff = new int[nProc];
	//Отправка
	if (rank == MainProc)
	{	
		for (int i = 0; i < arrSize; i++)
		{
				arr[i] = 1;
		}

		for (int i = 1; i < size; i++)
		{
			MPI_Send(&arr[i*nProc], nProc, MPI_INT, i, 0, MPI_COMM_WORLD);

		}

	}//принятие
	else 
		{   
			MPI_Recv(buff, nProc, MPI_INT, MainProc, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			std::cout << "Process #" << rank << '\n';
			std::cout << "nProc = " << nProc << std::endl;
			for (int i = 0; i < nProc; i++)
				partialSum += buff[i];
		}


	if (rank == MainProc)
	{	
		for (int i = 0; i < nProc; i++)
			partialSum += arr[i];
		arrSum += partialSum;
		for (int i = 1; i < size; i++)
		{
			MPI_Recv(&partialSum, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			arrSum += partialSum;
		}
		std::cout << arrSum << std::endl;
	}
	else
	{
		MPI_Send(&partialSum, 1, MPI_INT, MainProc, 0, MPI_COMM_WORLD);
	}


	//Освобождение памяти
	if (rank == 0)
	{ 
		delete[] arr;

	}
	status = MPI_Finalize();
	assert(status == MPI_SUCCESS);
	return 0;
}