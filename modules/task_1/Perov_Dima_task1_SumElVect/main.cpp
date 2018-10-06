#include <iostream>
#include <mpi.h>

const int SIZE = 10000;

using namespace std;

int main(int argc, char *argv[])
{
	double myVector[SIZE];
	int myId, numProcs;
	double startTime = 0;
	double endTime;

	for (int i = 0; i < SIZE; i++) {
		myVector[i] = 1;
	}

	int flag; // Флаг инциализации
	MPI_Status status;
	double mySum = 0;
	double myResult = 0;

	MPI_Init(&argc, &argv); // Инциализируем MPI-окружение
	MPI_Initialized(&flag); // Проверка
	if (!flag) {
		std::cout << "Error: MPI_Init";
	}

	// Описание коммуникатора
	// Коммуникатор управляет группами параллельных процеcсов
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs); // Определение количество процессов в группе
	MPI_Comm_rank(MPI_COMM_WORLD, &myId); // Определение ранга процесса в группе

	
	if (myId == 0) { 
	startTime = MPI_Wtime();
	}

	MPI_Bcast(&myVector[SIZE / numProcs * myId], SIZE / numProcs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (int i = SIZE / numProcs * myId; i < SIZE / numProcs + SIZE / numProcs * myId; i++) {
		mySum += myVector[i];
	}

	MPI_Reduce(&mySum, &myResult, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


	if (myId == 0) {
		cout << "Result: " << myResult << endl;
		endTime = MPI_Wtime();
		cout << "Time: " << (endTime - startTime) << " sec" << endl;
		//system("pause");
	}


	MPI_Finalize();

}
