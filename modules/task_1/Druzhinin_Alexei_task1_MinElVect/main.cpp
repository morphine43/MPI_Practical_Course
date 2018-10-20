#include <iostream>
#include <mpi.h>
#include <string>
#include <ctime>
#include <cstdlib>
#include <assert.h>

using namespace std;

int main(int argc, char *argv[])
{
	srand((int)time(0));
	
	double* Vector;
	double* myVector;
	int size;
	
	if (argc>1) {
		string sizeArg;
		sizeArg = argv[1];
		size = atoi(sizeArg.c_str());
	}
	else 
		size = 100;

	Vector = new double[size];
	
	// fill the vector
	for (int i = 0; i < size; i++) {
			Vector[i] = (double)(1+rand()%750) / ((double)(rand()%100) + 1);
 			//cout << "Vector[" << i << "]=" << Vector[i] << endl;
	}
	
	// variables for parallel block
	double Min = 0;  
	double startTime = 0;
	double endTime = 0;
	double Time = 0;
	double myMin = 0;
	int remainder = 0;
	int flag;
	int Id, numProcs;

	// variables for line block
	double Min_l = 0; 
	double startTime_l = 0;
	double endTime_l = 0;
	double Time_l = 0;
	
	// initialize MPI
	MPI_Init(&argc, &argv); 
	MPI_Initialized(&flag); // check
	if (!flag) {
		std::cout << "Error in MPI_Init!!!";
	}

	// communicators
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &Id);
	
	if (Id == 0) { 
	//LINE BLOCK
		startTime_l = MPI_Wtime();
		Min_l = Vector[0];
	    for (int i = 1; i < size; i++) {
			if (Min_l > Vector[i])
				Min_l = Vector[i];
		}
	//Line Results
	    endTime_l = MPI_Wtime();
		cout << endl << "===LINE===" << endl;
		cout << "Min: " << Min_l << endl;
		Time_l = endTime_l - startTime_l;
		cout << "Time: " << Time_l << " sec" << endl;
	//END LINE BLOCK

	//PARALLEL BLOCK
	startTime = MPI_Wtime(); // started counting time in 0 proccess
	}
	
	
	myVector = new double[size/numProcs];
	for (int i = size / numProcs * Id, j = 0; i < size / numProcs * Id + size / numProcs; i++, j++)
		myVector[j] = Vector[i];
	
	MPI_Bcast(myVector, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	myMin = myVector[0];
	for (int i = 0; i < size / numProcs; i++) {
		if (myMin > myVector[i])
			myMin = myVector[i];
	}

	MPI_Reduce(&myMin, &Min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    
	if (Id == 0) {
	     remainder = size % numProcs;
		 if (remainder) {
			for (int i = 0; i < remainder; i++)
				if (Min > Vector[size - remainder  + i])
					Min = Vector[size - remainder  + i];
		 }

	//Parallel Results 
		endTime = MPI_Wtime();
		cout << endl <<"===PARALLEL===" << endl;
		cout << "Min: " << Min << endl;
		Time = endTime - startTime;
		cout << "Time: " << Time << " sec" << endl << endl;
	}
	MPI_Finalize();
	//END PARALLEL BLOCK
	
	delete[]Vector;

	return 0;
}
