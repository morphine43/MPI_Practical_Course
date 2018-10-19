#include <math.h> 
#include <iostream>
#include <mpi.h>
#include <string>
#include <ctime>
#include <cstdlib>
#include <assert.h>

int main(int argc, char *argv[])
{ 
	
	srand((int)time(0));
	
	double* x;
	int N;
	
	if (argc>1) {
		N= atoi(argv[1]);
	}
	else 
		N = 50;

	x = new double[N]; 
	
	for(int i=0;i<N;i++) x[i] = (double)(1+rand()%7) / ((double)(rand()%10) + 1); 
	
	
	int ProcRank, ProcNum; 	
		
	MPI_Status Status; 

	MPI_Init(&argc,&argv); 
	MPI_Comm_size(MPI_COMM_WORLD,&ProcNum); 
	MPI_Comm_rank(MPI_COMM_WORLD,&ProcRank); 
	int ParallelSum = 0; 
	double ParallelTimeStart=0;
	int tempSum = 0; 
	if ( ProcRank == 0 ) 
	{ 
		//serial v 

		int SerialSum = 0; 
		double SerialTimeStart = MPI_Wtime(); 

		for(int i=0;i<N;i++) 
		{ 
			SerialSum+=x[i]; 
		} 
		double SerialMidscore = (double)SerialSum/N; 
		double SerialTime = MPI_Wtime() - SerialTimeStart; 

		printf("\nMidscore (serial) = %3.3f",SerialMidscore); 
		printf("\nTime (serial) = %3.6f",SerialTime); 

		//parallel v 

		ParallelTimeStart = MPI_Wtime(); 
		int p = N/ProcNum; 
		for(int i=1;i<ProcNum-1;i++) 
		{ 
			int jn = p*i; 
			int jk = p*(i+1); 
			if(i==ProcNum - 1) jk = N; 
			int blockSize = jk - jn; 
			MPI_Send(&blockSize, 1, MPI_INT, i, 0, MPI_COMM_WORLD); 
			MPI_Send(x+jn, blockSize, MPI_INT, i, 0, MPI_COMM_WORLD); 
		} 
		int blockSize = N - p*(ProcNum - 1); 
		MPI_Send(&blockSize, 1, MPI_INT, ProcNum-1, 0, MPI_COMM_WORLD);
		MPI_Send(x+p*(ProcNum - 1), blockSize, MPI_INT, ProcNum-1, 0, MPI_COMM_WORLD); 
		for(int i=0; i<N/ProcNum;i++) 
		{ 
			ParallelSum+=x[i]; 
		} 
	} 
	else 
	{ 
		int blockSize; 
		int *temp; 
		MPI_Recv(&blockSize, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &Status); 
		temp = new int[blockSize]; 
		MPI_Recv(temp, blockSize, MPI_INT, 0, 0, MPI_COMM_WORLD, &Status); 
		for(int i=0;i<blockSize;i++) 
		{ 
			tempSum+=temp[i]; 
		} 
		delete []temp; 
	} 
	int sum = 0; 
	MPI_Reduce(&tempSum, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
	ParallelSum+=sum; 
	if(ProcRank==0) 
	{ 
		double ParallelMidscore = (double)ParallelSum/N; 
		double ParallelTime = MPI_Wtime() - ParallelTimeStart; 
		printf("\nMidscore (parallel) = %3.3f",ParallelMidscore); 
		printf("\nTime (parallel) = %3.6f",ParallelTime); 
	} 
	MPI_Finalize(); 
	delete []x; 
}
