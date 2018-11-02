#include <mpi.h>
#include <assert.h>
#include <iostream>
#include <string>
#include <ctime>
#include <cstdlib>


using namespace std;

int main (int argc, char* argv[])
{
    srand((int)time(NULL));

    int rank, CountP, flag;
    if (argc < 2) 
    {
        cout << "Error";
        return -1;
    }

    int sizev = atoi(argv[1]);
    if (sizev <= 0)
    {
        cout << "Error";
        return -1;
    }

    double *v1,*v2;
    MPI_Status status; 
    double Line_MultVec = 0;  
    double MultVec = 0;
    double tempSum = 0;
    double Time_begin = 0;
    double Time_end = 0;

    MPI_Init (&argc, &argv);
    MPI_Initialized(&flag);
	if (!flag) 
    {
		cout << "Error";
        return -1;
	}
    MPI_Comm_size(MPI_COMM_WORLD, &CountP);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (rank == 0)
    {
        
        v1 = new double[sizev];
        v2 = new double[sizev];
        for (int i = 0;i<sizev;i++)
        {
            v1[i] = rand()%30-15;
            v2[i] = rand()%30-15;
        }

        if (sizev < 21)
        {
            cout << "v1:";
            for (int i = 0; i<sizev; i++)
            {
                cout <<" "<<v1[i];
            }
            cout << "\nv2:";
            for (int i = 0; i<sizev; i++)
            {
                cout <<" "<<v2[i];
            }
            cout << endl;
        }

        //  --------------  Line  --------------------
        Time_begin = MPI_Wtime();
        for (int i = 0; i<sizev;i++)
            Line_MultVec += v1[i]*v2[i];
        Time_end = MPI_Wtime();
        cout << "Line_Result = " << Line_MultVec << endl;
        cout << "Line_Time = " << Time_end - Time_begin << endl;
        
        // --------------  Parrallel  -------------
        if (CountP > 1){
            Time_begin = MPI_Wtime();
            for (int i = 1; i < CountP-1; i++)
            {
                MPI_Send(v1+i*(sizev/(CountP)), sizev/(CountP), MPI_DOUBLE,i,0,MPI_COMM_WORLD);
                MPI_Send(v2+i*(sizev/(CountP)), sizev/(CountP), MPI_DOUBLE,i,0,MPI_COMM_WORLD);
            }
            MPI_Send((v1+(CountP-1)*(sizev/(CountP))), (sizev - (CountP-1)*(sizev/(CountP))), MPI_DOUBLE, (CountP-1), 0, MPI_COMM_WORLD);
            MPI_Send((v2+(CountP-1)*(sizev/(CountP))), (sizev - (CountP-1)*(sizev/(CountP))), MPI_DOUBLE, (CountP-1), 0, MPI_COMM_WORLD);

            for (int i = 0; i < (sizev/CountP); i++)
                MultVec += v1[i]*v2[i];

            for (int i = 1; i < CountP; i++)
            {
                MPI_Recv(&tempSum,1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
                MultVec += tempSum;  
            }
            Time_end = MPI_Wtime();

            cout << "\nParallel_Result = " << MultVec << endl;
            cout << "Parallel_Time = " << Time_end - Time_begin << endl;
        }
    }
    else
        if (rank != CountP-1)
        {
            v1 = new double[sizev/(CountP)];
            v2 = new double[sizev/(CountP)];
            MPI_Recv(v1,sizev/(CountP),MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);  
            MPI_Recv(v2,sizev/(CountP),MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);
            for (int i = 0; i < sizev/(CountP); i++)
                tempSum += v1[i]*v2[i]; 
            MPI_Send(&tempSum, 1, MPI_DOUBLE,0,0,MPI_COMM_WORLD);
        } 
        else
        {
            v1 = new double[(sizev - (CountP-1)*(sizev/(CountP)))];
            v2 = new double[(sizev - (CountP-1)*(sizev/(CountP)))];
            MPI_Recv(v1,(sizev - (CountP-1)*(sizev/(CountP))),MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);  
            MPI_Recv(v2,(sizev - (CountP-1)*(sizev/(CountP))),MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);
            for (int i = 0; i < (sizev - (CountP-1)*(sizev/(CountP))); i++)
                tempSum += v1[i]*v2[i]; 
            MPI_Send(&tempSum, 1, MPI_DOUBLE,0,0,MPI_COMM_WORLD);
        }
    MPI_Finalize();

    delete[]v1;
    delete[]v2;
    
    return 0;
}
