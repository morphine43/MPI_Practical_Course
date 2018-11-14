#include <mpi.h>
#include <assert.h>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>


using namespace std;

int main (int argc, char* argv[])
{
    srand((int)time(NULL));

    int rank, CountP, flag;
    if (argc < 3) 
    {
        cout << "Error";
        return -1;
    }

    int rows = atoi(argv[1]);
	int cols = atoi(argv[2]);
    if (rows <= 0 || cols <= 0)
    {
        cout << "Error";
        return -1;
    }

    double *vect = NULL,
			*matrix = NULL;
    MPI_Status status; 
    double Line_Mult = 0;  
    double Time_begin = 0;
    double Time_end = 0;
	int Matrix_size = rows*cols;
	int Vect_size = cols;
	int Result_size = rows;
	double *Result = NULL;
	double *tempRes = NULL;
	

    MPI_Init (&argc, &argv);
    MPI_Initialized(&flag);
	if (!flag) 
    {
		cout << "Error";
        return -1;
	}


    MPI_Comm_size(MPI_COMM_WORLD, &CountP);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int submit_num = (int)ceil((double)cols / (double)CountP); // кол-во отправляемых столбцов
	int submit_num_last = cols - (CountP-1) * submit_num; // кол-во отправляемых столбцов для последнего процесса
	
	if (rank == 0)
	{
		vect = new double[Vect_size];
		matrix = new double[Matrix_size];
		Result = new double[Result_size];
		tempRes = new double[Result_size];

		// -----------------Initialization-----------------
		for (int i = 0; i < Matrix_size; i++)
		{
			if (i < Vect_size) vect[i] = rand() % 10 - 5;
			if (i < Result_size) Result[i] = 0;
			matrix[i] = rand() % 10 - 5;
		}
		// -----------Print-----------------
		if (rows < 11 && cols < 11)
		{
			cout << "Vector:" << endl;
			for (int i = 0; i < Vect_size; i++)
			{
				cout << " " << vect[i];
			}
			cout << "\nMatrix:" << endl;
			for (int i = 0; i < rows; i++)
			{
				for (int j = 0; j < cols; j++)
					cout << " " << matrix[i + rows * j];
				cout << endl;
			}
			cout << endl;
		}

		cout << "Matr size = " << Matrix_size << endl;
		cout << "Vector size = " << Vect_size << endl;
		cout << "res size = " << Result_size << endl;
		cout << "rows = " << rows << endl;
		cout << "cols = " << cols << endl;
		cout << "rows * cols = " << rows * cols << endl << endl;

		//  --------------  Line  --------------------
		Time_begin = MPI_Wtime();
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++)
			{
				Result[i] += matrix[i + rows * j] * vect[j];
			}
		}
		Time_end = MPI_Wtime();

		cout << "Line_Time = " << Time_end - Time_begin << endl;
		if (rows < 11 && cols < 11)
		{
			cout << "Line_Result: " << endl;
			for (int i = 0; i < Result_size; i++)
				cout << Result[i] << endl;
		}
		cout << endl;
		for (int i = 0; i < Result_size; i++)
		{		
			Result[i] = 0.0;
			tempRes[i] = 0.0;
		}
        
        // --------------  Parrallel  -------------
         if (CountP > 1)
		 {
			 
             Time_begin = MPI_Wtime();
             for (int i = 1; i < CountP-1; i++)
             {
                 MPI_Send(vect + i * submit_num, submit_num, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                 MPI_Send(matrix + i * submit_num * rows, submit_num * rows, MPI_DOUBLE,i,0,MPI_COMM_WORLD);
             }
			 MPI_Send(vect + (CountP-1) * submit_num, submit_num_last , MPI_DOUBLE, CountP-1, 0, MPI_COMM_WORLD);
			 MPI_Send(matrix + (CountP-1) * submit_num * rows, submit_num_last * rows, MPI_DOUBLE, CountP-1, 0, MPI_COMM_WORLD);
 
			 for (int k = 0; k < submit_num;k++)
				for (int j = 0; j < rows; j++)
					Result[j] += vect[k] * matrix[j+rows*k];
			  


             for (int i = 1; i < CountP; i++)
             {
                 MPI_Recv(tempRes,Result_size,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
				 for (int j = 0; j < Result_size; j++)
					 Result[j] += tempRes[j];
             }

             Time_end = MPI_Wtime();
 
			 cout << "parallel_time = " << Time_end - Time_begin << endl;			
			 if (rows < 11 && cols < 11)
			 {
				 cout << "Parallel_Result: " << endl;
				 for (int i = 0; i < Result_size; i++)
					 cout << Result[i] << endl;
			 }
			 cout << endl;	 
        }
    }
    else
        if (rank != CountP-1)
        {
            vect = new double[submit_num];
            matrix = new double[submit_num*rows];
			tempRes = new double[Result_size];
			for (int i = 0; i < Result_size; i++)
				tempRes[i] = 0.0;
            MPI_Recv(vect,submit_num,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);  
            MPI_Recv(matrix,submit_num*rows,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);

			for (int k = 0; k < submit_num; k++)
				for (int j = 0; j < Result_size; j++)
					tempRes[j] += vect[k] * matrix[j + rows * k];

            MPI_Send(tempRes, Result_size, MPI_DOUBLE,0,0,MPI_COMM_WORLD);
        } 
        else
        {
			vect = new double[submit_num_last];
			matrix = new double[submit_num_last*rows];
			tempRes = new double[Result_size];
			for (int i = 0; i < Result_size; i++)
				tempRes[i] = 0.0;

			MPI_Recv(vect, submit_num_last, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
			MPI_Recv(matrix, submit_num_last *rows, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

			for (int k = 0; k < submit_num_last; k++)
				for (int j = 0; j < Result_size; j++)
					tempRes[j] += vect[k] * matrix[j + rows * k];

            MPI_Send(tempRes, Result_size, MPI_DOUBLE,0,0,MPI_COMM_WORLD);
        } 

    MPI_Finalize();

    delete[]vect;
    delete[]matrix;
	delete[]tempRes;
	delete[]Result;
    
    return 0;
}
