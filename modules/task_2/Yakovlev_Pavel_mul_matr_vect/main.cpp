// copyright : (C) by Pahandrovich
#include <mpi.h>
#include <assert.h>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>

int main(int argc, char* argv[]) {
  std::srand(static_cast<int>(time(NULL)));

  int rank, CountP, flag;
  if (argc < 3) {
    std::cout << "Error";
    return -1;
  }

  int rows = atoi(argv[1]);
  int cols = atoi(argv[2]);
  if (rows <= 0 || cols <= 0) {
    std::cout << "Error";
    return -1;
  }
  double *vect = NULL,
    *matrix = NULL;
  MPI_Status status;
  double Time_begin = 0;
  double Time_end = 0;
  int Matrix_size = rows*cols;
  int Vect_size = cols;
  int Result_size = rows;
  double *Result = NULL;
  double *tempRes = NULL;

  MPI_Init(&argc, &argv);
  MPI_Initialized(&flag);
  if (!flag) {
    std::cout << "Error";
    return -1;
  }

  MPI_Comm_size(MPI_COMM_WORLD, &CountP);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int submit_num = static_cast<int>(
                ceil(static_cast<double>(cols) / static_cast<double>(CountP)));
  int submit_num_last = cols - (CountP-1) * submit_num;

  if (rank == 0) {
    vect = new double[Vect_size];
    matrix = new double[Matrix_size];
    Result = new double[Result_size];
    tempRes = new double[Result_size];

//  -----------------Initialization-----------------
    for (int i = 0; i < Matrix_size; i++) {
      if (i < Vect_size) vect[i] = std::rand() % 10 - 5;
      if (i < Result_size) Result[i] = 0;
      matrix[i] = std::rand() % 10 - 5;
    }
//  -----------Print-----------------
    if (rows < 11 && cols < 11) {
      std::cout << "Vector:" << std::endl;
      for (int i = 0; i < Vect_size; i++)
        std::cout << " " << vect[i];
      std::cout << "\nMatrix:" << std::endl;
      for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++)
          std::cout << " " << matrix[i + rows * j];
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }

    std::cout << "Matr size = " << Matrix_size << std::endl;
    std::cout << "Vector size = " << Vect_size << std::endl;
    std::cout << "res size = " << Result_size << std::endl;
    std::cout << "rows = " << rows << std::endl;
    std::cout << "cols = " << cols << std::endl;
    std::cout << "rows * cols = " << rows * cols << std::endl << std::endl;

//  --------------  Line  --------------------
    Time_begin = MPI_Wtime();
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++)
        Result[i] += matrix[i + rows * j] * vect[j];
    }
    Time_end = MPI_Wtime();

    std::cout << "Line_Time = " << Time_end - Time_begin << std::endl;
    if (rows < 11 && cols < 11) {
      std::cout << "Line_Result: " << std::endl;
      for (int i = 0; i < Result_size; i++)
        std::cout << Result[i] << std::endl;
    }
    std::cout << std::endl;
    for (int i = 0; i < Result_size; i++) {
      Result[i] = 0.0;
      tempRes[i] = 0.0;
    }

// --------------  Parrallel  -------------
    if (CountP > 1) {
      Time_begin = MPI_Wtime();
      for (int i = 1; i < CountP-1; i++) {
        MPI_Send(vect + i * submit_num,
                  submit_num,
                  MPI_DOUBLE,
                  i,
                  0,
                  MPI_COMM_WORLD);
        MPI_Send(matrix + i * submit_num * rows,
                  submit_num * rows,
                  MPI_DOUBLE,
                  i,
                  0,
                  MPI_COMM_WORLD);
      }
      MPI_Send(vect + (CountP-1) * submit_num,
                  submit_num_last,
                  MPI_DOUBLE,
                  CountP-1,
                  0,
                  MPI_COMM_WORLD);
      MPI_Send(matrix + (CountP-1) * submit_num * rows,
                  submit_num_last * rows,
                  MPI_DOUBLE,
                  CountP-1,
                  0,
                  MPI_COMM_WORLD);

      for (int k = 0; k < submit_num; k++)
        for (int j = 0; j < rows; j++)
          Result[j] += vect[k] * matrix[j+rows*k];

      for (int i = 1; i < CountP; i++) {
        MPI_Recv(tempRes,
                  Result_size,
                  MPI_DOUBLE,
                  i,
                  0,
                  MPI_COMM_WORLD,
                  &status);
        for (int j = 0; j < Result_size; j++)
          Result[j] += tempRes[j];
      }

      Time_end = MPI_Wtime();

      std::cout << "parallel_time = " << Time_end - Time_begin << std::endl;
      if (rows < 11 && cols < 11) {
        std::cout << "Parallel_Result: " << std::endl;
        for (int i = 0; i < Result_size; i++)
          std::cout << Result[i] << std::endl;
      }
      std::cout << std::endl;
    }
  } else {
    if (rank != CountP-1) {
      vect = new double[submit_num];
      matrix = new double[submit_num*rows];
      tempRes = new double[Result_size];
      for (int i = 0; i < Result_size; i++)
        tempRes[i] = 0.0;
      MPI_Recv(vect,
                submit_num,
                MPI_DOUBLE,
                0,
                0,
                MPI_COMM_WORLD,
                &status);
      MPI_Recv(matrix,
                submit_num*rows,
                MPI_DOUBLE,
                0,
                0,
                MPI_COMM_WORLD,
                &status);

      for (int k = 0; k < submit_num; k++)
        for (int j = 0; j < Result_size; j++)
          tempRes[j] += vect[k] * matrix[j + rows * k];

      MPI_Send(tempRes, Result_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    } else {
      vect = new double[submit_num_last];
      matrix = new double[submit_num_last*rows];
      tempRes = new double[Result_size];
      for (int i = 0; i < Result_size; i++)
        tempRes[i] = 0.0;

      MPI_Recv(vect,
                submit_num_last,
                MPI_DOUBLE,
                0,
                0,
                MPI_COMM_WORLD,
                &status);
      MPI_Recv(matrix,
                submit_num_last *rows,
                MPI_DOUBLE,
                0,
                0,
                MPI_COMM_WORLD,
                &status);

      for (int k = 0; k < submit_num_last; k++)
        for (int j = 0; j < Result_size; j++)
          tempRes[j] += vect[k] * matrix[j + rows * k];

      MPI_Send(tempRes, Result_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
  }

  MPI_Finalize();

  delete[]vect;
  delete[]matrix;
  delete[]tempRes;
  delete[]Result;

  return 0;
}
