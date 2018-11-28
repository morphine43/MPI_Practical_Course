/* Copyright Danilov Nikita */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <limits.h>
#include <memory.h>
#include <cstdlib>
#include <cstddef>
#include <ctime>
#include <iostream>
#include <string>

#define N 1000000
#define DEBUG 0
#define MASTER_PROCESS 0

double start_time, stop_time;

#if DEBUG
void print_info(int id, std::string msg) {
    std::cout << id << ": " << msg << (clock() - start_time) / CLOCKS_PER_SEC
              << " secs\n";
}
#endif

void print_vector(int *v, int vector_size) {
    int i;

    for (i = 0; i < vector_size; i++)
        printf("\tv[%2d] = %d\n", i, v[i]);
}

int *merge(int *v1, int n1, int *v2, int n2) {
    int *result;
    int i = 0,
        j = 0,
        k = 0;

    result = new int[n1 + n2];

    while (i < n1 && j < n2) {
        if (v1[i] < v2[j]) {
            result[k] = v1[i];
            i++; k++;
        } else {
            result[k] = v2[j];
            j++; k++;
        }
    }

    if (i == n1) {
        while (j < n2) {
            result[k] = v2[j];
            j++; k++;
        }
    } else {
        while (i < n1) {
            result[k] = v1[i];
            i++; k++;
        }
    }

    return result;
}

void my_swap(int *v, int i, int j) {
    int t;

    t = v[i];
    v[i] = v[j];
    v[j] = t;
}

void quick_sort(int *v, int left, int right) {
    int i, last;

    if (left >= right)
        return;

    my_swap(v, left, (left + right) / 2);
    last = left;
    for (i = left + 1; i <= right; i++)
        if (v[i] < v[left])
            my_swap(v, ++last, i);

    my_swap(v, left, last);
    quick_sort(v, left, last-1);
    quick_sort(v, last + 1, right);
}

int main(int argc, char **argv) {
    double seq_time = 0,
           par_time = 0;
    int m, vector_size = N;
    int id, proc_number;
    int *data   = nullptr,
        *datacp = nullptr;
    MPI_Status status;
    int chunk_size;
    int *chunk;
    int *other;
    int step;
    int i;

    if (argc > 1) {
        vector_size = atoi(argv[1]);
    }

    start_time = clock();

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_number);

#if DEBUG
    print_info(id, "MPI setup complete");
#endif

    if (id == MASTER_PROCESS) {
        int extra_elements;
        int data_size;

        std::srand(unsigned(std::time(NULL)));

        extra_elements = vector_size % proc_number;
        chunk_size = vector_size / proc_number;
        data_size = vector_size + ((extra_elements == 0) ? 0 :
                    proc_number - extra_elements);

        /* allocate memory for vector
         * increase the size of the array so that
         * it is divided without a balance by the number of processors
         */
        data   = new int[data_size];
        datacp = new int[vector_size];

        /* fill main elements of vector
         * rand() + rand() can give overflow int,
         * thus we get negative numbers in the vector.
         * it needs to be done because
         * rand() returns non-negative numbers
         */
        for (i = 0; i < vector_size; i++)
            data[i] = std::rand() + std::rand();

        /* copy data to another array for sequential sorting */
        memcpy(datacp, data, vector_size * sizeof(int));

        /* fill extra elements of vector
         * fill in the maximum number so that after sorting
         * you know which extra elements need to be thrown out
         */
        if (extra_elements != 0) {
            for (i = vector_size;
                 i < data_size;
                 i++)
                data[i] = INT_MAX;
            chunk_size++;
        }

        if (vector_size <= 20) {
            printf("Vector before sort: \n");
            print_vector(data, data_size);
        }

#if DEBUG
        print_info(id, "Generated the random numbers");
#endif

        MPI_Bcast(&chunk_size, 1, MPI_INT, MASTER_PROCESS, MPI_COMM_WORLD);
        chunk = new int[chunk_size];
        MPI_Scatter(data, chunk_size, MPI_INT, chunk, chunk_size, MPI_INT,
                    MASTER_PROCESS, MPI_COMM_WORLD);

#if DEBUG
        print_info(id, "Scattered data");
#endif

        quick_sort(chunk, 0, chunk_size - 1);

#if DEBUG
        print_info(id, "Sorted");
#endif

    } else {
        MPI_Bcast(&chunk_size, 1, MPI_INT, MASTER_PROCESS, MPI_COMM_WORLD);
        chunk = new int[chunk_size];
        MPI_Scatter(data, chunk_size, MPI_INT, chunk, chunk_size, MPI_INT,
                    MASTER_PROCESS, MPI_COMM_WORLD);

#if DEBUG
       print_info(id, "Got data");
#endif

        quick_sort(chunk, 0, chunk_size - 1);

#if DEBUG
        print_info(id, "Sorted");
#endif
    }

    step = 1;

    while (step < proc_number) {
        if (id % (2 * step) == 0) {
            if (id + step < proc_number) {
                MPI_Recv(&m, 1, MPI_INT, id + step, MASTER_PROCESS,
                         MPI_COMM_WORLD, &status);
                other = new int[m];
                MPI_Recv(other, m, MPI_INT, id + step, MASTER_PROCESS,
                         MPI_COMM_WORLD, &status);

#if DEBUG
                print_info(id, "Got merge data");
#endif

                chunk = merge(chunk, chunk_size, other, m);

#if DEBUG
                print_info(id, "Merged data");
#endif

                delete[] other;
                chunk_size += m;
            }
        } else {
            int near = id - step;
            MPI_Send(&chunk_size, 1, MPI_INT, near, MASTER_PROCESS,
                     MPI_COMM_WORLD);
            MPI_Send(chunk, chunk_size, MPI_INT, near, MASTER_PROCESS,
                     MPI_COMM_WORLD);

#if DEBUG
            print_info(id, "Sent merge data");
#endif

            break;
        }
        step *= 2;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (id == MASTER_PROCESS) {
        stop_time = clock();
        par_time = (stop_time - start_time) / CLOCKS_PER_SEC;
        printf("\nParallel algorithm:\n");
        printf("\nSort completed!\n\tVector size: %d\n\tNumber of processors:"
               " %d\n\tTime: %f secs\n\n", vector_size, proc_number,
               par_time);

        if (vector_size <= 20) {
            printf("Vector after sort: \n");
            print_vector(chunk, vector_size);
        }
    }

    delete[] chunk;

    MPI_Finalize();

    if (id == MASTER_PROCESS) {
        start_time = clock();
        quick_sort(datacp, 0, vector_size - 1);
        stop_time = clock();
        seq_time = (stop_time - start_time) / CLOCKS_PER_SEC;

        printf("\nSequential algorithm:\n");
        printf("\nSort completed!\n\tVector size: %d\n\t"
               "Time: %f secs\n\n", vector_size,
               seq_time);

        if (vector_size <= 20) {
            printf("Vector after sort: \n");
            print_vector(datacp, vector_size);
        }

        printf("\nAcceleration was: %f\n", (seq_time / par_time));

        delete[] datacp;
        delete[] data;
    }
    return 0;
}
