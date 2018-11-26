/* quicksort */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <limits.h>

#define N 1000000

void print_info(int id, char *m);
void print_vector(int *v, int vector_size, int id);
int *merge(int *v1, int n1, int *v2, int n2);
void swap(int *v, int i, int j);
void quick_sort(int *v, int left, int right);

double start_time, stop_time;

void print_info(int id, char *m)
{
    printf("%d: %s %f secs\n", id, m, (clock() - start_time) / CLOCKS_PER_SEC);
}

void print_vector(int *v, int vector_size, int id)
{
    int i;

    printf("%d: ", id);

    for (i = 0; i < vector_size; i++)
        printf("%d ", v[i]);
    putchar('\n');
}

int *merge(int *v1, int n1, int *v2, int n2)
{
    int *result;
    int i = 0,
        j = 0,
        k = 0;

    result = (int *)malloc((n1 + n2) * sizeof(int));

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

void swap(int *v, int i, int j)
{
    int t;

    t = v[i];
    v[i] = v[j];
    v[j] = t;
}

void quick_sort(int *v, int left, int right)
{
    int i,last;

    if (left >= right)
        return;

    swap(v, left, (left + right) / 2);
    last = left;
    for (i = left + 1; i <= right; i++)
        if(v[i] < v[left])
            swap(v, ++last, i);

    swap(v, left, last);
    quick_sort(v, left, last-1);
    quick_sort(v, last + 1, right);
}

int main(int argc, char **argv)
{
    MPI_Status status;
    int m, vector_size = N;
    int id, proc_number;
    int chunk_size;
    int *chunk;
    int *other;
    int *data;
    int step;
    int i;

    if (argc > 1) {
        vector_size = atoi(argv[1]);
    }

    start_time = clock();

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_number);

    print_info(id, "MPI setup complete");

    if (id == 0) {
        int extra_elements;
        srandom(clock());

        chunk_size = vector_size / proc_number;
        extra_elements = vector_size % proc_number;

        /* alloc memory for vector
         * increase the size of the array so that
         * it is divided without a balance by the number of processors
         */
        data = (int *)malloc((vector_size + chunk_size - extra_elements) * sizeof(int));

        /* fill main elements of vector
         * random() + random() can give overflow int,
         * thus we get negative numbers in the vector.
         * it needs to be done because
         * random() returns non-negative numbers
         */
        for (i = 0; i < vector_size; i++)
            data[i] = random() + random();

        /* fill extra elements of vector
         * fill in the minimum number so that after sorting
         * you know which extra elements need to be thrown out
         */
        if (extra_elements != 0) {
            for (i = vector_size; i < vector_size + chunk_size - extra_elements; i++)
                data[i] = INT_MIN;
            chunk_size++;
        }

        print_info(id, "generated the random numbers");

        MPI_Bcast(&chunk_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        chunk = (int *)malloc(chunk_size * sizeof(int));
        MPI_Scatter(data, chunk_size, MPI_INT, chunk, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

        print_info(id, "scattered data");

        quick_sort(chunk, 0, chunk_size - 1);

        print_info(id, "sorted");
    } else {
        MPI_Bcast(&chunk_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        chunk = (int *)malloc(chunk_size * sizeof(int));
        MPI_Scatter(data, chunk_size, MPI_INT, chunk, chunk_size, MPI_INT, 0, MPI_COMM_WORLD);

        print_info(id, "got data");

        quick_sort(chunk, 0, chunk_size - 1);

        print_info(id, "sorted");
    }

    step = 1;

    while (step < proc_number) {
        if (id % (2 * step) == 0) {
            if (id + step < proc_number) {
                MPI_Recv(&m, 1, MPI_INT, id + step, 0, MPI_COMM_WORLD, &status);
                other = (int *)malloc(m * sizeof(int));
                MPI_Recv(other, m, MPI_INT, id + step, 0, MPI_COMM_WORLD, &status);
                print_info(id, "got merge data");
                chunk = merge(chunk, chunk_size, other, m);
                print_info(id, "merged data");
                chunk_size += m;
            } 
        } else {
            int near = id - step;
            MPI_Send(&chunk_size, 1, MPI_INT, near, 0, MPI_COMM_WORLD);
            MPI_Send(chunk, chunk_size, MPI_INT, near, 0, MPI_COMM_WORLD);
            print_info(id, "sent merge data");
            break;
        }
        step *= 2;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (id == 0) {
        FILE *fout;

        stop_time = clock();
        printf("\nSort completed!\n\tVector size: %d\n\tNumber of processors: %d\n\tTime: %f secs\n\n", vector_size, proc_number, (stop_time - start_time) / CLOCKS_PER_SEC);

        print_info(id, "opening out file");
        fout = fopen("result", "w");
        for (i = 0; i < vector_size; i++)
            fprintf(fout, "%d\n", chunk[(chunk_size - vector_size) + i]);

        fclose(fout);
        print_info(id, "closed out file");
    }

    MPI_Finalize();

    return 0;
}