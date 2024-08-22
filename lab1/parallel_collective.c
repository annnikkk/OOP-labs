#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define N 120000

void filling(int *A, int *B, int size){
    for( int i = 0; i < size; ++i ){
        A[i] = (rand() % 1000);
        B[i] = (rand() % 1000);
    }
}

long long int mul(int *A, int *B, int size){
    long long int res = 0;
    for( int i = 0; i < size; ++i ){
        for( int j = 0; j < N; ++j ){
            res += A[i] * B[j];
        }
    }
    return res;
}

int main(int argc, char **argv){
    int rank;
    int number_of_processes;
    double start_time;
    double end_time;
    int* A;
    int* B;
    int* vec;
    MPI_Init(&argc, &argv);
    start_time = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);//текущий процесс
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);//кол-во процессов
    int offset = N / number_of_processes;
    long long int mul_res = 0;
    B = (int *) malloc(sizeof(int) * N);
    if( rank == 0 ){
        A = (int *) malloc(sizeof(int) * N);
        filling(A, B, N);
    }
    vec = (int *) malloc(sizeof(int) * N);

    MPI_Scatter(A, offset, MPI_INT, vec, offset, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(B, N, MPI_INT, 0, MPI_COMM_WORLD);

    long long int res = mul(vec, B, offset);

    MPI_Reduce(&res, &mul_res, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(rank == 0){
        end_time = MPI_Wtime();
        printf("Время выполнения функции: %f секунд\n", end_time-start_time);
        printf("Сумма: %lld\n", mul_res);
        free(A);
    }
    free(B);
    free(vec);
    MPI_Finalize();
    return 0;
}

