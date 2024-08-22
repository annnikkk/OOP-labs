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
    MPI_Init(&argc, &argv);
    start_time = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);//текущий процесс
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);//кол-во процессов
    int *A = (int *) malloc(sizeof(int) * N);
    int *B = (int *) malloc(sizeof(int) * N);
    int offset = N / number_of_processes;
    long long int mul_res;
    if( rank == 0 ){
        filling(A, B, N);
        for( int i = 1; i < number_of_processes; ++i ){
            MPI_Send(A + ( i * offset ), offset, MPI_INT, i, 333, MPI_COMM_WORLD);
            MPI_Send(B, N, MPI_INT, i, 333, MPI_COMM_WORLD);
        }
        mul_res = mul(A, B, offset);
        long long int res = 0;
        for( int i = 1; i < number_of_processes; ++i ){
            MPI_Recv(&res, 1, MPI_LONG_LONG_INT, MPI_ANY_SOURCE, 333, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            mul_res += res;
        }
    } else{
        MPI_Recv(A, offset, MPI_INT, 0, 333, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(B, N, MPI_INT, 0, 333, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        long long int res = mul(A, B, offset);
        MPI_Send(&res, 1, MPI_LONG_LONG_INT, 0, 333, MPI_COMM_WORLD);
    }
    if(rank == 0){
        end_time = MPI_Wtime();
        printf("Время выполнения функции: %f секунд\n", end_time-start_time);
        printf("Сумма: %lld\n", mul_res);
    }
    free(A);
    free(B);
    MPI_Finalize();
    return 0;
}

