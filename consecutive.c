#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

#define N 120000

void filling(int* A, int* B){
    for(int i = 0; i < N; ++i){
        A[i] = (rand() % 1000);
        B[i] = (rand() % 1000);
    }
}

long long int mul(int* A, int* B){
    long long int res = 0;
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < N; ++j){
            res += A[i] * B[j];
        }
    }
    return res;
}

int main(int argc, char **argv){
    double start_time;
    double end_time;
    MPI_Init(&argc, &argv);
    start_time = MPI_Wtime();
    int* A = (int*) malloc(sizeof(int)*N);
    int* B = (int*) malloc(sizeof(int)*N);
    filling(A, B);
    long long int mul_res = mul(A, B);
    end_time = MPI_Wtime();
    printf("Время выполнения функции: %f секунд\n", end_time-start_time);
    printf("Сумма: %lld\n", mul_res);
    free(A);
    free(B);
    MPI_Finalize();
    return 0;
}