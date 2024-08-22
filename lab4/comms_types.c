#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

void MatrixMul(double* local_A, double* local_B, double* local_C, int row_size, int n2, int column_size){
    for(int i = 0; i < row_size; i++){
        for(int j = 0; j < column_size; j++){
            local_C[i * column_size + j] = 0;
            for(int k = 0; k < n2; k++){
                local_C[i * column_size + j] += local_A[i * n2 + k] * local_B[k * column_size + j];
            }
        }
    }
}

void PrintMatrix(double* matrix, int rows, int columns){
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < columns; j++){
            printf("%f ", matrix[i*columns+j]);
        }
        printf("\n");
    }
}

void Filling(double* A, double* B, int N1, int N2, int N3){
    for(int i = 0; i < N1 * N2; ++i){
        A[i] = i + 1;
    }
    for(int i = 0; i < N2 * N3; ++i){
        B[i] = N3* N2 - i;
    }
}

int main(int argc, char** argv){
    double start_time;
    double end_time;
    MPI_Init(&argc, &argv);
    int rank;
    int number_of_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
    int rows = atoi(argv[1]);//кол-во частей в А
    int columns = atoi(argv[2]);//кол-во частей в В
    int N1 = 3456;
    int N2 = 3072;
    int N3 = 3456;
    int row_size = N1 / rows;//кол-во строк в одной части
    int column_size = N3 / columns;//кол-во столбцов в одной части
    double* A;
    double* B;
    double* C;
    double* local_A = (double*) malloc(sizeof(double) * row_size * N2);
    double* local_B = (double*) malloc(sizeof(double) * N2 * column_size);
    double* local_C = (double*) malloc(sizeof(double) * N1 * N3);
    if(rank == 0){
        start_time = MPI_Wtime();
        A = (double*) malloc(sizeof(double) * N1 * N2);
        B = (double*) malloc(sizeof(double) * N2 * N3);
        C = (double*) malloc(sizeof(double) * N1 * N3);
        Filling(A, B, N1, N2, N3);
    }
    int ndims = 2;
    int dims[ndims];
    dims[0] = rows;
    dims[1] = columns;
    int periods[ndims];
    int reorder = 0;
    MPI_Comm DEC_COMM;
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &DEC_COMM);
    
    int subdims[ndims];
    subdims[0] = 0;
    subdims[1] = 1;
    MPI_Comm ROW_COMM;
    MPI_Cart_sub(DEC_COMM, subdims, &ROW_COMM);
    subdims[1] = 0;
    subdims[0] = 1;
    MPI_Comm COLUMN_COMM;
    MPI_Cart_sub(DEC_COMM, subdims, &COLUMN_COMM);
    
    int proc_coords[ndims];
    MPI_Cart_coords(DEC_COMM, rank, ndims, proc_coords);
    if(proc_coords[1] == 0){
        MPI_Scatter(A, row_size*N2, MPI_DOUBLE, local_A, row_size*N2, MPI_DOUBLE, 0, COLUMN_COMM);
    }
    
    MPI_Datatype new_type;
    MPI_Type_vector(N2, column_size, N3, MPI_DOUBLE, &new_type);
    MPI_Type_commit(&new_type);
    if(rank == 0){
        for(int i = 1; i < columns; i++){
            MPI_Send(B + column_size * i, 1, new_type, i, 333, ROW_COMM);
        }
        int k = 0;
        for(int i = 0; i < N2; i++){
            for(int j = 0; j < column_size; j++){
                local_B[k] = B[N3 * i + j];
                k++;
            }
        }
    } else if(proc_coords[0] == 0){
        MPI_Recv(local_B, N2 * column_size, MPI_DOUBLE, 0, 333, ROW_COMM, MPI_STATUS_IGNORE);
    }
    MPI_Bcast(local_A, row_size * N2, MPI_DOUBLE, 0, ROW_COMM);
    MPI_Bcast(local_B, N2 * column_size, MPI_DOUBLE, 0, COLUMN_COMM);
    
    MatrixMul(local_A, local_B, local_C, row_size, N2, column_size);
    
    MPI_Datatype blocks;
    MPI_Type_vector(row_size, column_size, N3, MPI_DOUBLE, &blocks);
    MPI_Type_commit(&blocks);
    
    if(rank != 0){
        MPI_Send(local_C, row_size * column_size, MPI_DOUBLE, 0, 333, DEC_COMM);
    } else {
        for(int i = 0; i < row_size; i++){
            for(int j = 0; j < column_size; j++){
                C[i * N3 + j] = local_C[i * column_size + j];
            }
        }
        
        double* p;
        int this_coords[ndims];
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < columns; j++){
                if(i == 0 && j == 0){
                    continue;
                }
                p = C + (i * N3 * row_size + j * column_size);
                int this_rank = 0;
                this_coords[0] = i;
                this_coords[1] = j;
                MPI_Cart_rank(DEC_COMM, this_coords, &this_rank);
                MPI_Recv(p, 1, blocks, this_rank, 333, DEC_COMM, MPI_STATUS_IGNORE);
            }
        }
        //PrintMatrix(C, N1, N3);
    }
    if(rank == 0){
        end_time = MPI_Wtime();
        printf("Время выполнения функции: %f секунд\n", end_time - start_time);
        free(A);
        free(B);
        free(C);
    }
    free(local_A);
    free(local_B);
    free(local_C);
    MPI_Type_free(&new_type);
    MPI_Type_free(&blocks);
    MPI_Finalize();
    return 0;
}
