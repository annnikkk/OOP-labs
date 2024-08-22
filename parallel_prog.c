#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

void MatrixOnVector(double* matrix, double* vec, double* res, int size, int part_size){
    for(int i = 0; i < part_size; i++){
        res[i] = 0;
        for(int j = 0; j < size; j++){
            res[i] += matrix[i*size+j] * vec[j];
        }
    }
}

void VecMinusVec(double* vec1, double* vec2, double* res, int size){
    for(int i = 0; i < size; i++){
        res[i] = vec1[i] - vec2[i];
    }
}

double ScalarMul(double* vec1, double* vec2, int size){
    double res = 0;
    for(int i = 0; i < size; i++){
        res += vec1[i] * vec2[i];
    }
    return res;
}

void VectorOnScalar(double* vec, double* res, double scalar, int size){
    for(int i = 0; i < size; i++){
        res[i] = vec[i] * scalar;
    }
}

double Norm(double* vec, int size){
    double res = 0;
    for(int i = 0; i < size; i++){
        res += vec[i] * vec[i];
    }
    return res;
}

void TestFilling(double* matrix, double* vec, double* x_n, int size){
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            if(i == j){
                matrix[i*size+j] = 2.0;
            } else {
                matrix[i*size+j] = 1.0;
            }
        }
        vec[i] = size + 1;
        x_n[i] = 0;
    }
}

void Filling(double* matrix, double* vec, double* x_n, int size){
    for(int i = 0; i < size; i++){
        srand(i);
        for(int j = 0; j < size; j++){
            if(i == j){
                matrix[i*size+j] = (rand() % 1000) + size;
            } else {
                matrix[i*size+j] = size/2;
            }
        }
        x_n[i] = 0;
    }

    double* u = (double*)malloc(sizeof(double) * size);

    for(int i = 0; i < size; ++i) {
        u[i] = sin(2 * M_PI * i / size);
    }

    MatrixOnVector(matrix, u, vec, size, size);
    free(u);
}

int main(int argc, char **argv){
    double start_time;
    double end_time;
    MPI_Init(&argc, &argv);
    int rank;
    int number_of_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
    int N = 6144;
    int part_size = N/number_of_processes;
    double epsilon = 0.00001;
    //double o = 0.0000001;
    double tao = 0.1;
    double ending_coeff = 0;

    double* x_n = (double *)malloc(sizeof(double) * N);
    double* y = (double *)malloc(sizeof(double) * N);
    double* A;
    double* b;

    double* tmp = (double *)malloc(sizeof(double) * part_size);

    double* one = (double *)malloc(sizeof(double) * 2);
    double* two = (double *)malloc(sizeof(double) * 2);


    if(rank == 0){
        start_time = MPI_Wtime();
        A = (double*)malloc(sizeof(double) * N*N);
        b = (double *)malloc(sizeof(double) * N);
        Filling(A, b, x_n, N);
    }

    double* b_part = (double *)malloc(sizeof(double) * part_size);
    MPI_Scatter(b, part_size, MPI_DOUBLE, b_part, part_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double* A_part = (double *) malloc(sizeof(double) * part_size*N);
    MPI_Scatter(A, part_size * N, MPI_DOUBLE, A_part, part_size * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(x_n, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double b_norm = 0;
    double part_of_b_norm = Norm(b_part, part_size);
    MPI_Allreduce(&part_of_b_norm, &b_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    b_norm = sqrt(b_norm);

    for(int i = 0; i < 50000; i++){
        MatrixOnVector(A_part, x_n, tmp, N, part_size);
        VecMinusVec(tmp, b_part, tmp, part_size);
        double y_norm = 0;
        double part_of_y_norm = Norm(tmp, part_size);
        MPI_Allreduce(&part_of_y_norm, &y_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        y_norm = sqrt(y_norm);
        double index = y_norm / b_norm;
        /*if(rank == 0){
            printf("index = %lf\n", index);
        }*/
        if(index < epsilon){
            ending_coeff ++;
            if(ending_coeff == 3){
                if(rank == 0){
                    printf("iterations: %d\n", i);
                }
                break;
            }
        } else ending_coeff = 0;
        MPI_Allgather(tmp, part_size, MPI_DOUBLE, y, part_size, MPI_DOUBLE, MPI_COMM_WORLD);
        MatrixOnVector(A_part, y, tmp, N, part_size);
        double tao1 = 0;
        double tao2 = 0;
        double part1_tao = 0;
        double part2_tao = 0;
        part1_tao = ScalarMul(y+part_size*rank, tmp, part_size) ;
        part2_tao = ScalarMul(tmp, tmp, part_size);
        MPI_Allreduce(&part1_tao, &tao1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&part2_tao, &tao2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        tao = tao1/tao2;
        VectorOnScalar(y+part_size*rank, tmp, tao, part_size);
        VecMinusVec(x_n+part_size*rank, tmp, tmp, part_size);
        MPI_Allgather(tmp, part_size, MPI_DOUBLE, x_n, part_size, MPI_DOUBLE, MPI_COMM_WORLD);
    }
    if (rank == 0){
        end_time = MPI_Wtime();
        printf("Время выполнения функции: %f секунд\n", end_time-start_time);
        free(A);
        free(b);
    }
    free(b_part);
    free(A_part);
    free(x_n);
    free(y);
    free(tmp);
    MPI_Finalize();
    return 0;
}