#include <stdio.h>
#include <stdlib.h>

#include "mpi.h"

int NumberOfNeighbors(char* current_era, int coordX, int coordY, int width){
    int number_of_neighbors = 0;
    for(int i = -1; i <= 1; i++){
        for(int j = -1; j <= 1; j++){
            if(i == 0 && j ==0) continue;
            int new_coordX = coordX + j;
            int new_coordY = coordY + i;
            if(new_coordX < 0){
                new_coordX = width - 1;
            }
            if(new_coordX > width - 1){
                new_coordX = 0;
            }
            if(current_era[new_coordY*width + new_coordX] == 1){
                number_of_neighbors += 1;
            }
        }
    }
    return number_of_neighbors;
}

void MiddleCalc(char* current_era, char* old_era, int width, int height){//здесь height это высота всей части включая принимаемые строки
    for(int i = 2; i < height-2; i++){
        for(int j = 0; j < width; j++){
            int number_of_neighbors = NumberOfNeighbors(old_era, j, i, width);
            if(old_era[i*width+j] == 0){
                if(number_of_neighbors == 3){
                    current_era[i*width+j] = 1;
                } else {
                    current_era[i*width+j] = 0;
                }
            } else {
                if(number_of_neighbors < 2 || number_of_neighbors > 3){
                    current_era[i*width+j] = 0;
                } else {
                    current_era[i*width+j] = 1;
                }
            }
        }
    
    }
}

void TopCalc(char* current_era, char* old_era, int width){
    for(int j = 0; j < width; j++){
        int number_of_neighbors = NumberOfNeighbors(old_era, j, 1, width);
        if(old_era[width+j] == 0){
            if(number_of_neighbors == 3){
                current_era[width+j] = 1;
            } else {
                current_era[width+j] = 0;
            }
        } else {
            if(number_of_neighbors < 2 || number_of_neighbors > 3){
                current_era[width+j] = 0;
            } else {
                current_era[width+j] = 1;
            }
        }
    }
}

void LowCalc(char* current_era, char* old_era, int width, int height){
    for(int j = 0; j < width; j++){
        int number_of_neighbors = NumberOfNeighbors(old_era, j, height-2, width);
        int coord = width*(height-2)+j;
        if(old_era[coord] == 0){
            if(number_of_neighbors == 3){
                current_era[coord] = 1;
            } else {
                current_era[coord] = 0;
            }
        } else {
            if(number_of_neighbors < 2 || number_of_neighbors > 3){
                current_era[coord] = 0;
            } else {
                current_era[coord] = 1;
            }
        }
    }
}

void ComparingEras(char* eras[2048], char* current_era, char* matches, int current_number, int width, int height){
    char condition = 1;
    char* old_era;
    for(int k = 0; k <= current_number; k++){
        old_era = eras[k];
        for(int i = width; i < (height-1)*width; i++){
            if(old_era[i] != current_era[i]){
                condition = 0;
                break;
            }
        }
        matches[k] = condition;
    }
}

int main(int argc, char** argv){
    double start_time;
    double end_time;
    int width = atoi(argv[1]);
    int height = atoi(argv[2]);
    MPI_Init(&argc, &argv);
    int rank;
    int number_of_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
    int part_height = (height/number_of_processes) +2;
    char* part_of_current_era = (char*) calloc(width*part_height, 1);
    char* eras[2048];
    char* current_era = NULL;
    int top_process = rank-1;
    int low_process = rank + 1;
    if(rank == number_of_processes-1){
        low_process = 0;
    }
    if(rank == 0){
        top_process = number_of_processes-1;
        current_era = (char*) calloc(width*height, 1);
        current_era[width + 1] = 1;
        current_era[width * 2 + 2] = 1;
        current_era[3 * width] = 1;
        current_era[3 * width + 1] = 1;
        current_era[3 * width + 2] = 1;
        start_time = MPI_Wtime();
    }
    MPI_Scatter(current_era, (height/number_of_processes)*width, MPI_CHAR, part_of_current_era+width, (height/number_of_processes)*width,MPI_CHAR, 0, MPI_COMM_WORLD);
    
    

    char matches[2048];
    int i;
    for(i = 0; i < 2048; i++){
        MPI_Request r1, r2, r3, r4;
        MPI_Isend(part_of_current_era+width, width, MPI_CHAR, top_process, 333, MPI_COMM_WORLD,  &r1);//отправка второй
        MPI_Isend(part_of_current_era+width*(part_height-2), width, MPI_CHAR, low_process, 334, MPI_COMM_WORLD,  &r2);//отправка предпоследней
        
        MPI_Irecv(part_of_current_era, width, MPI_CHAR, top_process, 334, MPI_COMM_WORLD, &r3);//получаю первую
        MPI_Irecv(part_of_current_era+width*(part_height-1), width, MPI_CHAR, low_process, 333, MPI_COMM_WORLD, &r4);//получаю последнюю
        
        char* new_era = (char*) calloc(width*part_height, 1);
        MiddleCalc(new_era, part_of_current_era, width, part_height);
        
        MPI_Wait(&r3, MPI_STATUS_IGNORE);
        TopCalc(new_era, part_of_current_era, width);
        
        MPI_Wait(&r4, MPI_STATUS_IGNORE);
        LowCalc(new_era, part_of_current_era, width, part_height);
        
        MPI_Wait(&r1, MPI_STATUS_IGNORE);
        MPI_Wait(&r2, MPI_STATUS_IGNORE);
        
        eras[i] = part_of_current_era;
        char win_condition = 0;
        ComparingEras(eras, new_era, matches, i, width, part_height);
        MPI_Allreduce(MPI_IN_PLACE, &matches, i, MPI_CHAR, MPI_LAND, MPI_COMM_WORLD);
        for(int k = 0; k < i; k++){
            if(matches[k] == 1){
                win_condition = 1;
                break;
            }
        }
        if(win_condition == 1) break;
        part_of_current_era = new_era;
    }
    
    if(rank == 0){
        end_time = MPI_Wtime();
        printf("iterations: %d, time: %lf\n", i, end_time- start_time);
        free(current_era);
    }
    
    for(int k = 0; k <= i; k++){
        free(eras[k]);
    }
    
    MPI_Finalize();
}
