#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>

#define NUMBER_OF_TASKS 2048
#define NUMBER_OF_ITERATIONS 8
#define L 2048

#define NEED_NEW_TASK 333
#define NEED_NUMBER_OF_TASKS 444
#define SENDING_TASKS 555

#define FINISH_CONDITION 777

int rank;
int number_of_processes;

int remaining_tasks;
int done_tasks;

double global_res;
int* tasks;

pthread_mutex_t mutex;
pthread_t threads[2];

void FillingTasks(int iter_counter){
    for(int i = 0; i < NUMBER_OF_TASKS; i++){
        tasks[i] = abs(50 - i % 100) * abs(rank - (iter_counter % number_of_processes)) * L;
    }
}

void DoingTasks(){
    for(int i = 0; i < remaining_tasks; i++){
        pthread_mutex_lock(&mutex);
        int weight = tasks[i];
        pthread_mutex_unlock(&mutex);
        for(int j = 0; j < weight; j++){
            global_res += sin(j);
        }
        pthread_mutex_lock(&mutex);
        done_tasks++;
        pthread_mutex_unlock(&mutex);
    }
    pthread_mutex_lock(&mutex);
    remaining_tasks = 0;
    pthread_mutex_unlock(&mutex);
}

void* Work(){
    double local_start_time;
    double local_end_time;
    double iteration_time;
    double min_time, max_time;
    tasks = calloc(NUMBER_OF_TASKS, sizeof(int));
    for(int i = 0; i < NUMBER_OF_ITERATIONS; i++){
        pthread_mutex_lock(&mutex);
        FillingTasks(i);
        remaining_tasks = NUMBER_OF_TASKS;
        done_tasks = 0;
        pthread_mutex_unlock(&mutex);
        local_start_time = MPI_Wtime();
        DoingTasks();
        int number_of_extra_tasks = 0;
        int exit_flag = 0;
        while(1){
            if(exit_flag == 1) break;
            int number_of_free_ranks = 0;
            for(int j = 0; j < number_of_processes; j++){
                if(j == rank) continue;
                MPI_Send(&rank, 1, MPI_INT, j,  NEED_NEW_TASK, MPI_COMM_WORLD);
                MPI_Recv(&number_of_extra_tasks, 1, MPI_INT, j, NEED_NUMBER_OF_TASKS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if(number_of_extra_tasks == 0) number_of_free_ranks++;
                if(number_of_free_ranks == number_of_processes-1){
                    exit_flag = 1;
                    break;
                }
                if(number_of_extra_tasks == 0) continue;
                MPI_Recv(tasks, number_of_extra_tasks, MPI_INT, j, SENDING_TASKS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                pthread_mutex_lock(&mutex);
                remaining_tasks = number_of_extra_tasks;
                pthread_mutex_unlock(&mutex);
                DoingTasks();
            }
        }
        local_end_time = MPI_Wtime();
        iteration_time = local_end_time - local_start_time;
        MPI_Allreduce(&iteration_time,&max_time,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        MPI_Allreduce(&iteration_time,&min_time,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        if(rank == 0){
            double disbalance = max_time - min_time;
            double part_disbalance = (disbalance / max_time) * 100;
            printf("iteration: %d, iteration_time: %lf, max_time: %lf, min_time: %lf, disbalance: %lf, part_disbalance: %lf\n", i, iteration_time, max_time, min_time, disbalance, part_disbalance);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    int condition = FINISH_CONDITION;
    MPI_Send(&condition, 1, MPI_INT, rank, NEED_NEW_TASK, MPI_COMM_WORLD);
    free(tasks);
    pthread_exit(NULL);
}

void* Communication(){
    int sending_rank;
    int extra_tasks = 0;
    int tasks_to_send = 0;
    while(1){
        MPI_Recv(&sending_rank, 1, MPI_INT, MPI_ANY_SOURCE,  NEED_NEW_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(sending_rank == FINISH_CONDITION) break;
        pthread_mutex_lock(&mutex);
        extra_tasks = remaining_tasks - done_tasks;
        tasks_to_send = extra_tasks / (number_of_processes - 1);
        if(tasks_to_send > 0 && remaining_tasks > 0){
            remaining_tasks -= tasks_to_send;
            MPI_Send(&tasks_to_send, 1, MPI_INT, sending_rank, NEED_NUMBER_OF_TASKS, MPI_COMM_WORLD);
            MPI_Send(tasks+(remaining_tasks-tasks_to_send), tasks_to_send, MPI_INT, sending_rank, SENDING_TASKS, MPI_COMM_WORLD);
        } else {
            tasks_to_send = 0;
            MPI_Send(&tasks_to_send, 1, MPI_INT, sending_rank, NEED_NUMBER_OF_TASKS, MPI_COMM_WORLD);
        }
        pthread_mutex_unlock(&mutex);
    }
    pthread_exit(NULL);
}

int main(int argc, char* argv[]){
    double start_time;
    double end_time;
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if(provided != MPI_THREAD_MULTIPLE){
        perror("MPI implementation doesn't provide required threading level");
        MPI_Finalize();
        return 1;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
    
    pthread_mutex_init(&mutex, NULL);
    pthread_attr_t attrs;
    if (pthread_attr_init(&attrs) != 0){
        printf("error in creating attrs\n");
        MPI_Finalize();
        return 0;
    }
    start_time = MPI_Wtime();
    if (pthread_create(&threads[0], &attrs, Communication, NULL) != 0){
        printf("error in creating communicating thread\n");
        MPI_Finalize();
        return 0;
    }
    
    if (pthread_create(&threads[1], &attrs, Work, NULL) != 0){
        printf("error in creating working thread\n");
        MPI_Finalize();
        return 0;
    }
    for (int i = 0; i < 2; i++){
        if (pthread_join(threads[i],NULL) != 0){
            printf("Error joining thread %d\n",i);
            MPI_Finalize();
            return 0;
        }
    }
    end_time = MPI_Wtime();
    if(rank == 0){
        printf("TOTAL TIME: %lf\n", end_time - start_time);
    }
    pthread_attr_destroy(&attrs);
    pthread_mutex_destroy(&mutex);
    MPI_Finalize();
    return 0;
}
