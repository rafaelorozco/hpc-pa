#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD;

    int p, rank;
    int source_rank = 0;

    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    printf("Im rank %d\n",rank);

    //Get input arguments. 
    int  N; 
    long C;

    N = atoi(argv[1]);
    C = atoi(argv[2]);

    int proc_load; 
   
    //Dimension of hypercubic communication
    int d = log2(p);

    double scatter_values[N];
    double *local_values;
    //double local_values[4];


    //generate N/p random numbers
    srand48(C+rank);

    

    if(rank == 0) {
        for(int i = 0; i < N; i++) {
            scatter_values[i] = drand48();
        }

        proc_load = N/p;

        printf("all random numbers are: \n");
        fflush(stdout);
        for(int i = 0; i < N; i++) {
            printf("%f", scatter_values[i]);
        }
      
    }
    printf("\n");

    int n = N;

    //set my own local values if I am source_rank
    if(rank == source_rank) {
        local_values = (double *) malloc(proc_load*sizeof(double));
        memcpy(local_values, scatter_values, proc_load * sizeof(double));
    }

    MPI_Status stat;
    if(rank == source_rank) {
        for(int j = 1; j < p; j++) {
            //First send the amount of numbers to send. 
            MPI_Send(&proc_load, 1, MPI_INT, j, 11, comm);

            //First send the amount of numbers to send. 
            MPI_Send(&(scatter_values[j * proc_load]), proc_load, MPI_DOUBLE, j, 12, comm);
        }
    } else {
        //First rec the amount of numbers to rec. 
        MPI_Recv(&proc_load, 1, MPI_INT, source_rank, 11, comm, &stat);
        printf("Got proc load as %d\n",proc_load);
        printf("ITTm rank %d\n",rank);
        //Allocate the memory
        local_values = (double *) malloc(proc_load*sizeof(double));
        MPI_Recv(local_values, proc_load, MPI_DOUBLE, source_rank, 12, comm, &stat);

        printf("ITTm rank %d\n",rank);
        printf("my numbers are: \n");
        for(int i = 0; i < proc_load; i++) {
            printf("%f", local_values[i]);
        }
        printf("\n");

    }

    MPI_Barrier(MPI_COMM_WORLD);

    printf("proc load as %d\n",proc_load);

    printf("Im rank %d\n",rank);
    printf("my numbers are: \n");
    for(int i = 0; i < proc_load; i++) {
        printf("%f", local_values[i]);
    }
    printf("\n");
    

    MPI_Finalize();
}
