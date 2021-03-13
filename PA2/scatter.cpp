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

    int n_local; 
   
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

        printf("all random numbers are: \n");
        fflush(stdout);
        for(int i = 0; i < N; i++) {
            printf("%f", scatter_values[i]);
        }
      
    }
    printf("\n");

    int n = N;

    int num_res;
    int n_send;

    //set my own local values if I am source_rank
    if(rank == source_rank) {
        num_res = n%p;
        //Check if we are divisible
        if (num_res != 0) {
            n_local = floor(n/p) + 1;
            n_send  = n_local;
            num_res--;
        } else {
            n_local = floor(n/p);
            n_send  = n_local;
            num_res--;
        }          

        local_values = (double *) malloc(n*sizeof(double));
        memcpy(local_values, scatter_values, n_local * sizeof(double));
    }

    MPI_Status stat;
    if(rank == source_rank) {
        for(int j = 1; j < p; j++) {
            //Stop sending the extra numbers
            //printf("res%d\n",num_res);
            if (num_res == 0) {
                //printf("here");
                fflush(stdout);
                n_send--;
            } 
            num_res--;
            
            //printf("nsend %d\n",n_send);
            fflush(stdout);
            //First send the amount of numbers to send. 
            MPI_Send(&n_send, 1, MPI_INT, j, 11, comm);
            //First send the amount of numbers to send. 
            MPI_Send(&scatter_values[j * n_send], n_send, MPI_DOUBLE, j, 12, comm);
        }
    } else {
        //First rec the amount of numbers to rec. 
        MPI_Recv(&n_local, 1, MPI_INT, source_rank, 11, comm, &stat);
        //printf("says to receive %d",n_local);
        fflush(stdout);
        //Allocate the memory
        local_values = (double *) malloc(n*sizeof(double));
        MPI_Recv(local_values, n_local, MPI_DOUBLE, source_rank, 12, comm, &stat);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    printf("proc load as %d\n",n_local);
    fflush(stdout);
    printf("Im rank %d\n",rank);
    printf("my numbers are: \n");
    for(int i = 0; i < n_local; i++) {
        printf("%f", local_values[i]);
    }
    printf("\n");
    

    MPI_Finalize();
}
