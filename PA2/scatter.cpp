#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD;

    int p, rank;

    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    printf("all random numbers are: \n");
    fflush()
   
    //Get input arguments. 
    int  N; 
    long C;

    N = atoi(argv[1]);

    int proc_load = N/p;
   
    //Dimension of hypercubic communication
    int d = log2(p);

    double all_nums[N];
    double my_nums[proc_load];

    for(int i = 0; i < N; i++) {
        all_nums[i] = drand48();
    }

    if(rank == 0) {
        printf("all random numbers are: \n");
        fflush(stdout);
        for(int i = 0; i < N; i++) {
            printf("%f", all_nums[i]);
        }
        //Copy the n/p numbers into the rank 0 array
        memcpy(my_nums, all_nums, proc_load * sizeof(double));
    }
    printf("\n");

    MPI_Status stat;
    if(rank == 0) {
        for(int j = 1; j < p; j++) {
            //printf("send my numbers to %i\n",  j);
            MPI_Send(&all_nums[j * proc_load], proc_load, MPI_DOUBLE, j, 11, MPI_COMM_WORLD);
            
        }
    } else {
        //printf("receive my numbers from 0");
        MPI_Recv(&my_nums, proc_load, MPI_DOUBLE, 0, 11, MPI_COMM_WORLD,&stat);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    printf("Im rank %d\n",rank);
    printf("my numbers are: \n");
    for(int i = 0; i < proc_load; i++) {
        printf("%f", my_nums[i]);
    }
    printf("\n");
    

    MPI_Finalize();
}
