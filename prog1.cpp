#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]) {

    
    MPI_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD;

    int p, rank;

    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    printf(" hello from %i%i\n", rank, p);
    //Get input arguments. 
    int N; 
    long C;

    N = atoi(argv[1]);
    C = atol(argv[2]);

    int proc_load = N/p;
    int d = log2(p);

    //generate N/p random numbers
    srand48(C+rank);


    printf(" My random numbers are: \n");
    double sum = 0;
    double other_sum;

    for(int i = 0; i < proc_load; i++) {
        double curr_num = drand48();
        printf("%f", curr_num);
        sum += curr_num;
    }
    
    printf("My sum is %f\n",sum);
    printf("d is %i\n", d);
    fflush(stdout);

    MPI_Status stat;
    for(int j = 0; j < d; j++) {
        int b2_pow_j = 1 << j;
        if((rank & b2_pow_j) != 0) {
            printf("sebd my sum to %i\n", rank ^ (1 << j));
            MPI_Send(&sum, 1, MPI_DOUBLE, rank ^ b2_pow_j, 11, MPI_COMM_WORLD);
        }
        else {
            printf("receive my sum from %i\n", rank ^ (1 << j));
            MPI_Recv(&other_sum,1,MPI_DOUBLE,rank^b2_pow_j,11,MPI_COMM_WORLD,&stat);
            sum += other_sum;
        }

    }

    if(rank == 0) {
        printf("Final sum is %f\n",sum);
    }


    MPI_Finalize();
}
