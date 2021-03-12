/*
 * CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 2
 * 
 *  MPI polynomial evaluation algorithm function implementations go here
 * 
 */

#include "mpi_evaluator.h"
#include "const.h"

void scatter(const int n, double* scatter_values, int &n_local, double* &local_values, int source_rank, const MPI_Comm comm){

	int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    //set my own local values if I am source_rank
    if(rank == source_rank) {
    	n_local = n/p;
    	local_values = (double *) malloc(n*sizeof(double))
    	memcpy(local_values, scatter_values, n_local * sizeof(double));
    }

    MPI_Status stat;
    if(rank == source_rank) {
        for(int j = 1; j < p; j++) {
            //First send the amount of numbers to send. 
            MPI_Send(&n_local, 1, MPI_INT, j, 11, comm);

            //First send the amount of numbers to send. 
            MPI_Send(&scatter_values[j * n_local], n_local, MPI_DOUBLE, j, 12, comm);
        }
    } else {
        //First rec the amount of numbers to rec. 
        MPI_Recv(&n_local, 1, MPI_INT, source_rank, 11, comm, &stat);

        //Allocate the memory
        local_values = (double *) malloc(n*sizeof(double))
        MPI_Recv(local_values, n_local, MPI_DOUBLE, source_rank, 12, comm, &stat);
    }
}

double broadcast(double value, int source_rank, const MPI_Comm comm){
	int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    //Dimension of hypercubic communication
    int d = log2(p);

    //check if p = 2^d
    if((1 << (d)) != p ) {
    	d++;
    }

    MPI_Status stat;
    int flip = 1 << (d-1);
    int mask = flip - 1;
    for (int j = d-1; j>=0 ;j--) {
	    if((rank & mask) == 0) { //check if I am participating in this communic
	    	if((rank & flip) == 0) { //check if I send
	    		if((rank ^ flip) < p) {
	    			
	    			MPI_Send(&value, 1, MPI_DOUBLE, rank ^ flip, 11, comm);
	    		}
	    	} else {
	    		if((rank ^ flip) < p) {
	    		    //printf("receive my num from %i\n", rank ^ flip);
	    			MPI_Recv(&value, 1, MPI_DOUBLE, rank ^ flip, 11, comm, &stat);
	    		}
	    	}     
	    } 
	    mask = mask >> 1;
	    flip = flip >> 1;
	}
    return value;
}

void parallel_prefix(const int n, const double* values, double* prefix_results, const int OP, const MPI_Comm comm){
    //Implementation

}

double mpi_poly_evaluator(const double x, const int n, const double* constants, const MPI_Comm comm){
    //Implementation

    return 0;
}
