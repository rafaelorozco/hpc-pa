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
    	local_values = (double *) malloc(n*sizeof(double));
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
        local_values = (double *) malloc(n*sizeof(double));
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

    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);
    MPI_Barrier(comm);

    double local_prefix_end = (OP == PREFIX_OP_PRODUCT)? 1 : 0;

    double local_total_end = prefix_results[n-1];
    double conj_total;

    MPI_Status stat;

    for(int j = 0; j < log2(p); j++) {
            int conj_rank = rank ^ (1 << j);
            MPI_Request req;
            MPI_Irecv(&conj_total, 1, MPI_DOUBLE, conj_rank, 11, MPI_COMM_WORLD, &req);
            MPI_Send(&local_total_end, 1, MPI_DOUBLE, conj_rank, 11, MPI_COMM_WORLD);
            MPI_Wait(&req, &stat);
            local_total_end = (OP == PREFIX_OP_PRODUCT) ? (local_total_end * conj_total) : (local_total_end + conj_total);
            if(rank > conj_rank) local_prefix_end = (OP == PREFIX_OP_PRODUCT) ? (local_prefix_end * conj_total) : (local_prefix_end + conj_total);
    }

    for (int i=1; i<n; i++){
            prefix_results[i] =  (OP == PREFIX_OP_PRODUCT)? (prefix_results[i] * local_prefix_end) : (prefix_results[i] + local_prefix_end);
    }

}

double mpi_poly_evaluator(const double x, const int n, const double* constants, const MPI_Comm comm){
    //Implementation
        int p, rank;
        MPI_Comm_size(comm, &p);
        MPI_Comm_rank(comm, &rank);

        int source_rank;
        MPI_Comm_rank(comm, &source_rank);

        int p1;
        int n1=n;
        
        if (n%p != 0) n1 = ceil(n/p);           // in case n is not divisible by p
        p1 = pow(2,ceil(log2(double(p))));      // in case p is not a power of 2

        double* values = (double*) malloc(n*sizeof(double));
        double* local_prefix = (double*) malloc(n*sizeof(double));

        for (int i=0;i<n;i++){
            if (rank == 0 && i == 0) values[i] = 1;
            else values[i] = x;
        }

        // sequential prefix in each processor
        local_prefix[0] = values[0];
        for (int i=1; i<n; i++){
                local_prefix[i] =  local_prefix[i-1] * values[i];
        }
 
        parallel_prefix(n, values, local_prefix, PREFIX_OP_PRODUCT, comm);

        MPI_Barrier(MPI_COMM_WORLD);

        double sum = 0;
        double other_sum;
        for(int i = 0; i < n; i++) {
            sum += constants[i]*local_prefix[i];
        }
        MPI_Status stat;
        for(int j = 0; j < log2(p); j++) {
            int b2_pow_j = 1 << j;
            if((rank & b2_pow_j) != 0) {
                MPI_Send(&sum, 1, MPI_DOUBLE, rank ^ b2_pow_j, 11, MPI_COMM_WORLD);
            }
            else {
                MPI_Recv(&other_sum,1,MPI_DOUBLE,rank^b2_pow_j,11,MPI_COMM_WORLD,&stat);
                sum += other_sum;
            }
        }
        if (rank==0) return sum;
}
