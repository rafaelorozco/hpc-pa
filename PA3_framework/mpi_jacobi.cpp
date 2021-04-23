/**
 * @file    mpi_jacobi.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements MPI functions for distributing vectors and matrixes,
 *          parallel distributed matrix-vector multiplication and Jacobi's
 *          method.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include "mpi_jacobi.h"
#include "jacobi.h"
#include "utils.h"

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <vector>

/*
 * TODO: Implement your solutions here
 */

#define NDIM 2


void distribute_vector(const int n, double* input_vector, double** local_vector, MPI_Comm comm)
{
    // TODO
}


// gather the local vector distributed among (i,0) to the processor (0,0)
void gather_vector(const int n, double* local_vector, double* output_vector, MPI_Comm comm)
{
    // TODO
}

void distribute_matrix(const int n, double* input_matrix, double** local_matrix, MPI_Comm comm)
{
    // TODO
}


void transpose_bcast_vector(const int n, double* col_vector, double* row_vector, MPI_Comm comm)
{
    // TODO
}


void distributed_matrix_vector_mult(const int n, double* local_A, double* local_x, double* local_y, MPI_Comm comm)
{
    int p, rank;
    int wrank, crank;
    int my_row_rank;
    int my_col_rank;

    int coords[NDIM];
    int free_coords[NDIM];

    MPI_Comm comm1D_row;
    MPI_Comm comm1D_col;

    MPI_Comm_size(comm, &p);
    //MPI_Comm_rank(comm, &rank);

    MPI_Comm_rank(comm, &crank);
    MPI_Cart_coords(comm, crank, NDIM, coords);

   
    //How many elements from x are here?
    int msg_size;
    int q = sqrt(p);
    if(coords[0] < (n % (int) sqrt(p))){
        msg_size = ceil(n/q);
    } else {
        msg_size = floor(n/q);
    }

    //make row communicator 
    free_coords[0] = 0; /* rows */; free_coords[1] = 1; /* cols */
    MPI_Cart_sub(comm, free_coords, &comm1D_row);

    //make row communicator 
    free_coords[0] = 1; /* rows */; free_coords[1] = 0; /* cols */
    MPI_Cart_sub(comm, free_coords, &comm1D_col);

    //get my rank in this row
    MPI_Comm_rank(comm1D_row, &my_row_rank);

    //get my rank in this col
    MPI_Comm_rank(comm1D_col, &my_col_rank);

    printf("I'm %d %d and my row rank is %d and my msgsize %d and element \n", coords[0], coords[1], my_row_rank,msg_size);

    // //////Need to tranpose vector onto grid 
    //1.) (i,0) processors send their local element send to diagonal element. (i,0)
    MPI_Status stat;
    if(my_row_rank == 0) {
        //printf("I'm %d %d", coords[0], coords[1]);
        //printf("\n ");
        local_x[0] = 1.;
        printf("and I am sending to my element %f com rank %d",local_x,coords[0]);

        MPI_Send(&local_x, msg_size, MPI_DOUBLE, coords[0], 11, comm1D_row);
    }
    if(my_row_rank == coords[0]) {
        printf("and I am waiting ");
        MPI_Recv(&local_x, msg_size, MPI_DOUBLE, 0, 11, comm1D_row, &stat);

        //2.) (i,i) processors send their local element to all others in col communicator. 
    }
    if(coords[1] == coords[0]) {
        MPI_Bcast(&local_x, msg_size, MPI_DOUBLE, coords[0], comm1D_col);
    }
    //wait
    MPI_Barrier(comm);

    printf("I'm %d %d and my row rank is %d and my msgsize %d and element %f \n", coords[0], coords[1], my_row_rank,msg_size, local_x);




    //Do local multiplication using sequential matvec???
}

// Solves Ax = b using the iterative jacobi method
void distributed_jacobi(const int n, double* local_A, double* local_b, double* local_x,
                MPI_Comm comm, int max_iter, double l2_termination)
{
    // TODO
}


// wraps the distributed matrix vector multiplication
void mpi_matrix_vector_mult(const int n, double* A,
                            double* x, double* y, MPI_Comm comm)
{
    // distribute the array onto local processors!
    double* local_A = NULL;
    double* local_x = NULL;
    distribute_matrix(n, &A[0], &local_A, comm);
    distribute_vector(n, &x[0], &local_x, comm);

    // allocate local result space
    double* local_y = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_matrix_vector_mult(n, local_A, local_x, local_y, comm);

    // gather results back to rank 0
    gather_vector(n, local_y, y, comm);
}

// wraps the distributed jacobi function
void mpi_jacobi(const int n, double* A, double* b, double* x, MPI_Comm comm,
                int max_iter, double l2_termination)
{
    // distribute the array onto local processors!
    double* local_A = NULL;
    double* local_b = NULL;
    distribute_matrix(n, &A[0], &local_A, comm);
    distribute_vector(n, &b[0], &local_b, comm);

    // allocate local result space
    double* local_x = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_jacobi(n, local_A, local_b, local_x, comm, max_iter, l2_termination);

    // gather results back to rank 0
    gather_vector(n, local_x, x, comm);
}
