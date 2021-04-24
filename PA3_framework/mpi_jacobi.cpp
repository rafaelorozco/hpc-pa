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

    int dims[2], periods[2], coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);

    MPI_Comm col_comm;
    int keepdims[2] = {1, 0};
    MPI_Cart_sub(comm, keepdims, &col_comm);

    if (coords[1] == 0) return;
    
    int* count = new int[dims[0]];
    int* displs = new int[dims[0]];
    
    for (int i=0; i<dims[0]; i++) {
        count[i] = block_decompose(n, dims[0], i);
    }

    displs[0] = 0;
    for (int i = 1; i < dims[0]; i++) {
        displs[i] = displs[i-1] + count[i-1];
    }

    int local_size = count[coords[0]];
    (*local_vector) = new double[local_size];

    int root_rank;
    int root_coords[] = {0};
    MPI_Cart_rank(col_comm, root_coords, &root_rank);

    MPI_Scatterv(input_vector, count, displs, MPI_DOUBLE, *local_vector, local_size, MPI_DOUBLE, root_rank, col_comm);

    delete [] count;
    delete [] displs;

    return;
}


// gather the local vector distributed among (i,0) to the processor (0,0)
void gather_vector(const int n, double* local_vector, double* output_vector, MPI_Comm comm)
{
    // TODO

    int dims[2], periods[2], coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);

    MPI_Comm col_comm;
    int keepdims[2] = {1, 0};
    MPI_Cart_sub(comm, keepdims, &col_comm);

    if (coords[1] == 0) return;


    int* count = new int[dims[0]];
    int* displs = new int[dims[0]];

    for (int i=0; i<dims[0]; i++) {
        count[i] = block_decompose(n, dims[0], i);
    }

    displs[0] = 0;
    for (int i = 1; i < dims[0]; i++) {
        displs[i] = displs[i - 1] + count[i - 1];
    }

    int local_size = count[coords[0]];

    int root_rank;
    int root_coords[] = {0};
    MPI_Cart_rank(col_comm, root_coords, &root_rank);

    MPI_Gatherv(local_vector, local_size, MPI_DOUBLE, output_vector, count, displs, MPI_DOUBLE, root_rank, col_comm);

    delete [] count;
    delete [] displs;

    return;
}


void distribute_matrix(const int n, double* input_matrix, double** local_matrix, MPI_Comm comm)
{
    // TODO

    int dims[2], periods[2], coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);

    // Step 1. Distribute rows on processors in the first column
    // Create column communicators 
    MPI_Comm col_comm;
    int keepdims[2] = {1, 0};
    MPI_Cart_sub(comm, keepdims, &col_comm);

    // Create a temp matrix to store the values in the first column
    double* temp = NULL;

    if (coords[1] == 0) return;

    // Compute the number of elements to send to each processor
    int* count = new int[dims[0]];
    int* displs = new int[dims[0]];

    for (int i = 0; i < dims[0]; i++) {
        count[i] = n * block_decompose(n, dims[0], i);
    }

    displs[0] = 0;
    for (int i = 1; i < dims[0]; i++) {
        displs[i] = displs[i - 1] + count[i - 1];
    }

    int local_size = count[coords[0]];
    temp = new double[local_size];

    // Get the rank of root in the first column communicator
    int root_rank;
    int root_coords[] = {0};
    MPI_Cart_rank(col_comm, root_coords, &root_rank);

    // Scatter values to different processors
    MPI_Scatterv(input_matrix, count, displs, MPI_DOUBLE, temp, local_size, MPI_DOUBLE, root_rank, col_comm);

    // Step 2. Distribute submatrix on processors among each row
    // Create row communicators 
    MPI_Comm row_comm;
    keepdims[0] = 0; keepdims[1] = 1;
    MPI_Cart_sub(comm, keepdims, &row_comm);

    int num_row = block_decompose(n, dims[0], coords[0]);
    int num_col = block_decompose(n, dims[1], coords[1]);

    // Compute the number of elements to send to each processor

    for (int i = 0; i < dims[1]; i++) {
        count[i] = block_decompose(n, dims[1], i);
    }

    displs[0] = 0;
    for (int i = 1; i < dims[1]; i++) {
        displs[i] = displs[i - 1] + count[i - 1];
    }

    // Allocate spaces for local matrix
    (*local_matrix) = new double[num_row * num_col];

    // Get the rank of root in the first column
    MPI_Cart_rank(row_comm, root_coords, &root_rank);

    for (int i = 0; i < num_row; i++) {
        MPI_Scatterv((temp + i * n), count, displs, MPI_DOUBLE, (*local_matrix + i * num_col), num_col, MPI_DOUBLE, root_rank, row_comm);
    }

    delete [] count;
    delete [] displs;
    delete [] temp;

    return;
}


void transpose_bcast_vector(const int n, double* col_vector, double* row_vector, MPI_Comm comm)
{
    // TODO
    
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    //calculate the number of processors in each row and column in the grid
    int p, q;
    MPI_Comm_size(comm, &p);
    q = (int)sqrt(p);

    //get row and column subcommunicators
    MPI_Comm col_comm, row_comm;
    int keepdims[2] = {true, false};
    MPI_Cart_sub(comm, keepdims, &col_comm);
    keepdims[0] = false; keepdims[1] = true;
    MPI_Cart_sub(comm, keepdims, &row_comm);

    //get rank in the column and in the row
    int col_rank, row_rank;
    int coords[2];
    MPI_Cart_coords(comm, rank, 2, coords);
    row_rank = coords[0];
    col_rank = coords[1];

    if(rank == 0) {
        int count = block_decompose(n, q, rank);
        memcpy (row_vector, col_vector, count*sizeof(double));
    } else if(col_rank == 0) {
        int scount = block_decompose(n, q, row_rank);
        MPI_Send(col_vector, scount, MPI_DOUBLE, row_rank, 0, row_comm);
    } else if(row_rank == col_rank) {
        int rcount = block_decompose(n, q, row_rank);
        MPI_Recv(row_vector, rcount, MPI_DOUBLE, 0, 0, row_comm, MPI_STATUS_IGNORE);
    }

    //broadcast the vector within the column
    int bcount = block_decompose(n, q, col_rank);
    MPI_Bcast(row_vector, bcount, MPI_DOUBLE, col_rank, col_comm);

}


void distributed_matrix_vector_mult(const int n, double* local_A, double* local_x, double* local_y, MPI_Comm comm)
{
    int p, rank;
    int crank;
    int my_row_rank;
    int my_col_rank;

    int coords[NDIM];
    int free_coords[NDIM];

    MPI_Comm comm1D_row;
    MPI_Comm comm1D_col;

    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    MPI_Comm_rank(comm, &crank);
    MPI_Cart_coords(comm, crank, NDIM, coords);

   
    //How many elements from x are here?
    int msg_size;
    int q = sqrt(p);
    //printf("q %d",q);
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
    local_x[0] = 1;
    //printf("element %f",local_x[0]);
    fflush(stdout);
    //printf("I'm %d %d and my row rank is %d and my msgsize %d and element \n", coords[0], coords[1], my_row_rank,msg_size);

    // //////Need to tranpose vector onto grid 
    //1.) (i,0) processors send their local element send to diagonal element. (i,0)
    MPI_Status stat;
    if(coords[0] != 0) {
        if(my_row_rank == 0) {
            //printf("I'm %d %d", coords[0], coords[1]);
            //printf("\n ");
            //local_x[0] = 1;
            printf("and I am sending to my element %f com rank %d\n",local_x[0],coords[0]);
            fflush(stdout);

            MPI_Send(local_x, msg_size, MPI_DOUBLE, coords[0], 11, comm1D_row);
        }
        if(my_row_rank == coords[0]) {
            printf("I'm waiting and I am %d %d and my row rank is %d  \n", coords[0], coords[1], my_row_rank);

    
            fflush(stdout);
            MPI_Recv(local_x, msg_size, MPI_DOUBLE, 0, 11, comm1D_row, &stat);
            printf("I got element %f\n", local_x[0]);
            //2.) (i,i) processors send their local element to all others in col communicator. 
        }
    }
    if(coords[1] == coords[0]) {
        MPI_Bcast(&local_x, msg_size, MPI_DOUBLE, coords[0], comm1D_col);
    }
    //wait
    MPI_Barrier(comm);

    //printf("I'm %d %d and my row rank is %d and my msgsize %d and element %f \n", coords[0], coords[1], my_row_rank,msg_size, local_x);


    //Do local multiplication using sequential matvec???
    matrix_vector_mult(n, &A[0], &x[0], &y[0]);

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
