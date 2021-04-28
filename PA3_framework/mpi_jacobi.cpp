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
#include <iostream>
#include <algorithm>

/*
 * TODO: Implement your solutions here
 */


void distribute_vector(const int n, double* input_vector, double** local_vector, MPI_Comm comm) {

    // Get the Cartesian topology information
    int dims[2], periods[2], coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);

    // Create column communicators 
    MPI_Comm comm_col;
    int remain_dims[2] = {true, false};
    MPI_Cart_sub(comm, remain_dims, &comm_col);

    if (coords[1] == 0) {

        // std::cout << coords[0] << " " << coords[1] << std::endl;
        // if (coords[0] == 0 && coords[1] == 0) {
        //     for (int i = 0; i < n; i++) {
        //         std::cout << input_vector[i] << " ";
        //     }
        //     std::cout << std::endl;
        // }

        // Compute the number of elements to send to each processor
        int* sendcounts = new int[dims[0]];
        int* displs = new int[dims[0]];

        for (int i = 0; i < dims[0]; i++) {
            sendcounts[i] = block_decompose(n, dims[0], i);
        }

        displs[0] = 0;
        for (int i = 1; i < dims[0]; i++) {
            displs[i] = displs[i - 1] + sendcounts[i - 1];
        }

        // Compute local size for processors in the first column
        int local_size = sendcounts[coords[0]];
        (*local_vector) = new double[local_size];

        // Get the rank of root in the first column communicator
        int root_rank;
        int root_coords[] = {0};
        MPI_Cart_rank(comm_col, root_coords, &root_rank);

        // std::cout << root_rank << std::endl;
        // for (int i = 0; i < dims[0]; i++) {
        //     std::cout << sendcounts[i] << " ";
        // }
        // std::cout << std::endl;
        // for (int i = 0; i < dims[0]; i++) {
        //     std::cout << displs[i] << " ";
        // }
        // std::cout << std::endl;

        // Scatter values to different processors
        MPI_Scatterv(input_vector, sendcounts, displs, MPI_DOUBLE, *local_vector, local_size, MPI_DOUBLE, root_rank, comm_col);

        // for (int i = 0; i < local_size; i++) {
        //     std::cout << (*local_vector)[i] << " ";
        // }
        // std::cout << std::endl;

        delete [] sendcounts;
        delete [] displs;

    }

    // // Free the column communicator
    // MPI_Comm_free(&comm_col);

    return;
}


// gather the local vector distributed among (i,0) to the processor (0,0)
void gather_vector(const int n, double* local_vector, double* output_vector, MPI_Comm comm) {

    // Get the Cartesian topology information
    int dims[2], periods[2], coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);

    // Create column communicators 
    MPI_Comm comm_col;
    int remain_dims[2] = {true, false};
    MPI_Cart_sub(comm, remain_dims, &comm_col);

    if (coords[1] == 0) {

        // Compute the number of elements to send to each processor
        int* recvcounts = new int[dims[0]];
        int* displs = new int[dims[0]];

        for (int i = 0; i < dims[0]; i++) {
            recvcounts[i] = block_decompose(n, dims[0], i);
        }

        displs[0] = 0;
        for (int i = 1; i < dims[0]; i++) {
            displs[i] = displs[i - 1] + recvcounts[i - 1];
        }

        // Compute local size for processors in the first column
        int local_size = recvcounts[coords[0]];

        // Get the rank of root in the first column communicator
        int root_rank;
        int root_coords[] = {0};
        MPI_Cart_rank(comm_col, root_coords, &root_rank);

        // Gather values from different processors
        MPI_Gatherv(local_vector, local_size, MPI_DOUBLE, output_vector, recvcounts, displs, MPI_DOUBLE, root_rank, comm_col);

        delete [] recvcounts;
        delete [] displs;

    }

    // // Free the column communicator
    // MPI_Comm_free(&comm_col);

    return;
}

void distribute_matrix(const int n, double* input_matrix, double** local_matrix, MPI_Comm comm) {

    // Get the Cartesian topology information
    int dims[2], periods[2], coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);

    // if (coords[0] == 0 && coords[1] == 0) {
    //     std::cout << "===========================" << std::endl;
    //     std::cout << "Before distribution. A " << n << " * " << n << " matrix" << std::endl;
    //     std::cout << "Coordinates (" << coords[0] << " " << coords[1] << "): " << std::endl;;
    //     for (int i = 0; i < n; i++) {
    //         for (int j = 0; j < n; j++) {
    //             std::cout << *(input_matrix + i * n + j) << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << "===========================" << std::endl;
    // }

    // Step 1. Distribute rows on processors in the first column
    // Create column communicators 
    MPI_Comm comm_col;
    int remain_dims[2] = {true, false};
    MPI_Cart_sub(comm, remain_dims, &comm_col);

    // Create a temp matrix to store the values in the first column
    double* temp_matrix = NULL;

    if (coords[1] == 0) {

        // Compute the number of elements to send to each processor
        int* sendcounts = new int[dims[0]];
        int* displs = new int[dims[0]];

        for (int i = 0; i < dims[0]; i++) {
            sendcounts[i] = n * block_decompose(n, dims[0], i);
        }

        displs[0] = 0;
        for (int i = 1; i < dims[0]; i++) {
            displs[i] = displs[i - 1] + sendcounts[i - 1];
        }

        // Compute local size for processors in the first column
        int local_size = sendcounts[coords[0]];
        temp_matrix = new double[local_size];

        // Get the rank of root in the first column communicator
        int root_rank;
        int root_coords[] = {0};
        MPI_Cart_rank(comm_col, root_coords, &root_rank);

        // Scatter values to different processors
        MPI_Scatterv(input_matrix, sendcounts, displs, MPI_DOUBLE, temp_matrix, local_size, MPI_DOUBLE, root_rank, comm_col);

        delete [] sendcounts;
        delete [] displs;

    }

    // // Free the column communicator
    // MPI_Comm_free(&comm_col);

    // Step 2. Distribute submatrix on processors among each row
    // Create row communicators 
    MPI_Comm comm_row;
    remain_dims[0] = false; remain_dims[1] = true;
    MPI_Cart_sub(comm, remain_dims, &comm_row);

    int nrows = block_decompose(n, dims[0], coords[0]);
    int ncols = block_decompose(n, dims[1], coords[1]);

    // Compute the number of elements to send to each processor
    int* sendcounts = new int[dims[1]];
    int* displs = new int[dims[1]];

    for (int i = 0; i < dims[1]; i++) {
        sendcounts[i] = block_decompose(n, dims[1], i);
    }

    displs[0] = 0;
    for (int i = 1; i < dims[1]; i++) {
        displs[i] = displs[i - 1] + sendcounts[i - 1];
    }

    // Allocate spaces for local matrix
    (*local_matrix) = new double[nrows * ncols];

    // Get the rank of root in the first column
    int root_rank;
    int root_coords[] = {0};
    MPI_Cart_rank(comm_row, root_coords, &root_rank);

    for (int i = 0; i < nrows; i++) {
        MPI_Scatterv((temp_matrix + i * n), sendcounts, displs, MPI_DOUBLE, (*local_matrix + i * ncols), ncols, MPI_DOUBLE, root_rank, comm_row);
    }

    delete [] sendcounts;
    delete [] displs;
    delete [] temp_matrix;

    // std::cout << "===========================" << std::endl;

    // std::cout << "Coordinates (" << coords[0] << " " << coords[1] << "): " << std::endl;;
    // for (int i = 0; i < nrows; i++) {
    //     for (int j = 0; j < ncols; j++) {
    //         std::cout << *(*local_matrix + i * ncols + j) << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "===========================" << std::endl;

    // // Free the column communicator
    // MPI_Comm_free(&comm_row);

    return;
}


void transpose_bcast_vector(const int n, double* col_vector, double* row_vector, MPI_Comm comm) {

    //get the rank original rank
    int rank;
    MPI_Comm_rank(comm, &rank);
    
    //calculate the number of processors in each row and column in the grid
    int p, q;
    MPI_Comm_size(comm, &p);
    q = (int)sqrt(p);

    //get row and column subcommunicators
    MPI_Comm comm_row, comm_col;
    int remain_dims[2] = {true, false};
    MPI_Cart_sub(comm, remain_dims, &comm_col);
    remain_dims[0] = false; remain_dims[1] = true;
    MPI_Cart_sub(comm, remain_dims, &comm_row);

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
        MPI_Send(col_vector, scount, MPI_DOUBLE, row_rank, 0, comm_row);
    } else if(row_rank == col_rank) {
        int rcount = block_decompose(n, q, row_rank);
        MPI_Recv(row_vector, rcount, MPI_DOUBLE, 0, 0, comm_row, MPI_STATUS_IGNORE);
    }

    //broadcast the vector within the column
    int bcount = block_decompose(n, q, col_rank);
    MPI_Bcast(row_vector, bcount, MPI_DOUBLE, col_rank, comm_col);

}


void distributed_matrix_vector_mult(const int n, double* local_A, double* local_x, double* local_y, MPI_Comm comm) {

    // Get the Cartesian topology information
    int dims[2], periods[2], coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);

    // Compute the dimentions of local matrix
    int nrows = block_decompose(n, dims[0], coords[0]);
    int ncols = block_decompose(n, dims[1], coords[1]);

    // Get the local vector for matrix multiplication
    double* local_xx = new double[ncols];
    transpose_bcast_vector(n, local_x, local_xx, comm);

    // Get the local result after multiplication
    double* local_yy = new double[nrows];

    if (nrows == ncols) {
        matrix_vector_mult(nrows, local_A, local_xx, local_yy);
    } else {
        matrix_vector_mult(nrows, ncols, local_A, local_xx, local_yy);
    }

    // MPI_Reduce to store the result in the first column
    MPI_Comm comm_row;
    int remain_dims[2] = {false, true};
    MPI_Cart_sub(comm, remain_dims, &comm_row);

    // Get the rank of root in the first column
    int root_rank;
    int root_coords[2] = {0, 0};
    MPI_Cart_rank(comm_row, root_coords, &root_rank);

    MPI_Reduce(local_yy, local_y, nrows, MPI_DOUBLE, MPI_SUM, root_rank, comm_row);

    // // Free the column communicator
    // MPI_Comm_free(&comm_row);

    delete [] local_xx;
    delete [] local_yy;

    return;

}

// Solves Ax = b using the iterative jacobi method
void distributed_jacobi(const int n, double* local_A, double* local_b, double* local_x,
                MPI_Comm comm, int max_iter, double l2_termination) {

    // Get the Cartesian topology information
    int dims[2], periods[2], coords[2];
    MPI_Cart_get(comm, 2, dims, periods, coords);

    // Compute the dimentions of local matrix
    int nrows = block_decompose(n, dims[0], coords[0]);
    int ncols = block_decompose(n, dims[1], coords[1]);

    // Create row and column communicators
    MPI_Comm comm_col;
    int remain_dims[2] = {true, false};
    MPI_Cart_sub(comm, remain_dims, &comm_col);

    MPI_Comm comm_row;
    remain_dims[0] = false; remain_dims[1] = true;
    MPI_Cart_sub(comm, remain_dims, &comm_row);

    // Compute the diagonal elements and matrix R = A - D
    double* local_D = new double[nrows];
    double* local_R = new double[nrows * ncols];
    if (coords[0] == coords[1]) {
        std::copy(local_A, local_A + nrows * ncols, local_R);
        for (int i = 0; i < nrows; i++) {
            local_D[i] = local_A[i * (ncols + 1)];
            local_R[i * (ncols + 1)] = 0;
        }
        if (coords[0] != 0) {
            int dest_rank;
            int dest_coords[] = {0};
            MPI_Cart_rank(comm_row, dest_coords, &dest_rank);
            MPI_Send(local_D, nrows, MPI_DOUBLE, dest_rank, 0, comm_row);
        }
    } else {
        std::copy(local_A, local_A + nrows * ncols, local_R);
    }

    if (coords[1] == 0 && coords[0] != 0) {
        int source_rank;
        int source_coords[] = {coords[0]};
        MPI_Cart_rank(comm_row, source_coords, &source_rank);
        MPI_Recv(local_D, nrows, MPI_DOUBLE, source_rank, 0, comm_row, MPI_STATUS_IGNORE);
    }

    // Initialize x to zero
    if (coords[1] == 0) {
        for (int i = 0; i < nrows; i++) local_x[i] = 0;
    }

    double l2_norm_total;

    int iter = 1;
    while (iter <= max_iter) {

        double* local_y = new double[nrows];

        distributed_matrix_vector_mult(n, local_R, local_x, local_y, comm);

        if (coords[1] == 0) {
            for (int i = 0; i < nrows; i++) {
                local_x[i] = (local_b[i] - local_y[i]) / local_D[i];
            }
        }

        distributed_matrix_vector_mult(n, local_A, local_x, local_y, comm);

        if (coords[1] == 0) {
            double l2_norm = 0.0;
            for (int i = 0; i < nrows; i++) {
                l2_norm += (local_b[i] - local_y[i]) * (local_b[i] - local_y[i]);
            }
            MPI_Allreduce(&l2_norm, &l2_norm_total, 1, MPI_DOUBLE, MPI_SUM, comm_col);
        }

        delete [] local_y;

        l2_norm_total = sqrt(l2_norm_total);

        int root_rank;
        int root_coords[] = {0};
        MPI_Cart_rank(comm_row, root_coords, &root_rank);
        MPI_Bcast(&l2_norm_total, 1, MPI_DOUBLE, root_rank, comm_row);

        if (l2_norm_total <= l2_termination) {
            return;
        }

        iter++;
    }

    // // Free the column communicator
    // MPI_Comm_free(&comm_col);
    // MPI_Comm_free(&comm_row);

    delete [] local_D;
    delete [] local_R;
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
