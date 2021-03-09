/*
 * CSE 6220 Introduction to High Performance Computing
 * Programming Assignment 2
 *
 * MPI polynomial evaluation algorithm function definitions
 *
 */

/* 
 * File:   evaluator.h
 */

/*********************************************************************
 *                  !!  DO NOT CHANGE THIS FILE  !!                  *
 *********************************************************************/

#ifndef MPI_EVALUATOR_H
#define MPI_EVALUATOR_H

#include <mpi.h>

/**
 * 
 * Block decomposes an arrays of size n among processors. At the begining of the 
 * algorithm only the source_rank processor knows what elements it has, how many 
 * it has and which of them needs to be scattered. The rest of the processors
 * only knows the source_rank processor has something for them to receive. 
 * 
 * @param n                 Size of the array that needs to be scattered. This has
 *                          a valid value only in the source_rank processor.
 * @param scatter_values    Pointer to an array containing the data that needs 
 *                          scattering. This has a valid value only in the
 *                          source_rank processor.
 * @param n_local           The function will update this with the size of the 
 *                          local portion of the scattered array. Upon return of
 *                          this function each processor will know how much of 
 *                          elements it received from source_rank processor. 
 * @param local_values      Pointer to an array containing the local portion of
 *                          the scattered array. 
 *                          Note: 
 *                              The function implementation should allocate  
 *                              memory to this array before saving any data to 
 *                              this {eg: using malloc(...)}. This memory will be
 *                              free by the calling function, thus you do not
 *                              need to worry about cleaning up. 
 * @param source_rank       Rank of the processor containing the scatter_values.
 * @param comm              MPI communicator object
 */
void scatter(const int n, double* scatter_values, int &n_local, 
        double* &local_values, int source_rank, const MPI_Comm comm);

/**
 * 
 * Broadcast a single value among all other processors. 
 * 
 * @param value             Value to broadcast. Only the source_rank will have 
 *                          a valid value for this variable. Other processors 
 *                          should ignore the value of this variable at the start
 *                          of the algorithm.
 * @param source_rank       Processor containing the value to broadcast
 * @param comm              MPI communicator object
 * @return                  Value obtained through broadcast. At the end of the 
 *                          function each processor should have the value that 
 *                          has being sent via the broadcast algorithm. 
 */
double broadcast(double value, int source_rank, const MPI_Comm comm);

/**
 * 
 * Perform parallel prefix algorithm on an array of n elements with SUM or 
 * PRODUCT operation. This function will perform the relevant local prefix
 * operations as needed. 
 * 
 * eg:  OP             = PREFIX_OP_SUM
 *      values         = x_0, x_1, ..., x_{n-1}
 *      prefix_result  = (T + x_0), (T+ x_0 + x_1), ..., (T + x_0 + ... + x_{n-1})
 *                       {T: sum of all the numbers of "values" arrays of procs < my_rank}
 * 
 * @param n                 Size of the array containing the values. This is 
 *                          the local size of the array "values".
 * @param values            Pointer to the local input array containing values
 *                          which would participate in the parallel prefix
 *                          algorithm. 
 * @param prefix_results    Pointer to the array containing parallel prefix 
 *                          sum result
 * @param OP                PREFIX_OP_SUM | PREFIX_OP_PRODUCT 
 *                          (defined in const.h)
 * @param comm              MPI communicator object
 */
void parallel_prefix(const int n, const double* values, double* prefix_results, const int OP, const MPI_Comm comm);

/**
 * 
 * Perform parallel polynomial evaluation given the value for x and the constants
 * 
 * @param x                 Value of x to be applied to the polynomial function.
 *                          This is available for all processors. 
 * @param n                 No of constants. This is the number of 
 * @param constants         Pointer to an array containing the constants
 * @param comm              MPI communicator object
 * @return                  Result of evaluating the polynomial function
 */
double mpi_poly_evaluator(const double x, const int n, const double* constants, const MPI_Comm comm);

#endif /* MPI_EVALUATOR_H */

