/*
 *  CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 2
 * 
 *  Utility function definitions
 * 
 */

/* 
 * File:   utils.h
 */

#ifndef UTILS_H
#define UTILS_H

#include <mpi.h>

//------------------- Timer Functions (Do not change) -------------------//
/**
 * Update to current time
 * @param t         Time saved in this data structure
 * @param rank      Rank of the processor requesting time
 * @param comm      MPI communication object
 */
void set_time(double &t, const int rank, MPI_Comm comm);

/**
 * Find the time durations and return the the result in seconds
 * @param t_start   Duration start time
 * @param t_end     Duration end time
 * @return          Duration in seconds
 */
double get_duration(double &t_start, double &t_end);
//---------------------------------------------------------------------//

/*********************************************************************
 *                  DECLARE YOUR OWN FUNCTIONS HERE                  *
 *********************************************************************/


// ...

#endif /* UTILS_H */

