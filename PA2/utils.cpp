/*
 * CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 2
 * 
 *  Utility function implementation
 * 
 */

#include <mpi.h>
#include "utils.h"

//------------------- Timer Functions (Do not change) -------------------//
void set_time(double &t, const int rank, MPI_Comm comm){
    if (rank>=0) // Do not call barrier if rank is negative
        MPI_Barrier(comm);
    if (rank <= 0){ //only 1 processor will set the time
        t = MPI_Wtime();
    }
}

double get_duration(double &t_start, double &t_end){
    return (t_end - t_start);
}
//---------------------------------------------------------------------//

/*********************************************************************
 *                 Implement your own functions here                 *
 *********************************************************************/




// ...