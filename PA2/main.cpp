/*
 * CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 2
 * 
 *  Main executor for poly-eval
 * 
 */

/* 
 * File:   main.cpp
 */

/*********************************************************************
 *                  !!  DO NOT CHANGE THIS FILE  !!                  *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include "mpi_evaluator.h"
#include "evaluator.h"
#include "const.h"
#include "io.h"
#include "utils.h"

using namespace std;

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    int p, rank;
    const MPI_Comm comm = MPI_COMM_WORLD;
    
    //get the rank and number of processors
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);
    
    // Get the input data
    int n, m;
    vector<double> global_constants;
    vector<double> x;
    if (rank == 0) {    
        setup(argc, argv, n, global_constants, m, x);
    }
    
    double poly_evaltime;
    double t_start, t_end;
    
    if (p>1){
        //------------------- run the parallel version ---------------------//
        int local_n;
        double* local_values;
        scatter(n, &global_constants[0], local_n, local_values, 0, comm);
        MPI_Bcast(&m,1,MPI_INT, 0, comm);   //let all procs know how many evals
        double broadcast_time;
        for(int i=0; i<m; i++){
            double v = (rank == 0)? x[i] : 0;   //x value
            
            set_time(t_start, rank, comm);
            v = broadcast(v, 0, comm);          //send x value to all the procs
            set_time(t_end, rank, comm);
            broadcast_time = get_duration(t_start, t_end);

            set_time(t_start, rank, comm);
            double result = mpi_poly_evaluator(v, local_n, local_values, comm); //eval
            set_time(t_end, rank, comm);
            poly_evaltime = get_duration(t_start, t_end);

            print_output(v, result, broadcast_time, poly_evaltime, rank, p); //print result
        }
        free(local_values);
        //----------------------------------------------------------------//
    } else {
        //------------------- run the serial version ---------------------//
        for(int i=0; i<m; i++){
            set_time(t_start, -1, comm);
            double result = poly_evaluator(x[i], n, &global_constants[0]); //eval
            set_time(t_end, -1, comm);
            poly_evaltime = get_duration(t_start, t_end);
            print_output(x[i], result, -1, poly_evaltime, 0, 0);            //print result
        }
        //----------------------------------------------------------------//   
    }
    MPI_Finalize();
    return (EXIT_SUCCESS);
}


