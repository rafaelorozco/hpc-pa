/**
 * @file    jacobi.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements matrix vector multiplication and Jacobi's method.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */
#include "jacobi.h"

/*
 * TODO: Implement your solutions here
 */

// my implementation:
#include <iostream>
#include <math.h>

// Calculates y = A*x for a square n-by-n matrix A, and n-dimensional vectors x
// and y
void matrix_vector_mult(const int n, const double* A, const double* x, double* y)
{
    for(int i = 0; i < n; i++) {
    	y[i] = 0;
    	for(int j = 0; j < n; j++) {
    		//printf("%d %d\n",i,j);
    		y[i] += A[i*n + j]*x[j];
    	}
    }
}

// Calculates y = A*x for a n-by-m matrix A, a m-dimensional vector x
// and a n-dimensional vector y
void matrix_vector_mult(const int n, const int m, const double* A, const double* x, double* y)
{
    for(int i = 0; i < n; i++) {
    	y[i] = 0;
    	for(int j = 0; j < m; j++) {
    		y[i] += A[i*m + j]*x[j];
    	}
    }
}

double vector_l2(const int n, double* y)
{
	double sum = 0;
    for(int i = 0; i < n; i++) {
    	sum += pow(y[i], 2.);
    }
    return sum;
}

// implements the sequential jacobi method
void jacobi(const int n, double* A, double* b, double* x, int max_iter, double l2_termination)
{
	//make R and inv(D)
	double R[n*n];
	//double R[n*n] = {0};
	//double D_inv[n*n] = {0};
	double D_inv[n*n];
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			if(i != j){
				R[i*n + j] = A[i*n + j];
			} else {
				D_inv[i*n + j] = 1. / A[i*n + j];
			}
		}
    }

    for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			printf("%f ",D_inv[i*n + j]);
			
		}
		printf("\n ");
    }


	//initialize x
    //x = {0};
    //double temp[n] = {0};
    //double res[n] = {0};
	double temp[n];
	double res[n];

    for(int i = 0; i < n; i++) {
    	x[i] = 0;
   	}

    //for(int i = 0; i < n; i++) {
    //	printf("%f ", x[i]);
   	//}
	printf("\n");
   	//for(int i = 0; i < n; i++) {
    //	printf("%f ", temp[i]);
   	//}




    //matrix_vector_mult(n, &D_inv[0], &temp[0], &x[0]);
    double res_norm;


    for(int j = 0; j < n; j++) res[j] = -b[j];

    printf("iniital resnorm %f",vector_l2(n, res));
printf("ls teerm %f",l2_termination);

    while (vector_l2(n, res) > l2_termination) {
    //for(int i = 0; i < 50; i++) {
    	//temp = Rx
    	matrix_vector_mult(n, &R[0], &x[0], &temp[0]);

    	//temp = b - temp
    	// for(int i = 0; i < n; i++) {
    	// 	printf("%f ", temp[i]);
   		// }
   		// printf("\n");
   		// for(int i = 0; i < n; i++) {
    	// 	printf("%f ", b[i]);
   		// }
   		// printf("\n");
    	for(int j = 0; j < n; j++) temp[j] = b[j] - temp[j];
    	// for(int i = 0; i < n; i++) {
    	// 	printf("%f ", temp[i]);
   		// }
   		// printf("\n");

    	//x = inv(D)*temp

    	matrix_vector_mult(n, &D_inv[0], &temp[0], &x[0]);
    	// printf("x\n");
    	// for(int i = 0; i < n; i++) {
    	// 	printf("%f ", x[i]);
   		// }
   		// printf("\n");

   		//check residual
   		matrix_vector_mult(n, &A[0], &x[0], &temp[0]);

   		for(int j = 0; j < n; j++) res[j] = b[j] - temp[j];
   		res_norm = vector_l2(n, res);
   		// printf("resnorm: %f\n",res_norm);
    }

}
