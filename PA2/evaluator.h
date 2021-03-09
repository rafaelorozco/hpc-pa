/*
 * CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 2
 * 
 *  Serial polynomial evaluation algorithm function definitions
 * 
 */

/* 
 * File:   evaluator.h
 */

/*********************************************************************
 *                  !!  DO NOT CHANGE THIS FILE  !!                  *
 *********************************************************************/

#ifndef EVALUATOR_H
#define EVALUATOR_H

/**
 * 
 * Perform polynomial evaluation given the value for x and the constants
 * 
 * @param x                 Value of x to be applied to the polynomial function
 * @param n                 No of constants
 * @param constants         Pointer to an array containing the constants
 * @param comm              MPI communicator object
 * @return                  Result of evaluating the polynomial function
 */
double poly_evaluator(const double x, const int n, const double* constants);

#endif /* EVALUATOR_H */

