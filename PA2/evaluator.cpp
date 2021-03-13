/*
 *  CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 2
 * 
 *  Serial polynomial evaluation algorithm function implementations goes here
 * 
 */

double poly_evaluator(const double x, const int n, const double* constants){
    //Implementation

	double curr_result = constants[n-1];
	for(int i = n; i > 0; i--) {
		curr_result = constants[i-1] + x*curr_result;
	}
    return curr_result;
}
