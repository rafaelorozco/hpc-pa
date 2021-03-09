/*
 *  CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 2
 * 
 *  All constants defined in the application reside here
 * 
 */

/* 
 * File:   const.h
 */

/*********************************************************************
 *                  !!  DO NOT CHANGE THIS FILE  !!                  *
 *********************************************************************/

#ifndef CONST_H
#define CONST_H

#include <string>

//Parallel Prefix Operators
#define PREFIX_OP_SUM 1         //"+" operator 
#define PREFIX_OP_PRODUCT 2     //"*" operator 

//Random input generation default settings
#define X_LIMIT 2           //-X_LIMIT <= x <= X_LIMIT
#define CONST_LIMIT 100000  //-CONST_LIMIT <= constants <= CONST_LIMIT
#define RANDOM_SEED 0;      //Random seed value     

//Output generation settings {1 - show, 0 - hide}
#define SHOW_X 1
#define SHOW_Y 1
#define SHOW_BCAST_TIME 1
#define SHOW_POLY_EVAL_TIME 1


//Config file property names
const std::string XLIMIT="x-limit";
const std::string CONSTLIMIT="contant-limit";
const std::string RANDOMSEED="seed";
const std::string OUTPUT_SHOWX="show-x";
const std::string OUTPUT_SHOWY="show-y";
const std::string OUTPUT_SHOWBCAST_TIME="show-bcast-time";
const std::string OUTPUT_SHOWPOLYEVAL_TIME="show-poly-eval-time";

const std::string CONFIG_FILE="config.settings";


#endif /* CONST_H */

