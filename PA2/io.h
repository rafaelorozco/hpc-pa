/*
 * CSE 6220 Introduction to High Performance Computing
 *              Programming Assignment 2
 * 
 *  All IO related functions
 * 
 */

/* 
 * File:   io.h
 */

/*********************************************************************
 *                  !!  DO NOT CHANGE THIS FILE  !!                  *
 *********************************************************************/

#ifndef IO_H
#define IO_H

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <sstream>
#include <map>
#include <list>
#include "const.h"


#include <mpi.h>
#include <iostream>
#include <fstream>

/**
 * Help screen
 */
void print_usage(int argc, char** argv){
    int seed = RANDOM_SEED;
    std::cerr << "Usage: "<<std::endl;
    std::cerr << "    "<<argv[0]<<" <file1> <file2>" << std::endl;
    std::cerr << "                  Read the constants from file1 and values from file2" << std::endl;
    std::cerr << "                  and perform polynomial evaluation." << std::endl;
    std::cerr << "    "<<argv[0]<<" -n <n> -m <m> [-s <random-seed>]" << std::endl;
    std::cerr << "                  Creates n number of random constants and m number " << std::endl;
    std::cerr << "                  of x values to perform polynomial evaluation.{n,m>0}" << std::endl;
    std::cerr << "                  (optional, random-seed = "<<seed<<" [default])." << std::endl;
    std::cerr << "Output:"<<std::endl;
    std::cerr << "    <x_1> <y_1> <broadcast-time> <polynimial-eval-time>" << std::endl;
    std::cerr << "    <x_2> <y_2> <broadcast-time> <polynimial-eval-time>" << std::endl;
    std::cerr << "                          ..." << std::endl;
    std::cerr << "    <x_m> <y_m> <broadcast-time> <polynimial-eval-time>" << std::endl;
    std::cerr << "                  (broadcast-time shown only for parallel runs)" << std::endl;    
    std::cerr << std::endl;    
}

std::map<std::string, double> get_settings(){
    std::map<std::string, double> configuration;
    
    //setup default settings
    configuration[XLIMIT] = X_LIMIT;
    configuration[CONSTLIMIT] = CONST_LIMIT;
    configuration[RANDOMSEED] = RANDOM_SEED;
    configuration[OUTPUT_SHOWX] = SHOW_X;
    configuration[OUTPUT_SHOWY] = SHOW_Y;
    configuration[OUTPUT_SHOWBCAST_TIME] = SHOW_BCAST_TIME;
    configuration[OUTPUT_SHOWPOLYEVAL_TIME] = SHOW_POLY_EVAL_TIME;

    //get settings from the config file
    std::ifstream configfile(CONFIG_FILE.c_str());
    if (configfile.good()){
        std::string line;
        while (std::getline(configfile, line)){
            std::stringstream line_stream(line); //eval each line as key-value pair
            std::string key;
            if (std::getline(line_stream, key, '=')) {
                if (key[0] == '#')
                    continue;       //ignore comment lines
                
                std::string value;
                if (std::getline(line_stream, value)){
                    configuration[key] = atof(value.c_str());
                }
            }
        }
        configfile.close();
    }
    
    return configuration;
}

/**
 * Print the output based on the settings
 */
void print_output(double x, double y, double broadcast_time, 
        double poly_evaltime, int rank, int p){
    if (rank <= 0){
        std::map<std::string, double> settings = get_settings();
        bool show_output = settings[OUTPUT_SHOWX] || 
                            settings[OUTPUT_SHOWY] || 
                            (p>1 && settings[OUTPUT_SHOWBCAST_TIME]) || 
                            settings[OUTPUT_SHOWPOLYEVAL_TIME];
        if (show_output){
            std::string space="";
            if (settings[OUTPUT_SHOWX]){
                std::cout<<x;
                space="\t";
            }
            if (settings[OUTPUT_SHOWY]){
                std::cout<<space<<y;
                space="\t";
            }                
            if (p>1 && settings[OUTPUT_SHOWBCAST_TIME]){
                std::cout<<space<<broadcast_time;
                space="\t";
            }                
            if (settings[OUTPUT_SHOWPOLYEVAL_TIME]){
                std::cout<<space<<poly_evaltime;
                space="\t";
            } 
            std::cout<<std::endl;
        }
    }
}

/**
 * Read a file containing data type of T in to a vector. First value in the file
 * will be no of T values in the file.
 * @param filename
 * @param list
 * @return 
 */
template <typename T>
int read_file(std::string filename, std::vector<T> &list){
    //Check the file for error
    std::ifstream in(filename, std::ios::binary | std::ios::ate);
    if (!(in.good() && in.is_open()))
        throw std::runtime_error(std::string("Couldn't open file ") + filename);
    in.close();    
    
    //open the file for reading
    std::fstream datafile(filename, std::ios_base::in);
    int n;
    datafile >> n;
    for(int i=0; i<n; i++){
        T v;
        datafile >> v;
        list.push_back(v);
    }
    datafile.close();
    return n;
}

/**
 * 
 * Get Input Data
 * 
 * @param n                 Total no of polynomial constants
 * @param global_constants  Vector to store the polynomial constants
 * @param m                 Total no of x values to evaluate
 * @param x                 Vector to store the x values
 * @param rank              Rank of the current processor
 */
void setup(int argc, char** argv, int &n, std::vector<double> &global_constants, int &m, std::vector<double> &x){
    std::map<std::string, double> settings = get_settings();
    if (argc < 3){
        print_usage(argc, argv);    //not enough parameters to run the application
        exit(EXIT_FAILURE);
    }

    if (std::string(argv[1]) == "-n"){
        // randomly generate input
        
        if (argc<5){
            print_usage(argc, argv); //not enough parameters to run the application
            exit(EXIT_FAILURE);
        }
        n = atoi(argv[2]);
        if (!(n > 0)) {
            print_usage(argc, argv);    //n should be positive
            exit(EXIT_FAILURE);
        }
        if (std::string(argv[3]) == "-m"){
            m = atoi(argv[4]);
            if (!(m > 0)) {
                print_usage(argc, argv); //m should be positive
                exit(EXIT_FAILURE);
            }
        }
        if (argc==7){           
            if (std::string(argv[5]) == "-s"){ //if seed value provided
                settings[RANDOMSEED] = atoi(argv[6]);
            }            
        }
        
        srand(settings[RANDOMSEED]);

        // generate random constants
        double const_limit = settings[CONSTLIMIT];
        for(int i=0; i<n; i++){
            double randvalue = (1.0 * (2*const_limit) * rand() / RAND_MAX) - const_limit;
            global_constants.push_back(randvalue);
        }

        // generate random x values
        double x_limit = settings[XLIMIT];
        for(int i=0; i<m; i++){
            double randvalue = (1.0 * (2*x_limit) * rand() / RAND_MAX) - x_limit;
            x.push_back(randvalue);
        }
    }else{
        if (argc != 3){
            print_usage(argc, argv);
            exit(EXIT_FAILURE);
        }
        //read data from files
        n = read_file(std::string(argv[1]), global_constants);
        m = read_file(std::string(argv[2]), x);

        //make sure there are constants and x values for the algorithm to eval upon
        if (n<=0 || global_constants.size()==0){
           throw std::runtime_error("There should be at least one constant in the file!!!");
        }
        if (m<=0 || x.size()==0){
           throw std::runtime_error("There should be at least one value for x in the file!!!");
        }
    }
}

#endif /* IO_H */

