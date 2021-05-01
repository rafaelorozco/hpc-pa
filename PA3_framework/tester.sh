#!/bin/bash

for i in {0.0,0.5,1.0}
do
	echo "Run with  d=$i"
	for j in {1000..21000..5000}
	do 
	mpirun -np 16 ./jacobi -n $j -d $i
#
	done
done

#for i in {0.0,0.5,1.0}
#do
#	echo "Run with  d=$i"
#	for j in {1000..21000..5000}
#	do 
#  		./jacobi -n $j -d $i
#	done
#done


