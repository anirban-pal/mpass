#!/bin/bash
rm -f main
g++ -O3 -g -frounding-math -std=c++17 main.cxx -o main -lm -lgsl -lgslcblas -lcuba -lnlopt
for (( c=1; c<=50; c++ ))
do
  for i in {1..10}
  do
    ./main data.run $c
  done
done
