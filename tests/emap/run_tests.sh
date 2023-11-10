#!/bin/bash
rm -f main
g++ -O3 -g -frounding-math -std=c++17 main.cxx -o main -lm -lgsl -lgslcblas -lcuba -lnlopt
xrange=150
yrange=100

ystart=0
yend=100
#yend=`echo "$ystart+10" | bc`

for (( x=0; x<=$xrange; x++ ))
do
  for (( y=$ystart; y<=$yend; y++ ))
  do
    #echo $x $y
    
    x1=`bc <<< "scale=8; 1.5*$x/$xrange" | awk '{printf "%1.8f", $0}'`
    y1=`bc <<< "scale=8; 0.25*$y/$yrange" | awk '{printf "%1.8f", $0}'`
    
    cat data.run | sed -e "s/xx/$x1/g" | sed -e "s/yy/$y1/g" > data.curr
    
    sleep 1
    out=`./main data.curr 16`
    echo $x1 $y1 $out
    
  done
done
