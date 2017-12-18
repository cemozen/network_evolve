#!/bin/bash
# Make sure genera.py is called so that dynamix_plot_from_edges()  is active in the main code.
# To run this script, redirect a data file made of lines such as:
# [(0, 3, 'R'), (1, 2, 'R'), (1, 3, 'R'), (2, 0, 'R'), (2, 1, 'R'), (3, 0, 'R'), (3, 1, 'A'), (3, 2, 'R')]

# Make sure the phi, n parameters are consistent with the intended values as may have used before to generate 
# the graphs that is the input here.
phi=100.
n=3.
cnt=0
while read edge_data
do
((cnt++))
# first we output in block numbers the dynamics of all nodes
./genera.py "$edge_data"  $phi  $n > $$.dynamics.dat  
outfile="dynamics_$(printf "%05d%s\n" $cnt).eps"
# prepare a gnuplot script here
cat > dynamics.plt << EOF
set terminal postscript eps enhanced color font 'Helvetica,10'
set output "$outfile"
plot "$$.dynamics.dat" u 1:2 w l
EOF
gnuplot -p < dynamics.plt
rm $$.dynamics.dat
done
