#!/bin/bash
# Make sure genera.py is called so that network_plot_from_edges()  is active in the main code.
# To run this script, redirect a data file made of lines such as:
# [(0, 3, 'R'), (1, 2, 'R'), (1, 3, 'R'), (2, 0, 'R'), (2, 1, 'R'), (3, 0, 'R'), (3, 1, 'A'), (3, 2, 'R')]

cnt=0
while read edge_data
do
    ((cnt++))
    ./genera.py "$edge_data" network_$(printf "%05d%s\n" $cnt).eps
done
