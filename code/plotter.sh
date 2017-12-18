#!/bin/bash

ensemble_file=$1
phi=100.
n=3.



sed 1,/count/d $ensemble_file  | awk  '$1 !~ /Network/{print $0}' > $$.tmp1
awk '/Network/{$1=""; print $0}' $ensemble_file  > $$.tmp2
paste $$.tmp1 $$.tmp2 > $$.tmp3 
awk '{$1=""; $2="" ; print $0}' $$.tmp3 | sort -nr  | uniq > $$.tmp4

cnt=0
while read line
do
    ((cnt++))
    T=$(echo $line | awk '{print $1}')
    sig=$(echo $line | awk '{print $2}')
    edge_data="$(echo $line | sed "s/^.*\[/\[/")"

    # prepare the network diagram
    ./genera_network.py "$edge_data" network_cnt_${cnt}_T_${T}_sig_${sig}.eps
    
    # prepare the u_0 vs t dynamical data 
    ./genera_dynamics.py "$edge_data"  $phi  $n > $$.dynamics.dat
    outfile="dynamics_cnt_${cnt}_T_${T}_sig_${sig}.eps"
    
    # prepare a gnuplot session for dynamics plot
cat > dynamics.plt << EOF
set terminal postscript eps enhanced color font 'Helvetica,12' 
set output "$outfile"
set size 0.6,0.6
set xlabel "Time" font "Helvetica,18"
set ylabel "Concentration of u_0" font "Helvetica,18"
plot "$$.dynamics.dat" u 1:2 notitle w l
EOF
    gnuplot -p < dynamics.plt
done < $$.tmp4

# clean the directory of  all temporary files
rm $$.tmp1 $$.tmp2 $$.tmp3 $$.tmp4 $$.dynamics.dat dynamics.plt


