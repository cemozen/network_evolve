set terminal postscript eps enhanced color font 'Helvetica,12' 
set output "dynamics_cnt_36_T_16.1200_sig_5.3733.eps"
set size 0.6,0.6
set xlabel "Time" font "Helvetica,18"
set ylabel "Concentration of u_0" font "Helvetica,18"
plot "8248.dynamics.dat" u 1:2 notitle w l
