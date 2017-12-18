
#set terminal postscript eps enhanced

set xlabel "p_{0}" font "Helvetica,26"
set ylabel "p_{1}" font "Helvetica,26"
set zlabel "p_{2}" font "Helvetica,26"
set label "{/Symbol f}=9.7115" font "Helvetica,20" at 0.4, 0.35 tc lt 1
set view 55,139
#set output "lmt_cyc.eps"
splot "u_phi_9.7115.dat" u 1:2:3 notitle w l lc 3 

