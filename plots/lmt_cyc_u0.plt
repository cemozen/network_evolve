
#set terminal postscript eps enhanced

set xlabel "t" font "Helvetica,26"
set ylabel "p_{0}" font "Helvetica,26"
set label "{/Symbol f}=9.7115" font "Helvetica,20" at 0.4, 0.35 tc lt 1
#set output "lmt_cyc.eps"
plot [250:] "u0_phi_9.7115.dat"  notitle w l lc 3 

