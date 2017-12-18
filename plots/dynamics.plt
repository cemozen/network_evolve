set terminal postscript eps enhanced color font 'Helvetica,10'
set output "dynamics_00001.eps"
plot "3523.dynamics.dat" u 1:2 w l
