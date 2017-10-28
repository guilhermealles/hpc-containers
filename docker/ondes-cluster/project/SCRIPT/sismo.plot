set terminal png size 800,600
set output 'obs1.png'

set xlabel "Time"
set ylabel "Velocity"
set title "sismogram obs1.dat"
plot 'output/obs1.dat' using 1:2 with lines title 'Vx','output/obs1.dat' using 1:3 with lines title 'Vy','output/obs1.dat' using 1:3 with lines title 'Vz'
