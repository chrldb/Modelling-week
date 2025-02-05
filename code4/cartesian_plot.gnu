set terminal pngcairo enhanced size 800,600
set output "cartesian_plot.png"
set xlabel "Angle (°)"
set ylabel "Gain (dB)"
set grid
set title "Radiation pattern in dB"
set xrange [-100:100]
set yrange [-85:0]
plot "data_cartesian.dat" using 1:2 with lines lw 2 lc rgb "blue" title " "
