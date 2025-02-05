set terminal pngcairo enhanced size 800,600
set output "radiation_pattern.png"
set xlabel "Theta (radians)"
set ylabel "Power (dB)"
set title "Radiation pattern in dB"
set grid
plot "radiation_pattern.dat" using 1:2 with lines lw 2 title "Curve 1", "radiation_pattern.dat" using 1:3 with lines lw 2 title "Curve 2", "radiation_pattern.dat" using 1:4 with lines lw 2 title "Curve 3", "radiation_pattern.dat" using 1:5 with lines lw 2 title "Curve 4", "radiation_pattern.dat" using 1:6 with lines lw 2 title "Curve 5", "radiation_pattern.dat" using 1:7 with lines lw 2 title "Curve 6", "radiation_pattern.dat" using 1:8 with lines lw 2 title "Curve 7", "radiation_pattern.dat" using 1:9 with lines lw 2 title "Curve 8", "radiation_pattern.dat" using 1:10 with lines lw 2 title "Curve 9", "radiation_pattern.dat" using 1:11 with lines lw 2 title "Curve 10", "radiation_pattern.dat" using 1:12 with lines lw 2 title "Curve 11", "radiation_pattern.dat" using 1:13 with lines lw 2 title "Curve 12", "radiation_pattern.dat" using 1:14 with lines lw 2 title "Curve 13", "radiation_pattern.dat" using 1:15 with lines lw 2 title "Curve 14", "radiation_pattern.dat" using 1:16 with lines lw 2 title "Curve 15", "radiation_pattern.dat" using 1:17 with lines lw 2 title "Curve 16", "radiation_pattern.dat" using 1:18 with lines lw 2 title "Curve 17", "radiation_pattern.dat" using 1:19 with lines lw 2 title "Curve 18", "radiation_pattern.dat" using 1:20 with lines lw 2 title "Curve 19", "radiation_pattern.dat" using 1:21 with lines lw 2 title "Curve 20", "radiation_pattern.dat" using 1:22 with lines lw 2 title "Curve 21", "radiation_pattern.dat" using 1:23 with lines lw 2 title "Curve 22", "radiation_pattern.dat" using 1:24 with lines lw 2 title "Curve 23", "radiation_pattern.dat" using 1:25 with lines lw 2 title "Curve 24", "radiation_pattern.dat" using 1:26 with lines lw 2 title "Curve 25", "radiation_pattern.dat" using 1:27 with lines lw 2 title "Curve 26", "radiation_pattern.dat" using 1:28 with lines lw 2 title "Curve 27"
