set key outside
set title ""
set xlabel "t" 
set ylabel "MSD"
 
set style line 1 lc rgb '#D93030' pt 0 ps 1 lw 1
set style line 2 lc rgb '#CA612D' pt 0 ps 1 lw 1
set style line 3 lc rgb '#BB8B2B' pt 0 ps 1 lw 1
set style line 4 lc rgb '#ACAD28' pt 0 ps 1 lw 1
set style line 5 lc rgb '#769E25' pt 0 ps 1 lw 1
set style line 6 lc rgb '#479023' pt 0 ps 1 lw 1

set logscale xy
set format x "10^{%T}"
set format y "10^{%T}"

delta = '0.001 0.003 0.01 0.03 0.1 0.3'

plot for [n=0:5] "data_1.dat" i n u 1:2 w linespoints ls n+1 title sprintf("delta = ".word(delta,n+1))
pause -1