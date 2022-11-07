
set size ratio 1

set style line 1 lc rgb 'red' pt 0 ps 1 lw 1

set title ""
set xlabel "IMC" 
set ylabel "MSD"

set logscale xy
set format x "10^{%T}"
set format y "10^{%T}"

plot 'data_3.dat' using 1:2 w linespoints ls 1 notitle

pause -1