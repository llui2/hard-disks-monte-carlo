
set size ratio 1

set style line 1 lc rgb 'red' pt 7 ps 1 lw 2

set logscale xy
set format x "10^{%T}"
set format y "10^{%T}"

set title ""
set xlabel "delta" 
set ylabel "D"

plot 'data_2.dat' using 1:2 w linespoints ls 1 notitle

pause -1