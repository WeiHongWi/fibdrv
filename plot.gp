set title "perf"
set xlabel "nth fibonacci"
set ylabel "second(ns)"
set xtics 0,100,2000
set ytics 0,1000,18000
set key left
set term png enhanced font 'Verdana,10'
set output 'runtime.png'
plot "datav1.txt" using 1:2 with linespoints linewidth 2 title "fast doubling v4","datav2.txt" using 1:2 with linespoints linewidth 2 title "fast doubling v3"
