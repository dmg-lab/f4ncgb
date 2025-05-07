# Use with gnuplot -e "name='braid3.ms'" index_map.gnuplot

set autoscale
set terminal pngcairo size 3000,2000
set output sprintf("%s_monomial_store_accesses.png", name)

set title "Monomial Store Accesses with operator[]"
set xlabel "Monomial Idx"
set ylabel "Time [ms]"
set pointsize 0.5

plot sprintf("<zstdcat %s_monomial_store.zst", name) using 1:2 w p notitle

set output sprintf("%s_polynomial_store_accesses.png", name)
set title "Polynomial Store Accesses with operator[]"
set xlabel "Polynomial Idx"

plot sprintf("<zstdcat %s_polynomial_store.zst", name) using 1:2 w p notitle
