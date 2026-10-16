# onda.gnuplot
reset

# Formato de salida: PDF vectorial
set terminal pdfcairo enhanced color \
    size 12cm, 8cm font "Arial,11"
set output "onda.pdf"

# Aspecto del gráfico
set title "Onda amortiguada"
set xlabel "Posición x"
set ylabel "Amplitud"

set xrange [0:10]
set yrange [-1.2:1.2]
set border linewidth 1.2
set tics out
set grid back lc rgb "#dddddd" lw 0.7

# Funciones
f(x) = sin(2*pi*x) * exp(-0.12*x)
g(x) = sin(2*pi*x)

# Curvas
set style line 1 lc rgb "#0072B2" lw 2.2
set style line 2 lc rgb "#D55E00" lw 1.5 dt 2

plot f(x) with lines ls 1 title "Amortiguada", \
     g(x) with lines ls 2 title "Sin amortiguamiento"

unset output
