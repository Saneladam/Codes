
### circulo.gnuplot
reset

# Salida vectorial para publicación
set terminal pdfcairo enhanced color \
    font "Liberation Serif,12" linewidth 1.2 \
    size 12cm,12cm
set output "circulo.pdf"

# Ventana y escala física
set xrange [-1.3:1.3]
set yrange [-1.3:1.3]
set size ratio -1

# Tipografía y ejes
set xlabel "x" font "Liberation Serif,12"
set ylabel "y" font "Liberation Serif,12"
set xtics 0.5
set ytics 0.5
set tics out
set border linewidth 1.0
set zeroaxis linewidth 0.8 dashtype 2

# Estilo de la curva
set style line 1 lc rgb "#1F4E79" lw 1.8

# Círculo parametrizado
set parametric
set trange [0:2*pi]
set samples 400

plot cos(t), sin(t) with lines ls 1 \
    notitle

unset output
