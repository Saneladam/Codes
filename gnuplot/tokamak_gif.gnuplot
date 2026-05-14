set parametric
set hidden3d
set view equal xyz
set isosamples 40,30

set pm3d depthorder

unset key
unset grid
unset border
unset tics

set xrange [-8.8:8.8]
set yrange [-8.8:8.8]
set zrange [-3:3]

set style line 1 lc rgb "red" lw 1.5

set term gif animate delay 4 size 600,400
set term gif animate transparent

set output "tokamak.gif"

do for [i=0:360:1] {
    set view 60, i, 2, 2
    splot [-pi:pi][-pi:pi] \
        (6.2 + (2 + 0.6*cos(v))*cos(v))*cos(u), (6.2 + (2 + 0.6*cos(v))*cos(v))*sin(u), 2.8*sin(v) with lines ls 1
}
