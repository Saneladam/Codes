set nokey
set parametric
set hidden3d
set view 30
set isosamples 30,20
splot [-pi:pi][-pi:pi] cos(u)*(cos(v)+3), sin(u)*(cos(v)+3), sin(v)
