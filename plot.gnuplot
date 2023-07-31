# gfortran -O3 -march=native --fast-math  main_bu.f90 && ./a.out && gnuplot plot_gif.gnuplot && firefox foobar.gif
set view map

set term gif size 1500, 1500 animate delay 10 
set output 'foobar.gif'
do for [i=1:800:1]{
	set multiplot layout 3,1 rowsfirst
#	set zrange [-3:3]; set cbrange [-3:3];
	splot sprintf("VelX_%04.0f.out",i) matrix with image notitle
	splot sprintf("VelY_%04.0f.out",i) matrix with image notitle
#	set zrange [0:4.5]; set cbrange [0:4.5];
	splot sprintf("Smoke_%04.0f.out",i) matrix with image notitle
}
