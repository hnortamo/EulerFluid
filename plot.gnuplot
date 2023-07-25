# gfortran -O3 -march=native --fast-math  main_bu.f90 && ./a.out && gnuplot plot_gif.gnuplot && firefox foobar.gif
set view map

set term gif size 1500, 1500 animate delay 10 
set output 'foobar.gif'
#set zrange [-10:10]; set cbrange [-10:10];
do for [i=1:500:1]{
	set multiplot layout 3,1 rowsfirst
	splot sprintf("VelX_%04.0f.out",i) matrix with image notitle
	splot sprintf("VelY_%04.0f.out",i) matrix with image notitle
	splot sprintf("Smoke_%04.0f.out",i) matrix with image notitle
}
