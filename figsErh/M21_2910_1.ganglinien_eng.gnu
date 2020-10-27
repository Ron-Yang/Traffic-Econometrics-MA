# (pointstyle (s.u.):p 5 = closed square, p 9 = closed triangle)
#######################################
set style line 1 lt 7 lw 6 pt 9  ps 1.7 #7=schwarz 
set style line 2 lt 3 lw 6 pt 10 ps 1.7 #3=blau
set style line 3 lt 5 lw 6 pt 11 ps 1.7 #5=tuerkisgruen
set style line 4 lt 2 lw 6 pt 12 ps 1.7 #2=gruen
set style line 5 lt 6 lw 6 pt 5  ps 1.7 #6=gelb 
set style line 6 lt 8 lw 6 pt 6  ps 1.7 #8 art blassrot
set style line 7 lt 1 lw 6 pt 7  ps 1.7 #1=rot
set style line 8 lt 4 lw 6 pt 8  ps 1.7 #4=lila

set style line 10 lt 7 lw 2 pt 7  ps 1.7 #7=schwarz 
set style line 20 lt 3 lw 2 pt 10 ps 1.7 #3=blau
set style line 30 lt 5 lw 2 pt 11 ps 1.7 #5=tuerkisgruen
set style line 40 lt 2 lw 2 pt 12 ps 1.7 #2=gruen
set style line 50 lt 6 lw 2 pt 5  ps 1.7 #6=gelb 
set style line 60 lt 8 lw 2 pt 6  ps 1.7 #8 art blassrot
set style line 70 lt 1 lw 2 pt 7  ps 1.7 #1=rot
set style line 80 lt 4 lw 2 pt 8  ps 1.7 #4=lila
#######################################

#set nogrid
set grid

set term post eps enhanced color solid "Helvetica" 24
set out "M21_2910_1.fund.eps"
print "plotting M21_2910_1.fund.eps"

set xtics 20
set ytics 500
set size 0.8,1
set xrange [0:80]
set xlabel "Density (vehicles/km/lane)"
set ylabel "Flow (vehicles/h/lane)" offset 0.5,0

plot "M21_2910_1.dat" u ($2):($2*$3) w linesp ls 10


set out "M21_2910_1.ganglinie.v.eps"
print "plotting M21_2910_1.ganglinie.v.eps"

set size 1,0.6

set xlabel "t (h)"
set xrange [0:24]
set xtics 4

set ylabel "V (km/h)"
set auto y
set ytics 20
plot "M21_2910_1.dat" u ($1/60):($3) w l ls 70

#####################################

set out "M21_2910_1.ganglinieDetail.v.eps"
print "plotting M21_2910_1.ganglinieDetail.v.eps"

set term post eps enhanced color solid "Helvetica" 36
set size 1,0.9

set xrange [6.8:10.5]
set xtics 1
plot "M21_2910_1.dat" u ($1/60):($3) w l ls 7

#####################################

set out "M21_2910_1.ganglinie.Q.eps"
print "plotting M21_2910_1.ganglinie.Q.eps"

set term post eps enhanced color solid "Helvetica" 24
set size 1,0.6

set xrange [0:24]
set xtics 4

set ylabel "Fluss (Fz/h/Spur)"
set auto y
set ytics 500

plot "M21_2910_1.dat" u ($1/60):($2*$3) w l ls 20
