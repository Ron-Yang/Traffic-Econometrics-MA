



set encoding iso_8859_1
set style line 1 lt 1 lw 4 pt 7 ps 1.9 lc rgb "#000000"
set style line 2 lt 1 lw 4 lc rgb "#ee0000"
set style line 3 lt 1 lw 2 pt 8 ps 1.9 lc rgb "#0000cc"

set term post eps enhanced color solid "Helvetica" 30
set nokey
set param

set isosample 60,60

set palette defined ( 0  "#bb0000", 5 "red",\
       10 "orange", 20 "yellow", 35 "#66ff66",\
       60 "#6666ff", 100 "#ffffff")
unset surface

set pm3d       # bugfrei/am sichersten:  "set pm3d", DANN "set pm3d map" 
set pm3d; set pm3d map 
set contour surface

set cntrparam bspline 

unset clabel  # dann lauter gleiche Kontourlinien; 
                     # Farbe/Typ mit "w l ls" beim splot-Kommando

set grid   # grid/nogrid: whether 2D grid on xy plane 
set size 1.05,1.5



#################################
set out "freefall_objLandscape.eps"
print "plotting freefall_objLandscape.eps"
#################################

SSEmin=sqrt(320)
SSEmax=sqrt(20000)
dSSE=10  #1000
set xlabel "{/Symbol b}_0 (terminal speed at y=0)[km/h]"
set ylabel "{/Symbol b}_1 (height for 1/e decay)[km]"

set cbrange [SSEmin:SSEmax]
set zrange [:1.1*SSEmax]
#set cbtics dSSE
set cntrparam levels incr SSEmin,dSSE,0.7*SSEmax


set label 1 "sqrt(SSE)" at screen 0.88,1.31 front
splot "simObjfun.dat" u (3.6*$1):(0.001*$2):(sqrt($3))  w l ls 1  


#################################
set out "freefall_vdataSim.eps"
print "plotting freefall_vdataSim.eps"
#################################

unset colorbox
set size 1.5,1.5
unset label 1

set xlabel "Zeit seit Absprung [s]"
set xrange [0:265]
set ylabel "Fallgeschwindigkeit [km/h]" offset 1.5,0
set key

plot[t=0:1]\
  "simResults.dat" u 1:(-3.6*$3) t "Simulation" w l ls 2,\
  "dataBaumgartner.txt" u 1:2 t "Daten" w p ls 1,\
  50,1342 t "max Speed?" w p ls 3,\
  40+t*20,1342 t "" w l ls 3

#################################
set out "freefall_vdata.eps"
print "plotting freefall_vdata.eps"
#################################
plot[t=0:1]\
  "dataBaumgartner.txt" u 1:2 t "Daten" w p ls 1,\
  50,1342 t "max Speed 1342 km/h?" w p ls 3,\
  40+t*20,1342 t "" w l ls 3

#################################
set out "freefall_ydataSim.eps"
print "plotting freefall_ydataSim.eps"
#################################

unset colorbox
set size 1.5,1.5
unset label 1

set xlabel "Zeit seit Absprung [s]"
set xrange [0:265]
set ylabel "Hoehe [km]" offset 1.5,0
set key

plot\
  "simResults.dat" u 1:(0.001*$2) t "Simulation" w l ls 2,\
  0,39.045   t "Daten" w p ls 1,\
  260,2.516  t "" w p ls 1

#################################
set out "freefall_ydata.eps"
print "plotting freefall_ydata.eps"
#################################

plot\
  0,39.045   t "Daten" w p ls 1,\
  260,2.516  t "" w p ls 1
