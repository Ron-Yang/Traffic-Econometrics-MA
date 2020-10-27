
set encoding iso_8859_1   # dann äöüßÄÖÜ durch woertl. Eingabe korrekt
#\304     "A
#\326     "O
#\334     "U
#\344     "a
#\366     "o
#\374     "u
#\337     &szlig;


set style line 99 lt 1 lw 3 pt 4 ps 1.5 linecolor rgb "#1100EE" #blau, solid, open box

set style line 1 lt 1 lw 2 pt 1 ps 1.5  lc rgb "#000000" #schwarz,solid,plus sign
set style line 2 lt 7 lw 8 pt 2 ps 1.5  lc rgb "#CC0022" #rot, dash, Kreuz
set style line 3 lt 8 lw 5 pt 3 ps 1.5 #blassrot, offener Stern
set style line 4 lt 6 lw 5 pt 4 ps 1.5 #gelb, offenes Quadrat
set style line 5 lt 1 lw 8 pt 5 ps 1.5  lc rgb "#00DD22"  #gruen,solid,ClosedBox
set style line 6 lt 5 lw 5 pt 6 ps 1.5 #blasstuerkisblau, offener Kreis
set style line 7 lt 1 lw 2 pt 7 ps 1.5  lc rgb "#1100AA"  #blau,solid,Bullet
set style line 8 lt 4 lw 2 pt 8 ps 1.5 #lila, aufrechtes geschloss. Dreieck
set style line 9 lt 7 lw 2 pt 9 ps 1.5 #schwarz, aufrechtes geschl. Dreieck
set style line 10 lt 1 lw 2 pt 10 ps 1.5 #rot, upside-down offenes Dreieck
set style line 11 lt 7 lw 10 pt 11 ps 1.5 #schwarz, upside-down geschl. Dreieck
set style line 12 lt 1 lw 2 pt 12 ps 1.5 #rot, offene Raute
set style line 13 lt 8 lw 2 pt 13 ps 1.5 #blassrot, geschl. Raute
set style line 14 lt 6 lw 2 pt 14 ps 1.5 #gelb, "Sonderzeichen" ...
set style line 15 lt 2 lw 2 pt 15 ps 1.5 #gruen, geschl. Fuenfeck
set style line 16 lt 5 lw 2 pt 16 ps 1.5 #blasstuerkisblau, 
set style line 17 lt 3 lw 2 pt 17 ps 1.5 #blau, 
set style line 18 lt 4 lw 2 pt 18 ps 1.5 #lila,
set style line 19 lt 7 lw 2 pt 19 ps 1.5
set style line 20 lt 1 lw 2 pt 20 ps 1.5



############### Beispiele fuer Funktionen ####################

limitedGrowth(t)=ys/(1+ (ys/y0-1)*exp(-t/tau))

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x


set term post eps enhanced color solid "Helvetica" 24

set out "limitedGrowth2.eps"
print "plotting limitedGrowth2.eps ..."
set noparam
set xlabel "exogene Variable (Zeit) x" offset 0,0.5
set ylabel "endogene Variable Y"
set xrange [0:3]
set yrange [0:]
y0=0.5
ys=5.
tau=0.5
plot limitedGrowth(x) w l ls 2
