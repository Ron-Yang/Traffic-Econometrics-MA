# EVA-Funktion nach Lohse und von mir modifiziert

#======== Input  ==============
E=3  # Potenz des asymptotischen Powerlaws B(w) propto 1/w^E
F=5 # Steilheit des Uebergangs von Potenz=0 auf Potenz=E in der Naehe
       # von w0 (je hoeher, desto steiler => desto geringer |B'(w0)|)
w0=23.5 # Widerstand, bei dem die Potenzfunktion Wendepunkt hat (pot=E/2)
#===================================

eva(E,F,w0,w)=(1+w)**(-phiEva(E,F,w0,w))
phiEva(E,F,w0,w)=E/(1+exp(F*(1.-w/w0)))

evaMod(E,F,w0,w)=eva(E,F,1,w/w0)
powerLaw(E,w0,w)=(w<w0) ? 1 : (w0/w)**E

#geordnet nach hue (ps)
# set style line: linetype lt, point type pt

set encoding iso_8859_1 # dann -bäöüßÄÖÜ durch woertl. Eingabe korrekt-A

set style line 1 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000" #schwarz,solid,bullet
set style line 2 lt 1 lw 2 pt 2 ps 1.5  lc rgb "#CC0022" #rot, dash, Kreuz
set style line 3 lt 8 lw 2 pt 4 ps 1.2 #blassrot, offenes Quadrat
set style line 4 lt 6 lw 2 pt 4 ps 1.5  lc rgb "#FFAA00"  #gelb, offenes Quadrat
set style line 5 lt 1 lw 2 pt 5 ps 1.5  lc rgb "#00DD22"  #gruen,solid,ClosedBox
set style line 6 lt 5 lw 2 pt 4 ps 1.5  lc rgb "#00AAAA" #offenes Quadrat
set style line 7 lt 1 lw 2 pt 4 ps 2.0  lc rgb "#1100AA"  #blau,offenes Quadrat
set style line 8 lt 4 lw 2 pt 8 ps 1.5 #lila, aufrechtes geschloss. Dreieck
set style line 9 lt 7 lw 2 pt 9 ps 1.5  lc rgb "#999999"  #grau, aufr. gschl. Dreieck


max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x

set term post eps enhanced color solid "Helvetica" 20
set out "EVA_plot1.eps"
print "plotting EVA_plot1.eps"

set xlabel "w (min)"
set ylabel "B(w) in [0,1]"

set noparam
set key
set title "E=3, F=5, w0=23.5 Minuten"
plot[w=0:80]\
  eva(E,F,w0,w) t "EVA" w l ls 2,\
  evaMod(E,F,w0,w) t "EVAmod" w l ls 3,\
  powerLaw(E,0.7*w0,w) t "Power law (0.7*w0)" w l ls 4,\
  powerLaw(E,w0,w) t "Power law (w0)" w l ls 5

set out "EVA_plot2.eps"
print "plotting EVA_plot2.eps"
set title "E=3, F=2, w0=23.5 Minuten"
F=2
plot[w=0:80]\
  eva(E,F,w0,w) t "EVA" w l ls 2,\
  evaMod(E,F,w0,w) t "EVAmod" w l ls 3,\
  powerLaw(E,0.7*w0,w) t "Power law (0.7*w0)" w l ls 4,\
  powerLaw(E,w0,w) t "Power law (w0)" w l ls 5

set out "EVA_plot3.eps"
print "plotting EVA_plot3.eps"
set title "E=3, F=5, w0=1/2 Stunde"
set xlabel "w (h)"
E=3
F=5
w0=0.5
plot[w=0:3]\
  eva(E,F,w0,w) t "EVA" w l ls 2,\
  evaMod(E,F,w0,w) t "EVAmod" w l ls 3,\
  powerLaw(E,0.7*w0,w) t "Power law (0.7*w0)" w l ls 4,\
  powerLaw(E,w0,w) t "Power law (w0)" w l ls 5
