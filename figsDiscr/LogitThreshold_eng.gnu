


# siehe ~/info/gnuplot, ~/info/gnuTemplate.gnu

##########################################################
#geordnet nach hue (ps)
# set style line: linetype lt, point type pt

set encoding iso_8859_1 # dann -bäöüßÄÖÜ durch woertl. Eingabe korrekt-A

set style line 1 lt 7 lw 6 pt 1 ps 1.2 #schwarz, plus sign
set style line 2 lt 1 lw 8 pt 4 ps 1.2 #rot, open box
set style line 3 lt 8 lw 4 pt 5 ps 1.2 #blassrot, closed square
set style line 4 lt 6 lw 4 pt 6 ps 1.2 #gelb, open circle
set style line 5 lt 2 lw 4 pt 8  ps 1.2 #gruen, open triangle
set style line 6 lt 5 lw 4 pt 9  ps 1.2 #blasstuerkisblau, closed triangle
set style line 7 lt 3 lw 4 pt 10 ps 1.2 #blau, upside-down open triangle
set style line 8 lt 4 lw 4 pt 11 ps 1.2 #lila, upside-down closed triangle

set style line 11 lt 7 lw 1 pt 4 ps 1.2 #rot, open box
set style line 12 lt 1 lw 1 pt 4 ps 1.2 #rot, open box
set style line 13 lt 8 lw 1 pt 5 ps 1.5 #blassrot, closed square
set style line 14 lt 6 lw 1 pt 6 ps 1.5 #gelb, open circle
set style line 15 lt 2 lw 1 pt 8  ps 1.5 #gruen, open triangle
set style line 16 lt 5 lw 1 pt 9  ps 1.8 #blasstuerkisblau, closed triangle
set style line 17 lt 3 lw 1 pt 10 ps 1.5 #blau, upside-down open triangle
set style line 18 lt 4 lw 1 pt 11 ps 1.8 #lila, upside-down closed triangle



max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x


#########################################################
# Logit-Wahrsch. bei zwei Alt, eine exog. Var Reisezeitdiff T, nl Nutzenfkt
#########################################################


V1(dT,b0,b1,b2) =b0+b1*(dT-b2*tanh(dT/b2))  # Delta V, V2=0, b2=plateau width
P1(dT,b0,b1,b2)=1./(exp(-V1(dT,b0,b1,b2))+1)
P2(dT,b0,b1,b2)=1. -P1(dT,b0,b1,b2)

#########################################################
# Daten (symmetrische Auswahlwahrsch. um dT=0)
#########################################################

n=100. # n unabh. Entscheidungen fuer jede Zeitdifferenz
dT1=0;   y11=50; # Zahl der Ja-Antworten beim ersten Choice set (y21=n-y11 usw)
dT2=5;   y21=50; # falls 30, kein Threshold!
dT3=10;  y31=42; # falls 15, kein Threshold!
dT4=15;  y41=10;
dT5=20;  y51=6;
dT6=25;  y61=4;
dT7=30;  y71=3;
dT8=40;  y81=1;
dT9=50;  y91=0;

lnLpart(dT,y1,b0,b1,b2)=y1*log(P1(dT,b0,b1,b2))+(n-y1)*log(P2(dT,b0,b1,b2))\
  + (n-y1)*log(P1(-dT,b0,b1,b2))+y1*log(P2(-dT,b0,b1,b2))

#########################################################
# Log-Likelihood
#########################################################

lnL(b0,b1,b2)=\
   lnLpart(dT1,y11,b0,b1,b2)\
  +lnLpart(dT2,y21,b0,b1,b2)\
  +lnLpart(dT3,y31,b0,b1,b2)\
  +lnLpart(dT4,y41,b0,b1,b2)\
  +lnLpart(dT5,y51,b0,b1,b2)\
  +lnLpart(dT6,y61,b0,b1,b2)\
  +lnLpart(dT7,y71,b0,b1,b2)\
  +lnLpart(dT8,y81,b0,b1,b2)\
  +lnLpart(dT9,y91,b0,b1,b2)




#########################################################
# Plotten
#########################################################

#################################
set out "LogitThreshold_lnL.eps"
print "LogitThreshold_lnL.eps"
#################################


set term post eps enhanced color solid "Helvetica" 30
set nokey
set noparam

set isosample 60,60

set palette defined ( 0 "#ffffff", 20 "#6666ff", \
       40 "#66ff66", 70 "yellow", 90 "orange", 100 "#bb0000")
unset surface

set pm3d       # bugfrei/am sichersten:  "set pm3d", DANN "set pm3d map" 
set pm3d; set pm3d map 
set contour surface

set cntrparam bspline 

unset clabel  # dann lauter gleiche Kontourlinien; 
                     # Farbe/Typ mit "w l ls" beim splot-Kommando

set grid   # grid/nogrid: whether 2D grid on xy plane 
set size 1.5,1.5


b1min=-0.5
b1max=0.
b2min=0.
b2max=20

set label 1 "ln(L)" at screen 1.3,1.31 front
set xlabel "{/Symbol b}_1 [minutes^{-1}]"
set xrange [b1min:b1max]
set ylabel "{/Symbol b}_2 [minutes]"
set yrange [b2min:b2max]

set cbrange [-700:]
unset cntrparam
set cntrparam levels incr -650,10,0
set label 2 at screen 0.25,1.4 front
set label 2 "Utility V1(dT)={/Symbol b}_1*[dT-{/Symbol b}_2*tanh(dT/{/Symbol b}_2)]"
set label 3 at screen 0.25,1.32 front
set label 3 "Threshold in data: f(5min)=50/100, f(-5min)=1-f(5min)"

set label 4 at screen 0.3,0.4 front
set label 4 "n=1800 choices"

splot lnL(0,x,y) w l lt 6 lw 4 

#################################
set out "LogitNoThreshold_lnL.eps"
print "LogitNoThreshold_lnL.eps"
#################################

set label 3 "No Threshold in data: f(5min)=30/100, f(10min)=15/100"
y21=30;
y31=15;
splot lnL(0,x,y) w l lt 6 lw 4 


#################################
set out "LogitNoMaximum_lnL.eps"
print "LogitNoMaximum_lnL.eps"
#################################

set label 3 "Data without unique maximum: f(dT>=30min)=0/100"
y21=50;
y31=50;
y71=0;
y81=0;
splot lnL(0,x,y) w l lt 6 lw 4 


#####################################################################

#################################
set out "LogitThrQuasilin_lnL.eps"
print "LogitThrQuasilin_lnL.eps"
#################################

set ylabel "{/Symbol b}_2 [minutes^{-1}]"

y31=42;
y71=3;
y81=1;

cutofflin(x,a)=(abs(x)>a) ? x/abs(x)*(abs(x)-a) : 0

b1min=-0.3
b1max=0.3
b2min=-0.5
b2max=0.
set xrange [b1min:b1max]
set yrange [b2min:b2max]

V1(dT,b0,b1,b2)=b0+b1*dT+b2*cutofflin(dT,8)

set label 2 "Quasilinear utility {/Symbol b}_1*dT+{/Symbol b}_2*plateauLin(dT=10 min)"
set label 3 "Standard data"
splot lnL(0,x,y) w l lt 6 lw 4 

#################################
set out "LogitNoThrQuasilin_lnL.eps"
print "LogitNoThrQuasilin_lnL.eps"
#################################

set label 3 "No Threshold in data: f(5min)=30/100, f(10min)=15/100"
y21=30;
y31=15;
splot lnL(0,x,y) w l lt 6 lw 4 

#################################
set out "LogitThrQuasilinTestNoMax_lnL.eps"
print "LogitThrQuasilinTestNoMax_lnL.eps"
#################################

y31=50;
y71=0;
y81=0;
set label 3 "Data where tanh function showed no unique maximum"
splot lnL(0,x,y) w l lt 6 lw 4 

#################################
set out "LogitThrQuasilinTestScaleTooSmall_lnL.eps"
print "LogitThrQuasilinTestScaleTooSmall_lnL.eps"
#################################

y31=42;
y71=3;
y81=1;

V1(dT,b0,b1,b2)=b0+b1*dT+b2*cutofflin(dT,0.1)

set label 2 "Test too small plateau scale: 0.1 min instead of 10 min"
set label 3 "Standard data"
splot lnL(0,x,y) w l lt 6 lw 4 


#################################
set out "LogitThrQuasilinTestScaleTooLarge_lnL.eps"
print "LogitThrQuasilinTestScaleTooLarge_lnL.eps"
#################################


V1(dT,b0,b1,b2)=b0+b1*dT+b2*cutofflin(dT,100)

set label 2 "Test too large plateau scale: 100 min instead of 10 min"
splot lnL(0,x,y) w l lt 6 lw 4 
