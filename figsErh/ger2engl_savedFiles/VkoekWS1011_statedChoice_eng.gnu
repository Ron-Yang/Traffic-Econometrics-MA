

# Original gnuplot File in ~/vorlesungen/Verkehrsoekonometrie_Ma/discrChoice_cc
# in anderen Dirs nur kopiert!!


##########################################################
#geordnet nach hue (ps)
# set style line: linetype lt, point type pt

set encoding iso_8859_1 # dann -bäöüßÄÖÜ durch woertl. Eingabe korrekt-A
#\304     "A
#\326     "O
#\334     "U
#\344     "a
#\366     "o
#\374     "u
#\337     &szlig;


set style line 99 lt 1 lw 3 pt 4 ps 1.5 linecolor rgb "#1100EE" #blau, solid, open box

set style line 1 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000" #schwarz,solid,bullet
set style line 2 lt 1 lw 2 pt 5 ps 1.5  lc rgb "#CC0022" #rot, dash, Quadrat
set style line 3 lt 8 lw 2 pt 4 ps 1.2 #blassrot, offenes Quadrat
set style line 4 lt 6 lw 2 pt 4 ps 1.5  lc rgb "#FFAA00"  #gelb, offenes Quadrat
set style line 5 lt 1 lw 2 pt 5 ps 1.5  lc rgb "#00DD22"  #gruen,solid,ClosedBox
set style line 6 lt 5 lw 2 pt 4 ps 1.5  lc rgb "#00AAAA" #offenes Quadrat
set style line 7 lt 1 lw 2 pt 4 ps 2.0  lc rgb "#1100AA"  #blau,offenes Quadrat
set style line 8 lt 4 lw 2 pt 8 ps 1.5 #lila, aufrechtes geschloss. Dreieck
set style line 9 lt 7 lw 2 pt 9 ps 1.5  lc rgb "#999999"  #grau, aufr. gschl. Dreieck



set style line 11 lt 1 lw 6 pt 7 ps 1.9  lc rgb "#000000" #schwarz,solid,bullet
set style line 12 lt 1 lw 6 pt 2 ps 1.5  lc rgb "#CC0022" #rot, dash, Kreuz
set style line 13 lt 8 lw 6 pt 3 ps 1.2 #blassrot, offener Stern
set style line 14 lt 6 lw 6 pt 4 ps 1.5  lc rgb "#FFAA00"  #gelb, offenes Quadrat
set style line 15 lt 1 lw 6 pt 5 ps 1.5  lc rgb "#00DD22"  #gruen,solid,ClosedBox
set style line 16 lt 5 lw 6 pt 7 ps 1.5  lc rgb "#00CCCC" #offener Kreis
set style line 17 lt 1 lw 6 pt 7 ps 1.5  lc rgb "#1100AA"  #blau,solid,Bullet
set style line 18 lt 4 lw 6 pt 8 ps 1.5 #lila, aufrechtes geschloss. Dreieck
set style line 19 lt 7 lw 6 pt 9 ps 1.5  lc rgb "#999999"  #grau, aufr. gschl. Dreieck

#Sinnvolle point typen (pt)
# 1=Plus,2=Kreuz,4=openQuare,5=closedSquare, 6=openCirc,7=closedCirc,
# 9-11=triangles, 12-13=Rauten


############### Beispiele fuer Funktionen ####################


max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x

#########################################################
# Input
#########################################################

#Daten aus  VkoekWS1011_statedChoice.userInput
# Achtung! Punkt be gnuplot noetig!!

nPers=19

  h00=3.;     h10=19-h00;
  h01=6.;       h11=19-h01;
  h02=19.;       h12=19-h02;
  h03=19.;       h13=19-h03;
  h04=17.;       h14=19-h04;
  h05=19.;       h15=19-h05;
  h06=19.;       h16=19-h06;
  h07=14.;       h17=19-h07;
  h08=19.;       h18=19-h08;

  T00=30.;       T10=30.;   
  T01=30.;       T11=40.;  
  T02=30.;       T12=50.;
  T03=30.;       T13=60.; 
  T04=30.;       T14=30.;   
  T05=30.;       T15=30.;   
  T06=30.;       T16=30.;   
  T07=20.;       T17=30.;   
  T08=15.;       T18=30.;   

  K00=0.;     K10=0.; 
  K01=0.;     K11=0.; 
  K02=0.;     K12=0.; 
  K03=0.;     K13=0.; 
  K04=0.;     K14=1.; 
  K05=0.;     K15=2.; 
  K06=0.;     K16=3.; 
  K07=0.;     K17=0.; 
  K08=0.;     K18=0.; 

#beta und Funktionalmatrix J aus Output von discrChoice_kalib

beta0=1.464305; stddev0=1.205988
beta1=-0.408988; stddev1=0.056505
beta2=-0.267504; stddev2=0.042130
beta3=-4.941125; stddev3=0.720795

 J00=12.716348;J01=340.397125;J02=-456.941154;J03=-1.799572;
 J10=340.397125;J11=9436.519835;J12=-12475.434902;J13=-53.987146;
 J20=-456.941154;J21=-12475.434902;J22=17036.847934;J23=53.987146;
 J30=-1.799572;J31=-53.987146;J32=53.987146;J33=1.842485;

#########################################################
# Ende Input
#########################################################


# Modell: Wahrscheinlichkeit fuer Alternative 1 (von 2)

V1(T1,K1)=beta0+beta1*T1
V2(T2,K2)=beta2*T2+beta3*K2
P1(T1,T2,K1,K2)=1./(1+exp(V2(T2,K2)-V1(T1,K1)))
lnP1P2(T1,K1,T2,K2)=V1(T1,K1)-V2(T2,K2)

#########################################################
# Ergebnisse plotten
#########################################################

set term post eps enhanced color solid "Helvetica" 24

set out "VkoekWS1011_statedChoice_1.eps"
print "plotting VkoekWS1011_statedChoice_1.eps"

set xlabel "T_{\326V} (min)"
Tmin=25.
Tmax=55.
set xrange [Tmin:Tmax]

set ylabel "ln(P_{Fu\337,Rad} / P_{\326V} )"
set auto y

set key
plot[t=Tmin:Tmax]\
 t, lnP1P2(T00, 0,t, 0) t "Logit-Modell" w l ls 1,\
 T10, log(h00/h10) t "Umfrage" w p ls 2 ,\
 T11, log(h01/h11) t "" w p ls 2 ,\
 T12, log(h02/h12) t "" w p ls 2 ,\
 T13, log(h03/h13) t "" w p ls 2

###############################
set out "VkoekWS1011_statedChoice_2.eps"
print "plotting VkoekWS1011_statedChoice_2.eps"
###############################


set xlabel "K_{\326V} (Euro)"
Kmin=-0.5
Kmax=3.
set xrange [Kmin:Kmax]

plot[t=Kmin:Kmax]\
 t, lnP1P2(T00, 0,T10, t) t "Logit-Modell" w l ls 1,\
 K10, log(h00/h10) t "Umfrage" w p ls 2 ,\
 K15, log(h04/h14) t "" w p ls 2 ,\
 K16, log(h05/h15) t "" w p ls 2 

###############################
set out "VkoekWS1011_statedChoice_3.eps"
print "plotting VkoekWS1011_statedChoice_3.eps"
###############################

set xlabel  "T_{Fu\337,Rad}(min)"
Tmin=15.
Tmax=35.
set xrange [Tmin:Tmax]

plot[t=Tmin:Tmax]\
 t, lnP1P2(t, 0,T10, 0) t "Logit-Modell" w l ls 1,\
 T00, log(h00/h10) t "Umfrage" w p ls 2 ,\
 T07, log(h07/h17) t "" w p ls 2,\
 T08, log(h08/h18) t "" w p ls 2


###############################
set out "VkoekWS1011_statedChoice_4.eps"
print "plotting VkoekWS1011_statedChoice_4.eps"
###############################

set nolabel
set auto
set xlabel "Situations-Nummer"
set ylabel "Zahl der Fu\337 - und Radnutzer"
hFuss(altern,h)=(altern<0.5) ? h : nPers-h
set key bottom

plot\
  "VkoekWS1011_statedChoice.dat" u ($1+1):(hFuss($2,$5)) t "Daten" w lp ls 11,\
  "VkoekWS1011_statedChoice.dat" u ($1+1):(hFuss($2,$6)) t "Modell" w lp ls 12


###############################
set out "VkoekWS1011_statedChoice_5.eps"
print "plotting VkoekWS1011_statedChoice_5.eps"
###############################
set term post eps enhanced color solid "Helvetica" 20
set key at screen 0.95,0.20
T1(T2,K2,c)=(log(c)+beta2*T2+beta3*K2-beta0)/beta1

set xlabel  "T_{\326V}(min)"
set ylabel  "T_{Fu\337,Rad}(min)"
Tmin=0.
Tmax=60.
set xrange [Tmin:Tmax]
set yrange [0:50]
A1=0.1
A2=0.3
A3=0.5
A4=0.7
A5=0.9

plot[t=Tmin:Tmax]\
 t, T1(t,0,A1/(1-A1)) t "10% Fu\337,Rad" w l ls 2,\
 t, T1(t,0,A2/(1-A2)) t "30% Fu\337,Rad" w l ls 3,\
 t, T1(t,0,A3/(1-A3)) t "50% Fu\337,Rad" w l ls 15,\
 t, T1(t,0,A4/(1-A4)) t "70% Fu\337,Rad" w l ls 6,\
 t, T1(t,0,A5/(1-A5)) t "90% Fu\337,Rad" w l ls 7,\
 t, T1(t,2,A3/(1-A3)) t "50% Fu\337,Rad, 2 Euro \326V-Kosten" w l ls 1






print ""

set term post eps enhanced color solid "Helvetica" 28
set nokey
set noparam

set isosample 60,60

set palette defined ( 0 "#ffffff", 5 "blue", \
      15 "green", 25 "yellow", 40 "orange", 60 "red", 100 "#aa0066")
unset surface

set pm3d       # bugfrei/am sichersten:  "set pm3d", DANN "set pm3d map" 
set pm3d  map 
set contour surface

set cntrparam bspline 
set cntrparam levels 30                          # n equidistant levels from min-max
#set cntrparam levels discrete 0.5,1,1.5,2,4,6,8,10,12,14,16,18,20

unset clabel  # dann lauter gleiche Kontourlinien; 
                     # Farbe/Typ mit "w l ls" beim splot-Kommando

set grid   # grid/nogrid: whether 2D grid on xy plane 
set size 1.5,1.5


##############################
set out "VkoekWS1011_statedChoice_Lbeta1beta2.eps"
print "plotting VkoekWS1011_statedChoice_Lbeta1beta2.eps"

logL(b0,b1,b2,b3)=-(\
 J00*(b0-beta0)**2+J11*(b1-beta1)**2+J22*(b2-beta2)**2+J33*(b3-beta3)**2\
 +2*J01*(b0-beta0)*(b1-beta1)\
 +2*J02*(b0-beta0)*(b2-beta2)\
 +2*J03*(b0-beta0)*(b3-beta3)\
 +2*J12*(b1-beta1)*(b2-beta2)\
 +2*J13*(b1-beta1)*(b3-beta3)\
 +2*J23*(b2-beta2)*(b3-beta3)\
)

db=2; 
set xlabel "{/Symbol b}_1 (Zeit Sensitivit\344t Fu\337/Rad)"
set xrange [beta1-db*stddev1:beta1+db*stddev1]
set ylabel "{/Symbol b}_2 (Zeit-Sensitivit\344t \326V)"
set yrange [beta2-db*stddev2:beta2+db*stddev2]

set noparam
splot exp(logL(beta0,x,y,beta3)) w l ls 1


##############################
set out "VkoekWS1011_statedChoice_Lbeta1beta3.eps"
print "plotting VkoekWS1011_statedChoice_Lbeta1beta3.eps"
set ylabel "{/Symbol b}_3 (Geld-Sensitivit\344t \326V)"
set yrange [beta3-db*stddev3:beta3+db*stddev3]

splot exp(logL(beta0,x,beta2,y)) w l ls 1

##############################
set out "VkoekWS1011_statedChoice_Lbeta2beta3.eps"
print "plotting VkoekWS1011_statedChoice_Lbeta2beta3.eps"
set xlabel "{/Symbol b}_2 (Zeit-Sensitivit\344t \326V)"
set xrange [beta2-db*stddev2:beta2+db*stddev2]
splot exp(logL(beta0,beta1,x,y)) w l ls 1

##############################
set out "VkoekWS1011_statedChoice_Lbeta0beta1.eps"
print "plotting VkoekWS1011_statedChoice_Lbeta0beta1.eps"
set xlabel "{/Symbol b}_0 (Globale Bevorzugung Fu\337/Rad)"
set xrange [beta0-db*stddev0:beta0+db*stddev0]
set ylabel "{/Symbol b}_1 (Zeit Sensitivit\344t Fu\337/Rad)"
set yrange [beta1-db*stddev1:beta1+db*stddev1]
splot exp(logL(x,y,beta2,beta3)) w l ls 1


##############################
set xlabel "{/Symbol b}_2 (Geld-Vorfaktor)"
set xrange [beta2-1:beta2+1]
set ylabel "{/Symbol b}_3 (Rad/Fuss-Bevorzugung)"
set yrange [beta3-2:beta3+2]

print "Unschaerfe Zufallsnutzen in Radminuten = -1/beta1= ",-1/beta1
print "Unschaerfe Zufallsnutzen in OEPNV-Minuten = -1/beta2= ",-1/beta2
print "Unschaerfe Zufallsnutzen in OEPNV-Euro = -1/beta3= ",-1/beta3
print "Bonus Fuss/Rad in OEPNV-Minuten = -beta0/beta2= ", -beta0/beta2
print "Zeit-Geld-relation OEPNV in Euro/h = 60*beta2/beta3= ",60*beta2/beta3










