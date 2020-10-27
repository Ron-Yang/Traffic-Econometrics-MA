
#\304     "A
#\326     "O
#\334     "U
#\344     "a
#\366     "o
#\374     "u
#\337     &szlig;

# Original gnuplot File in ~/vorlesungen/Verkehrsoekonometrie_Ma/discrChoice_cc
# !!!!!!!!!

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x
delta (i,j)=(i==j) ? 1 : 0

# 0=Fuss, 1=Rad, 2=OEV, 3=MIV, 1=Referenzalternative
beta0=       1.17794; sig0=0.34344
beta1=       2.28040; sig1=0.36586
beta2=       1.56634; sig2=0.33746
beta3=       -0.20034; sig3=0.02742
beta4=       -2.49928; sig4=0.61021
beta5=       -2.22591; sig5=0.54908
lnLmax=-428.555

V(i,dT1, dT2, dT3, dK2, dK3, WetterSchoen)=\
   beta0*delta(i,1)+ beta1*delta(i,2)+ beta2*delta(i,3)\
 + beta3*dT1*delta(i,1)+ beta3*dT2*delta(i,2)+ beta3*dT3*delta(i,3)\
 + beta4*dK2*delta(i,2)+ beta4*dK3*delta(i,3)\
 + beta5*WetterSchoen*(delta(i,2)+delta(i,3))

denom(dT1, dT2, dT3, dK2, dK3, WetterSchoen)=\
   exp(V(0,dT1, dT2, dT3, dK2, dK3, WetterSchoen))\
 + exp(V(1,dT1, dT2, dT3, dK2, dK3, WetterSchoen))\
 + exp(V(2,dT1, dT2, dT3, dK2, dK3, WetterSchoen))\
 + exp(V(3,dT1, dT2, dT3, dK2, dK3, WetterSchoen))

Prob(i,dT1, dT2, dT3, dK2, dK3, WetterSchoen)=\
   exp(V(i,dT1, dT2, dT3, dK2, dK3, WetterSchoen))\
 / denom(dT1, dT2, dT3, dK2, dK3, WetterSchoen)

print "AC Rad-Fuss in min = beta0/beta3= ",beta0/beta3
print "AC OEV-Fuss in min = beta1/beta3= ",beta1/beta3
print "AC MIV-Fuss in min = beta2/beta3= ",beta2/beta3
print "Zeit-Geld-Relation in Euro/h = 60*beta3/beta4= ",60*beta3/beta4
print "Schoen vs Schlechtwetter in min = beta5/beta3= ",beta5/beta3
print "Schoen vs Schlechtwetter in Euro = beta5/beta4= ",beta5/beta4


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
set style line 5 lt 1 lw 2 pt 8 ps 1.8  lc rgb "#00DD22"  #gruen,solid,ClosedBox
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






set term post eps enhanced color solid "Helvetica" 24

###############################
set out "statedChoiceWS1213_fProb.eps"
print "plotting statedChoiceWS1213_fProb.eps"
###############################

set nolabel
set auto
set xlabel "Choice Set"
set ylabel "Relative H\344ufigkeiten"

set key top

set yrange [0:1.3]
set xrange [0.9:]
plot\
  "statedChoiceWS1213.outData" u ($1):($2/$10) t "Daten Fu\337" w p ls 7,\
  "statedChoiceWS1213.outData" u ($1):($3/$10) t "Daten Rad" w p ls 5,\
  "statedChoiceWS1213.outData" u ($1):($4/$10) t "Daten \326V" w p ls 2,\
  "statedChoiceWS1213.outData" u ($1):($6/$10) t "Modell Fu\337" w l ls 7,\
  "statedChoiceWS1213.outData" u ($1):($7/$10) t "Modell Rad" w l ls 5,\
  "statedChoiceWS1213.outData" u ($1):($8/$10) t "Modell \326V" w l ls 2



set term post eps enhanced color solid "Helvetica" 20

###############################
set out        "statedChoiceWS1213_ProbOEV_T.eps"
print "plotting statedChoiceWS1213_ProbOEV_T.eps"
###############################

dT1=0; dT3=0; dK3=0;
dT2min=-15
dT2max=15
set xlabel "T_{\326V}-T_{Rad,Fu\337,MIV} [min]"
set xrange [dT2min:dT2max]
set ylabel "P_{\326V}
set yrange [0:1.2]

plot[dT2=dT2min:dT2max]\
  dT2, Prob(2,dT1, dT2, dT3, 0, dK3, 0)\
    t "Wetter novembergrau" w l ls 16,\
  dT2, Prob(2,dT1, dT2, dT3, 1, dK3, 0)\
    t "Wetter novembergrau, 1 Euro {\326V}-Kosten" w l ls 17,\
  dT2, Prob(2,dT1, dT2, dT3, 0, dK3, 1)\
    t "Wetter sommerlich schoen" w l ls 2,\
  dT2, Prob(2,dT1, dT2, dT3, 1, dK3, 1)\
    t "Wetter sommerlich schoen, 1 Euro {\326V}-Kosten" w l ls 3

###############################
set out        "statedChoiceWS1213_Prob4Modi_TOEV.eps"
print "plotting statedChoiceWS1213_Prob4Modi_TOEV.eps"
###############################

plot[dT2=dT2min:dT2max]\
  dT2, Prob(0,dT1, dT2, dT3, 0, dK3, 0)\
    t "Fuss, novembergrau, keine OEV-Kosten" w l ls 17,\
  dT2, Prob(0,dT1, dT2, dT3, 0, dK3, 0)+Prob(1,dT1, dT2, dT3, 0, dK3, 0)\
    t "Fuss+Rad" w l ls 15,\
  dT2, 1-Prob(3,dT1, dT2, dT3, 0, dK3, 0)\
    t "Fuss+Rad+OEV" w l ls 12

###############################
set out        "statedChoiceWS1213_Prob4Modi_Entf_wetterSchoen.eps"
print "plotting statedChoiceWS1213_Prob4Modi_Entf_wetterSchoen.eps"
###############################
wetter=1
TruestOEV=0.
TruestMIV=0.
v0_kmh=5.
v1_kmh=15.; dTentf1(r)=60*r*(1./v1_kmh-1./v0_kmh)
v2_kmh=20.; dTentf2(r)=60*r*(1./v2_kmh-1./v0_kmh)+TruestOEV
v3_kmh=25.; dTentf3(r)=60*r*(1./v3_kmh-1./v0_kmh)+TruestMIV
rmin=0.
rmax=10.

set label 1 "Angenommene Geschwindikeiten:" at screen 0.5,0.75
set label 2 "Fuss 5 km/h, Rad 15 km/h" at screen 0.6,0.70
set label 3 "OEV 20 km/h, MIV 25 km/h" at screen 0.6,0.65

set xlabel "Entfernung (km):
set xrange [rmin:rmax]


plot[r=rmin:rmax]\
  r, Prob(0,dTentf1(r), dTentf2(r), dTentf3(r), 0, dK3, wetter)\
    t "Fuss, Schoenwetter, keine OEV-Kosten" w l ls 17,\
  r, Prob(0,dTentf1(r), dTentf2(r), dTentf3(r), 0, dK3, wetter)\
   + Prob(1,dTentf1(r), dTentf2(r), dTentf3(r), 0, dK3, wetter)\
    t "Fuss+Rad" w l ls 15,\
  r, 1-Prob(3,dTentf1(r), dTentf2(r), dTentf3(r), 0, dK3, wetter)\
    t "Fuss+Rad+OEV" w l ls 12

###############################
set out        "statedChoiceWS1213_Prob4Modi_Entf_wetterSchlecht.eps"
print "plotting statedChoiceWS1213_Prob4Modi_Entf_wetterSchlecht.eps"
###############################
wetter=0
plot[r=rmin:rmax]\
  r, Prob(0,dTentf1(r), dTentf2(r), dTentf3(r), 0, dK3, wetter)\
    t "Fuss, Novembergrau, keine OEV-Kosten" w l ls 17,\
  r, Prob(0,dTentf1(r), dTentf2(r), dTentf3(r), 0, dK3, wetter)\
   + Prob(1,dTentf1(r), dTentf2(r), dTentf3(r), 0, dK3, wetter)\
    t "Fuss+Rad" w l ls 15,\
  r, 1-Prob(3,dTentf1(r), dTentf2(r), dTentf3(r), 0, dK3, wetter)\
    t "Fuss+Rad+OEV" w l ls 12

unset label

###############################
set out        "statedChoiceWS1213_AequiRatio_Rad_OEV.eps"
print "plotting statedChoiceWS1213_AequiRatio_Rad_OEV.eps"
###############################

dT2min=-25.
dT2max=15.
Kmin=0.
Kmax=3.

OEVcost(PRadPOEV, dT2,W)=1./beta4*(beta0-beta1-beta3*dT2-beta5*W-log(PRadPOEV))
set xlabel "T_{\326V}-T_{Rad,Fu\337,MIV} [min]"
set xrange [dT2min:dT2max]

set ylabel "AD-Hoc-Kosten OEV
set yrange [Kmin:Kmax]

plot [dT2=dT2min:dT2max]\
dT2, OEVcost(0.1, dT2, 0) t "Rad:OEV=1:10" w l ls 12,\
dT2, OEVcost(0.1, dT2, 1) t "Rad:OEV=1:10, Sommerwetter" w l ls 2,\
dT2, OEVcost(0.25, dT2, 0) t "Rad:OEV=1:4" w l ls 13,\
dT2, OEVcost(1, dT2, 0) t "Rad:OEV=1:1" w l ls 15,\
dT2, OEVcost(4, dT2, 0) t "Rad:OEV=4:1" w l ls 17,\
dT2, OEVcost(4, dT2, 1) t "Rad:OEV=4:1, Sommerwetter" w l ls 7






########################################################
# Contour plotting
########################################################


set term post eps enhanced color solid "Helvetica" 28
set nokey
set noparam
set grid   # grid/nogrid: whether 2D grid on xy plane 
set size 1.5,1.5
set autoscale

################## Labels #############################
beta0label="{/Symbol b}_1 (AC Rad-Fu\337)"
beta1label="{/Symbol b}_2 (AC \326-Fu\337)"
beta2label="{/Symbol b}_3 (AC MIV-Fu\337)"
beta3label="{/Symbol b}_4 (Zeitsensitivit\344t)"
beta4label="{/Symbol b}_5 (Kostensensitivit\344t)"
beta5label="{/Symbol b}_6 (Dummy Wetter)"
################################################################


set palette defined ( 0 "#ffffff", 30 "#bbbbff", \
      50 "#44ff44", 68 "yellow", 85 "orange", 97 "red", 100 "#aa0066")

unset surface

set pm3d       # bugfrei/am sichersten:  "set pm3d", DANN "set pm3d map" 
set pm3d  map 
set contour surface
set cntrparam bspline 
#set cntrparam levels 30                          # n equidistant levels from min-max
#set cntrparam levels discrete 0.5,1,1.5,2,4,6,8,10,12,14,16,18,20
unset clabel  # dann lauter gleiche Kontourlinien; 
                     # Farbe/Typ mit "w l ls" beim splot-Kommando




##############################
set out "statedChoiceWS1213_L_beta0_beta1.eps"
print "plotting statedChoiceWS1213_L_beta0_beta1.eps"
##############################

set label 1 "L" at screen 1.3,1.33
set xlabel beta0label
set ylabel beta1label
set xrange [beta0-2*sig0:beta0+2*sig0]
set yrange [beta1-2*sig1:beta1+2*sig1]
set cbrange [0:1.01]
set zrange [0:1.01]
set cntrparam levels incr 0.01,0.05,1

splot "statedChoiceWS1213_beta0_beta1" u 1:2:(exp($3-lnLmax))  w l ls 1


set auto x
set auto y

set label 1 "ln L" at screen 1.3,1.33
set cbrange [lnLmax-20:lnLmax+0.1]
set zrange [lnLmax-50:lnLmax+0.1]
set cntrparam levels incr lnLmax-12,1,lnLmax


##############################
set out "statedChoiceWS1213_lnL_beta0_beta1.eps"
print "plotting statedChoiceWS1213_lnL_beta0_beta1.eps"
##############################
set xlabel beta0label
set ylabel beta1label
splot "statedChoiceWS1213_beta0_beta1" u 1:2:($3)  w l ls 1

#quit falls Mparam=2


##############################
set out "statedChoiceWS1213_lnL_beta0_beta2.eps"
print "plotting statedChoiceWS1213_lnL_beta0_beta2.eps"
##############################
set xlabel beta0label
set ylabel beta2label
splot "statedChoiceWS1213_beta0_beta2" u 1:2:($3)  w l ls 1

##############################
set out "statedChoiceWS1213_lnL_beta1_beta2.eps"
print "plotting statedChoiceWS1213_lnL_beta1_beta2.eps"
##############################
set xlabel beta1label
set ylabel beta2label
splot "statedChoiceWS1213_beta1_beta2" u 1:2:($3)  w l ls 1

#quit falls Mparam=3


##############################
set out "statedChoiceWS1213_lnL_beta0_beta3.eps"
print "plotting statedChoiceWS1213_lnL_beta0_beta3.eps"
##############################
set xlabel beta0label
set ylabel beta3label
splot "statedChoiceWS1213_beta0_beta3" u 1:2:($3)  w l ls 1

##############################
set out "statedChoiceWS1213_lnL_beta1_beta3.eps"
print "plotting statedChoiceWS1213_lnL_beta1_beta3.eps"
##############################
set xlabel beta1label
set ylabel beta3label
splot "statedChoiceWS1213_beta1_beta3" u 1:2:($3)  w l ls 1

##############################
set out "statedChoiceWS1213_lnL_beta2_beta3.eps"
print "plotting statedChoiceWS1213_lnL_beta2_beta3.eps"
##############################
set xlabel beta2label
set ylabel beta3label
splot "statedChoiceWS1213_beta2_beta3" u 1:2:($3)  w l ls 1

#quit falls Mparam=4



##############################
set out "statedChoiceWS1213_lnL_beta0_beta4.eps"
print "plotting statedChoiceWS1213_lnL_beta0_beta4.eps"
##############################
set xlabel beta0label
set ylabel beta4label
splot "statedChoiceWS1213_beta0_beta4" u 1:2:($3)  w l ls 1

##############################
set out "statedChoiceWS1213_lnL_beta1_beta4.eps"
print "plotting statedChoiceWS1213_lnL_beta1_beta4.eps"
##############################
set xlabel beta1label
set ylabel beta4label
splot "statedChoiceWS1213_beta1_beta4" u 1:2:($3)  w l ls 1

##############################
set out "statedChoiceWS1213_lnL_beta2_beta4.eps"
print "plotting statedChoiceWS1213_lnL_beta2_beta4.eps"
##############################
set xlabel beta2label
set ylabel beta4label
splot "statedChoiceWS1213_beta2_beta4" u 1:2:($3)  w l ls 1

##############################
set out "statedChoiceWS1213_lnL_beta3_beta4.eps"
print "plotting statedChoiceWS1213_lnL_beta3_beta4.eps"
##############################
set xlabel beta3label
set ylabel beta4label
splot "statedChoiceWS1213_beta3_beta4" u 1:2:($3)  w l ls 1

#quit falls Mparam=5


##############################
set out "statedChoiceWS1213_lnL_beta0_beta5.eps"
print "plotting statedChoiceWS1213_lnL_beta0_beta5.eps"
##############################
set xlabel beta0label
set ylabel beta5label
splot "statedChoiceWS1213_beta0_beta5" u 1:2:($3)  w l ls 1

##############################
set out "statedChoiceWS1213_lnL_beta1_beta5.eps"
print "plotting statedChoiceWS1213_lnL_beta1_beta5.eps"
##############################
set xlabel beta1label
set ylabel beta5label
splot "statedChoiceWS1213_beta1_beta5" u 1:2:($3)  w l ls 1

##############################
set out "statedChoiceWS1213_lnL_beta2_beta5.eps"
print "plotting statedChoiceWS1213_lnL_beta2_beta5.eps"
##############################
set xlabel beta2label
set ylabel beta5label
splot "statedChoiceWS1213_beta2_beta5" u 1:2:($3)  w l ls 1

##############################
set out "statedChoiceWS1213_lnL_beta3_beta5.eps"
print "plotting statedChoiceWS1213_lnL_beta3_beta5.eps"
##############################
set xlabel beta3label
set ylabel beta5label
splot "statedChoiceWS1213_beta3_beta5" u 1:2:($3)  w l ls 1

##############################
set out "statedChoiceWS1213_lnL_beta4_beta5.eps"
print "plotting statedChoiceWS1213_lnL_beta4_beta5.eps"
##############################
set xlabel beta4label
set ylabel beta5label
splot "statedChoiceWS1213_beta4_beta5" u 1:2:($3)  w l ls 1

#quit falls Mparam=6











