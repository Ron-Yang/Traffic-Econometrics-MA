print "hint: Update lnLmax and beta from SC_WS1819_3alt_altspecT_weather.outData"

######################################################
# Input aus SC_WS1819_3alt_altspecT_weather.outData
######################################################


lnLinit=-176.877
lnLmax=-120.516
N=161        #(number of single choices)
Mparam=7
beta0=1.03293; sig0=0.735792
beta1=0.655507; sig1=0.399473
beta2=-0.530675; sig2=0.250633
beta3=-0.138102; sig3=0.0335814
beta4=-0.108261; sig4=0.0294114
beta5=-0.0592279; sig5=0.0261952
beta6=3.57852; sig6=1.14594

w_stddev=3  # number of stddev plotted in contour plots


######################################################
# Formuliere Model
######################################################

# 0=Ped/Bike, 1=PT/Car

beta0label="{/Symbol b}_0 (AC Ped - PT/Car)"
beta1label="{/Symbol b}_1 (AC Bike - PT/Car)"
beta2label="{/Symbol b}_2 (Cost sensitivity)"
beta3label="{/Symbol b}_3 (Time sensitivity Ped)"
beta4label="{/Symbol b}_4 (Time sensitivity Bike)"
beta5label="{/Symbol b}_5 (Time sensitivity PT/Car)"
beta6label="{/Symbol b}_6 (Sensitivity weather)"

#W=0: schoen, W=1: schlecht

V(i,T0,T1,T2,C,W)=\
   beta0*delta(i,0)+ beta1*delta(i,1)\
 + beta2*C*delta(i,2)\
 + beta3*T0*delta(i,0)+beta4*T1*delta(i,1)+beta5*T2*delta(i,2)\
 + beta6*W*delta(i,2)

denom(T0,T1,T2,C,W)=exp(V(0,T0,T1,T2,C,W))+exp(V(1,T0,T1,T2,C,W))+exp(V(2,T0,T1,T2,C,W))

Prob(i,T0,T1,T2,C,W)=exp(V(i,T0,T1,T2,C,W))/denom(T0,T1,T2,C,W)

cost(P2duchrP1,T1,T2,W)=1./beta2*(log(P2duchrP1)-(beta3*T2-beta3*T1)+beta1-beta4*W)

print "AC Ped-PT/Car in ped-min  = -beta0/beta3= ",-beta0/beta3
print "AC Bike-PT/Car  in ped-min  = -beta1/beta3= ",-beta1/beta3
print "AC Ped-PT/Car in Euro = -beta0/beta2= ",-beta0/beta2
print "AC Bike-PT/Car  in Euro = -beta1/beta2= ",-beta1/beta2
print "PT-Zeitwert in Euro/h = 60*beta5/beta2= ",60*beta5/beta2
print "Wert weatherdummy PT/Car vs Ped/Bike in Euro = -beta6/beta2= ",-beta6/beta2

v0_kmh=5.;  T0entf(r)=60*r/v0_kmh
v1_kmh=15.; T1entf(r)=60*r/v1_kmh
v2_kmh=25.; T2entf(r)=60*r/v2_kmh

# Ende User-Iput
######################################################
 
####################################################################
set encoding iso_8859_1 
#\304     "A
#\326     "O
#\334     "U
#\344     "a
#\366     "o
#\374     "u
#\337     &szlig;

# Original gnuplot File in ~/vorlesungen/Verkehrsoekonometrie_Ma/discrChoice_cc_Levmar
# !!!!!!!!!

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x
delta (i,j)=(i==j) ? 1 : 0




set style line 99 lt 1 lw 1 pt 4 ps 1.5 linecolor rgb "#1100EE" #blau, solid, open box

set style line 1 lt 1 lw 1 pt 7 ps 1.9  lc rgb "#000000" #schwarz,solid,bullet
set style line 2 lt 1 lw 1 pt 5 ps 1.5  lc rgb "#CC0022" #rot, dash, Quadrat
set style line 3 lt 8 lw 1 pt 4 ps 1.2 #blassrot, offenes Quadrat
set style line 4 lt 6 lw 1 pt 4 ps 1.5  lc rgb "#FFAA00"  #gelb, offenes Quadrat
set style line 5 lt 1 lw 1 pt 8 ps 1.8  lc rgb "#00DD22"  #gruen,solid,ClosedBox
set style line 6 lt 5 lw 1 pt 4 ps 1.5  lc rgb "#00AAAA" #offenes Quadrat
set style line 7 lt 1 lw 1 pt 4 ps 2.0  lc rgb "#1100AA"  #blau,offenes Quadrat
set style line 8 lt 4 lw 1 pt 8 ps 1.5 #lila, aufrechtes gebadoss. Dreieck
set style line 9 lt 7 lw 1 pt 9 ps 1.5  lc rgb "#999999"  #grau, aufr. gbad. Dreieck



set style line 11 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000" #schwarz,solid,bullet
set style line 12 lt 1 lw 2 pt 2 ps 1.5  lc rgb "#CC0022" #rot, dash, Kreuz
set style line 13 lt 8 lw 2 pt 3 ps 1.2 #blassrot, offener Stern
set style line 14 lt 6 lw 2 pt 4 ps 1.5  lc rgb "#FFAA00"  #gelb, offenes Quadrat
set style line 15 lt 1 lw 2 pt 5 ps 1.5  lc rgb "#00DD22"  #gruen,solid,ClosedBox
set style line 16 lt 5 lw 2 pt 7 ps 1.5  lc rgb "#00CCCC" #offener Kreis
set style line 17 lt 1 lw 2 pt 7 ps 1.5  lc rgb "#1100AA"  #blau,solid,Bullet
set style line 18 lt 4 lw 2 pt 8 ps 1.5  lc rgb "#FF00FF"
set style line 19 lt 7 lw 2 pt 9 ps 1.5  lc rgb "#999999"  #grau, aufr. gbad. Dreieck

#Sinnvolle point typen (pt)
# 1=Plus,2=Kreuz,4=openQuare,5=closedSquare, 6=openCirc,7=closedCirc,
# 9-11=triangles, 12-13=Rauten






set term pngcairo enhanced color notransparent crop font "Helvetica, 14"

###############################
set out "SC_WS1819_3alt_altspecT_weather_fProb.png"
print "plotting SC_WS1819_3alt_altspecT_weather_fProb.png"
###############################

set nolabel
set auto
set xlabel "Choice Set"
set ylabel "Modal split"

set key top

set yrange [0:1.25]
set xrange [0.9:]
plot[t=0:1]\
  "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather.outData" u ($1):($2/$8) t "Data Ped" w p ls 17,\
  "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather.outData" u ($1):($5/$8) t "Model Ped" w l ls 17,\
  "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather.outData" u ($1):(($2+$3)/$8)\
     t "Data Ped+Bike" w p ls 2,\
  "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather.outData" u ($1):(($5+$6)/$8)\
     t "Model Ped+Bike" w l ls 12,\
  12*t, 1 t "" w l ls 11


set term pngcairo enhanced color notransparent crop font "Helvetica, 14"

###############################
set out        "SC_WS1819_3alt_altspecT_weather_ProbPT_T.png"
print "plotting SC_WS1819_3alt_altspecT_weather_ProbPT_T.png"
###############################


dTmin=-30
dTmax=30
T0=30.
T1=30.
set title sprintf("T_{Ped}=%2.0f, T_{Bike}=%2.0f", T0, T1)


set xlabel "T_{PT/Car} [min]"
set xrange [T0+dTmin:T0+dTmax]
set ylabel "Prob(PT/Car)"
set yrange [0:1.2]

plot[T2=T0+dTmin:T0+dTmax]\
  T2, Prob(2,T0,T1,T2,0,0) t  "0 Euro {PT}-costs" w l ls 17,\
  T2, Prob(2,T0,T1,T2,1,0) t  "1 Euro {PT}-costs" w l ls 15,\
  T2, Prob(2,T0,T1,T2,2,0) t  "2 Euro {PT}-costs" w l ls 12,\
  T2, Prob(2,T0,T1,T2,2,1) t  "2 Euro, bad weather" w l ls 18


set out        "SC_WS1819_3alt_altspecT_weather_ProbPT_T1.png"
print "plotting SC_WS1819_3alt_altspecT_weather_ProbPT_T1.png"
plot[T2=T0+dTmin:T0+dTmax]\
  T2, Prob(2,T0,T1,T2,0,0) t  "0 Euro {PT}-costs" w l ls 17

set out        "SC_WS1819_3alt_altspecT_weather_ProbPT_T2.png"
print "plotting SC_WS1819_3alt_altspecT_weather_ProbPT_T2.png"
plot[T2=T0+dTmin:T0+dTmax]\
  T2, Prob(2,T0,T1,T2,0,0) t  "0 Euro {PT}-costs" w l ls 17,\
  T2, Prob(2,T0,T1,T2,1,0) t  "1 Euro {PT}-costs" w l ls 15

set out        "SC_WS1819_3alt_altspecT_weather_ProbPT_T3.png"
print "plotting SC_WS1819_3alt_altspecT_weather_ProbPT_T3.png"
plot[T2=T0+dTmin:T0+dTmax]\
  T2, Prob(2,T0,T1,T2,0,0) t  "0 Euro {PT}-costs" w l ls 17,\
  T2, Prob(2,T0,T1,T2,1,0) t  "1 Euro {PT}-costs" w l ls 15,\
  T2, Prob(2,T0,T1,T2,2,0) t  "2 Euro {PT}-costs" w l ls 12


###############################
set out        "SC_WS1819_3alt_altspecT_weather_ProbPT_Entf.png"
print "plotting SC_WS1819_3alt_altspecT_weather_ProbPT_Entf.png"
###############################


C_PT=1.0
str_vel=sprintf("v_{Ped}=%1.1f km/h, v_{Bike}=%1.1f km/h, v_{PT}=%1.1f km/h",\
   v0_kmh,v1_kmh,v2_kmh)

str_title=sprintf("PT-costs %1.1f Euro", C_PT)


rmin=0.
rmax=6.

set title str_title
# set label 1 str_vel at screen 0.2,0.66

set xlabel "Distance [km]"
set xrange [rmin:rmax]
set yrange [0:1.5]
set ylabel "Probability"


plot[r=rmin:rmax]\
  r, Prob(0,T0entf(r),T1entf(r),T2entf(r),C_PT,0) t "Ped" w l ls 17,\
  r, Prob(0,T0entf(r),T1entf(r),T2entf(r),C_PT,1) t "Ped, bad weather" w l ls 16,\
  r, 1-Prob(2,T0entf(r),T1entf(r),T2entf(r),C_PT,0) t "Ped+Bike" w l ls 12,\
  r, 1-Prob(2,T0entf(r),T1entf(r),T2entf(r),C_PT,1) t "Ped+Bike, bad weather" w l ls 13,\
  r, 1 t "" w l ls 11


unset label
unset title

###############################
set out        "SC_WS1819_3alt_altspecT_weather_AequiRatio_PT_Bike.png"
print "plotting SC_WS1819_3alt_altspecT_weather_AequiRatio_PT_Bike.png"
###############################

T0=30. # nicht verwendet
T1=30.
T2=30. # nicht verwendet
dTmin=-20.
dTmax=40.
Cmin=0.
Cmax=3.

set title sprintf("Linien  P_{PT}/P_{Bike}=const bei T_{Bike}=%2.0f Minuten",T1);
set xlabel "Differenz der komplexe Reisezeit T_{PT}-T_{Bike} [min]"
set xrange [dTmin:dTmax]

set ylabel "AD-Hoc-costs PT
set yrange [Cmin:Cmax]

plot [dT=dTmin:dTmax]\
dT, cost(0.1,T1, T1+dT,0) t "P_{PT}/P_{Bike}=1:10" w l ls 12,\
dT, cost(1,  T1, T1+dT,0) t "P_{PT}/P_{Bike}=1:1"  w l ls 15,\
dT, cost(10, T1, T1+dT,0) t "P_{PT}/P_{Bike}=10:1" w l ls 17,\
dT, cost(1,  T1, T1+dT,1) t "P_{PT}/P_{Bike}=1:1, weather badecht" w l ls 5

###############################
set out        "SC_WS1819_3alt_altspecT_weather_AequiRatio_PT_Bike_highT.png"
print "plotting SC_WS1819_3alt_altspecT_weather_AequiRatio_PT_Bike_highT.png"
###############################
T1=60.
set title sprintf("Linien  P_{PT}/P_{Bike}=const bei T_{Bike}=%2.0f Minuten",T1);
plot [dT=dTmin:dTmax]\
dT, cost(0.1,T1, T1+dT,0) t "P_{PT}/P_{Bike}=1:10" w l ls 12,\
dT, cost(1,  T1, T1+dT,0) t "P_{PT}/P_{Bike}=1:1"  w l ls 15,\
dT, cost(10, T1, T1+dT,0) t "P_{PT}/P_{Bike}=10:1" w l ls 17,\
dT, cost(1,  T1, T1+dT,1) t "P_{PT}/P_{Bike}=1:1, weather badecht" w l ls 5

unset title

########################################################
print "\nGoodnessof-Fit measures"

AIC=-2*(lnLmax+Mparam*N/(N-Mparam-1))
BIC=-2*lnLmax+Mparam*log(N)
rho2=1-lnLmax/lnLinit
adjrho2=1-(lnLmax-Mparam)/lnLinit
print "AIC=",AIC
print "BIC=",BIC
print "rho2=",rho2
print "adjrho2=",adjrho2
print ""

########################################################
# Contour plotting
########################################################


set term pngcairo enhanced color notransparent crop font "Helvetica, 14"
set nokey
set noparam
set grid   # grid/nogrid: whether 2D grid on xy plane 
set size 1.0,1.0
set autoscale


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
set out "SC_WS1819_3alt_altspecT_weather_L_beta2_beta3.png"
print "plotting SC_WS1819_3alt_altspecT_weather_L_beta2_beta3.png"
##############################

set label 1 "L/L_{max}" at screen 1.3,1.33
set xlabel beta2label; set xrange[beta2-w_stddev*sig2:beta2+w_stddev*sig2]
set ylabel beta3label; set yrange[beta3-w_stddev*sig3:beta3+w_stddev*sig3]
set cbrange [0:1.01]
set zrange [0:1.01]
set cntrparam levels incr 0.01,0.05,1

splot "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather_beta2_beta3" u 1:2:(exp($3-lnLmax))  w l ls 1


set auto x
set auto y

set label 1 "ln L" at screen 1.3,1.33
set cbrange [lnLmax-20:lnLmax+0.1]
set zrange [lnLmax-50:lnLmax+0.1]
set cntrparam levels incr lnLmax-12,1,lnLmax


##############################
set out "SC_WS1819_3alt_altspecT_weather_lnL_beta0_beta1.png"
print "plotting SC_WS1819_3alt_altspecT_weather_lnL_beta0_beta1.png"
##############################
set xlabel beta0label; set xrange[beta0-w_stddev*sig0:beta0+w_stddev*sig0]
set ylabel beta1label; set yrange[beta1-w_stddev*sig1:beta1+w_stddev*sig1]
splot "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather_beta0_beta1" u 1:2:($3)  w l ls 1

#quit falls Mparam=2


##############################
set out "SC_WS1819_3alt_altspecT_weather_lnL_beta0_beta2.png"
print "plotting SC_WS1819_3alt_altspecT_weather_lnL_beta0_beta2.png"
##############################
set xlabel beta0label; set xrange[beta0-w_stddev*sig0:beta0+w_stddev*sig0]
set ylabel beta2label; set yrange[beta2-w_stddev*sig2:beta2+w_stddev*sig2]
splot "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather_beta0_beta2" u 1:2:($3)  w l ls 1

##############################
set out "SC_WS1819_3alt_altspecT_weather_lnL_beta1_beta2.png"
print "plotting SC_WS1819_3alt_altspecT_weather_lnL_beta1_beta2.png"
##############################
set xlabel beta1label; set xrange[beta1-w_stddev*sig1:beta1+w_stddev*sig1]
set ylabel beta2label; set yrange[beta2-w_stddev*sig2:beta2+w_stddev*sig2]
splot "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather_beta1_beta2" u 1:2:($3)  w l ls 1

#quit falls Mparam=3


##############################
set out "SC_WS1819_3alt_altspecT_weather_lnL_beta0_beta3.png"
print "plotting SC_WS1819_3alt_altspecT_weather_lnL_beta0_beta3.png"
##############################
set xlabel beta0label; set xrange[beta0-w_stddev*sig0:beta0+w_stddev*sig0]
set ylabel beta3label; set yrange[beta3-w_stddev*sig3:beta3+w_stddev*sig3]
splot "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather_beta0_beta3" u 1:2:($3)  w l ls 1

##############################
set out "SC_WS1819_3alt_altspecT_weather_lnL_beta1_beta3.png"
print "plotting SC_WS1819_3alt_altspecT_weather_lnL_beta1_beta3.png"
##############################
set xlabel beta1label; set xrange[beta1-w_stddev*sig1:beta1+w_stddev*sig1]
set ylabel beta3label; set yrange[beta3-w_stddev*sig3:beta3+w_stddev*sig3]
splot "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather_beta1_beta3" u 1:2:($3)  w l ls 1

##############################
set out "SC_WS1819_3alt_altspecT_weather_lnL_beta2_beta3.png"
print "plotting SC_WS1819_3alt_altspecT_weather_lnL_beta2_beta3.png"
##############################
set xlabel beta2label; set xrange[beta2-w_stddev*sig2:beta2+w_stddev*sig2]
set ylabel beta3label; set yrange[beta3-w_stddev*sig3:beta3+w_stddev*sig3]
splot "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather_beta2_beta3" u 1:2:($3)  w l ls 1

#quit falls Mparam=4

#quit

##############################
set out "SC_WS1819_3alt_altspecT_weather_lnL_beta0_beta4.png"
print "plotting SC_WS1819_3alt_altspecT_weather_lnL_beta0_beta4.png"
##############################
set xlabel beta0label; set xrange[beta0-w_stddev*sig0:beta0+w_stddev*sig0]
set ylabel beta4label; set yrange[beta4-w_stddev*sig4:beta4+w_stddev*sig4]
splot "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather_beta0_beta4" u 1:2:($3)  w l ls 1

##############################
set out "SC_WS1819_3alt_altspecT_weather_lnL_beta1_beta4.png"
print "plotting SC_WS1819_3alt_altspecT_weather_lnL_beta1_beta4.png"
##############################
set xlabel beta1label; set xrange[beta1-w_stddev*sig1:beta1+w_stddev*sig1]
set ylabel beta4label; set yrange[beta4-w_stddev*sig4:beta4+w_stddev*sig4]
splot "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather_beta1_beta4" u 1:2:($3)  w l ls 1

##############################
set out "SC_WS1819_3alt_altspecT_weather_lnL_beta2_beta4.png"
print "plotting SC_WS1819_3alt_altspecT_weather_lnL_beta2_beta4.png"
##############################
set xlabel beta2label; set xrange[beta2-w_stddev*sig2:beta2+w_stddev*sig2]
set ylabel beta4label; set yrange[beta4-w_stddev*sig4:beta4+w_stddev*sig4]
splot "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather_beta2_beta4" u 1:2:($3)  w l ls 1

##############################
set out "SC_WS1819_3alt_altspecT_weather_lnL_beta3_beta4.png"
print "plotting SC_WS1819_3alt_altspecT_weather_lnL_beta3_beta4.png"
##############################
set xlabel beta3label; set xrange[beta3-w_stddev*sig3:beta3+w_stddev*sig3]
set ylabel beta4label; set yrange[beta4-w_stddev*sig4:beta4+w_stddev*sig4]
splot "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather_beta3_beta4" u 1:2:($3)  w l ls 1

#quit falls Mparam=5
quit

##############################
set out "SC_WS1819_3alt_altspecT_weather_lnL_beta0_beta5.png"
print "plotting SC_WS1819_3alt_altspecT_weather_lnL_beta0_beta5.png"
##############################
set xlabel beta0label; set xrange[beta0-w_stddev*sig0:beta0+w_stddev*sig0]
set ylabel beta5label; set yrange[beta5-w_stddev*sig5:beta5+w_stddev*sig5]
splot "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather_beta0_beta5" u 1:2:($3)  w l ls 1

##############################
set out "SC_WS1819_3alt_altspecT_weather_lnL_beta1_beta5.png"
print "plotting SC_WS1819_3alt_altspecT_weather_lnL_beta1_beta5.png"
##############################
set xlabel beta1label; set xrange[beta1-w_stddev*sig1:beta1+w_stddev*sig1]
set ylabel beta5label; set yrange[beta5-w_stddev*sig5:beta5+w_stddev*sig5]
splot "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather_beta1_beta5" u 1:2:($3)  w l ls 1

##############################
set out "SC_WS1819_3alt_altspecT_weather_lnL_beta2_beta5.png"
print "plotting SC_WS1819_3alt_altspecT_weather_lnL_beta2_beta5.png"
##############################
set xlabel beta2label; set xrange[beta2-w_stddev*sig2:beta2+w_stddev*sig2]
set ylabel beta5label; set yrange[beta5-w_stddev*sig5:beta5+w_stddev*sig5]
splot "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather_beta2_beta5" u 1:2:($3)  w l ls 1

##############################
set out "SC_WS1819_3alt_altspecT_weather_lnL_beta3_beta5.png"
print "plotting SC_WS1819_3alt_altspecT_weather_lnL_beta3_beta5.png"
##############################
set xlabel beta3label; set xrange[beta3-w_stddev*sig3:beta3+w_stddev*sig3]
set ylabel beta5label; set yrange[beta5-w_stddev*sig5:beta5+w_stddev*sig5]
splot "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather_beta3_beta5" u 1:2:($3)  w l ls 1

##############################
set out "SC_WS1819_3alt_altspecT_weather_lnL_beta4_beta5.png"
print "plotting SC_WS1819_3alt_altspecT_weather_lnL_beta4_beta5.png"
##############################
set xlabel beta4label; set xrange[beta4-w_stddev*sig4:beta4+w_stddev*sig4]
set ylabel beta5label; set yrange[beta5-w_stddev*sig5:beta5+w_stddev*sig5]
splot "calibDiscrChoice/SC_WS1819_3alt_altspecT_weather_beta4_beta5" u 1:2:($3)  w l ls 1

#quit falls Mparam=6













