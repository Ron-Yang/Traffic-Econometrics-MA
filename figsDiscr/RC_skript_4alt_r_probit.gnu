print "hint: Update lnLmax and beta from RC_skript_4alt_r_probit.outData"

######################################################
# Input aus RC_skript_4alt_r_probit.outData
######################################################

lnLmax=-243.356
beta0=-0.924458; sig0=0.159352
beta1=-0.335206; sig1=0.0592727
beta2=-0.0738022; sig2=0.0358347
beta3=2.70827; sig3=0.395281
beta4=2.39815; sig4=0.328548
beta5=1.97387; sig5=0.292098

logitFact=pi/sqrt(6)
w_stddev=3*logitFact  # number of stddev plotted in contour plots


######################################################
# Formuliere Modell
######################################################

# 0=Fuss/Rad, 1=OEV/MIV

beta0label="{/Symbol b}_1 (Differenzielle Entfernungssensitivit\344t Fu\337-MIV)"
beta1label="{/Symbol b}_2 (Differenzielle Entfernungssensitivit\344t Rad-MIV)"
beta2label="{/Symbol b}_3 (Differenzielle Entfernungssensitivit\344t \326V-MIV)"
beta3label="{/Symbol b}_4 (AC Fu\337 - MIV)"
beta4label="{/Symbol b}_5 (AC Rad - MIV)"
beta5label="{/Symbol b}_6 (AC \326V - MIV)"


V(i,r)=\
 + r*(beta0*delta(i,0) +  beta1*delta(i,1) +  beta2*delta(i,2))\
 + beta3*delta(i,0) +  beta4*delta(i,1) +  beta5*delta(i,2)


print "AC Alt 0 vs Alt 3 in Fusskm - MIVkm = -beta3/beta0= ",-beta3/beta0
print "AC Alt 1 vs Alt 3 in Radkm - MIVkm = -beta4/beta1= ",-beta4/beta1
print "Diff. Entfernung Fuss-MIV [km] in sigeps = -1/beta0= ",-1/beta0
print "Diff. Entfernung Rad-MIV [km] in sigeps = -1/beta1= ",-1/beta1


# Ende User-Iput
######################################################
 

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x
delta (i,j)=(i==j) ? 1 : 0


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






#set term post eps enhanced color solid "Helvetica" 24
set term pngcairo enhanced color notransparent crop font "Helvetica, 14"

###############################
set out "RC_skript_4alt_r_probit_fProb.png"
print "plotting RC_skript_4alt_r_probit_fProb.png"
###############################

set nolabel
set auto
set xlabel "Personengruppe"
set ylabel "Relative H\344ufigkeiten"

set key top

set yrange [0:1.3]
set xrange [0.9:]
plot\
  "calibDiscrChoice/RC_skript_4alt_r_probit.outData" u ($1):($2/$10) t "Daten Fu\337" w p ls 7,\
  "calibDiscrChoice/RC_skript_4alt_r_probit.outData" u ($1):($3/$10) t "Daten Rad" w p ls 5,\
  "calibDiscrChoice/RC_skript_4alt_r_probit.outData" u ($1):($4/$10) t "Daten \326V" w p ls 2,\
  "calibDiscrChoice/RC_skript_4alt_r_probit.outData" u ($1):($6/$10) t "Modell Fu\337" w l ls 7,\
  "calibDiscrChoice/RC_skript_4alt_r_probit.outData" u ($1):($7/$10) t "Modell Rad" w l ls 5,\
  "calibDiscrChoice/RC_skript_4alt_r_probit.outData" u ($1):($8/$10) t "Modell \326V" w l ls 2




########################################################
# Contour plotting
########################################################



set nokey
set noparam
set grid   # grid/nogrid: whether 2D grid on xy plane 
set size 1.5,1.5
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
set out "RC_skript_4alt_r_probit_L_beta0_beta1.png"
print "plotting RC_skript_4alt_r_probit_L_beta0_beta1.png"
##############################

set label 1 "L/L_{max}" at screen 1.3,1.33
set xlabel beta0label; set xrange[logitFact*beta0-w_stddev*sig0:logitFact*beta0+w_stddev*sig0]
set ylabel beta1label; set yrange[logitFact*beta1-w_stddev*sig1:logitFact*beta1+w_stddev*sig1]
set cbrange [0:1.01]
set zrange [0:1.01]
set cntrparam levels incr 0.01,0.05,1

splot "calibDiscrChoice/RC_skript_4alt_r_probit_beta0_beta1" u (logitFact*$1):(logitFact*$2):(exp($3-lnLmax))  w l ls 1


set auto x
set auto y

set label 1 "ln L" at screen 1.3,1.33
set cbrange [lnLmax-20:lnLmax+0.1]
set zrange [lnLmax-50:lnLmax+0.1]
set cntrparam levels incr int(lnLmax-12),1,lnLmax


##############################
set out "RC_skript_4alt_r_probit_lnL_beta0_beta1.png"
print "plotting RC_skript_4alt_r_probit_lnL_beta0_beta1.png"
##############################
set xlabel beta0label; set xrange[logitFact*beta0-w_stddev*sig0:logitFact*beta0+w_stddev*sig0]
set ylabel beta1label; set yrange[logitFact*beta1-w_stddev*sig1:logitFact*beta1+w_stddev*sig1]
splot "calibDiscrChoice/RC_skript_4alt_r_probit_beta0_beta1" u (logitFact*$1):(logitFact*$2):($3)  w l ls 1

#quit falls Mparam=2



##############################
set out "RC_skript_4alt_r_probit_lnL_beta0_beta2.png"
print "plotting RC_skript_4alt_r_probit_lnL_beta0_beta2.png"
##############################
set xlabel beta0label; set xrange[logitFact*beta0-w_stddev*sig0:logitFact*beta0+w_stddev*sig0]
set ylabel beta2label; set yrange[logitFact*beta2-w_stddev*sig2:logitFact*beta2+w_stddev*sig2]
splot "calibDiscrChoice/RC_skript_4alt_r_probit_beta0_beta2" u (logitFact*$1):(logitFact*$2):($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_probit_lnL_beta1_beta2.png"
print "plotting RC_skript_4alt_r_probit_lnL_beta1_beta2.png"
##############################
set xlabel beta1label; set xrange[logitFact*beta1-w_stddev*sig1:logitFact*beta1+w_stddev*sig1]
set ylabel beta2label; set yrange[logitFact*beta2-w_stddev*sig2:logitFact*beta2+w_stddev*sig2]
splot "calibDiscrChoice/RC_skript_4alt_r_probit_beta1_beta2" u (logitFact*$1):(logitFact*$2):($3)  w l ls 1

#quit falls Mparam=3


##############################
set out "RC_skript_4alt_r_probit_lnL_beta0_beta3.png"
print "plotting RC_skript_4alt_r_probit_lnL_beta0_beta3.png"
##############################
set xlabel beta0label; set xrange[logitFact*beta0-w_stddev*sig0:logitFact*beta0+w_stddev*sig0]
set ylabel beta3label; set yrange[logitFact*beta3-w_stddev*sig3:logitFact*beta3+w_stddev*sig3]
splot "calibDiscrChoice/RC_skript_4alt_r_probit_beta0_beta3" u (logitFact*$1):(logitFact*$2):($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_probit_lnL_beta1_beta3.png"
print "plotting RC_skript_4alt_r_probit_lnL_beta1_beta3.png"
##############################
set xlabel beta1label; set xrange[logitFact*beta1-w_stddev*sig1:logitFact*beta1+w_stddev*sig1]
set ylabel beta3label; set yrange[logitFact*beta3-w_stddev*sig3:logitFact*beta3+w_stddev*sig3]
splot "calibDiscrChoice/RC_skript_4alt_r_probit_beta1_beta3" u (logitFact*$1):(logitFact*$2):($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_probit_lnL_beta2_beta3.png"
print "plotting RC_skript_4alt_r_probit_lnL_beta2_beta3.png"
##############################
set xlabel beta2label; set xrange[logitFact*beta2-w_stddev*sig2:logitFact*beta2+w_stddev*sig2]
set ylabel beta3label; set yrange[logitFact*beta3-w_stddev*sig3:logitFact*beta3+w_stddev*sig3]
splot "calibDiscrChoice/RC_skript_4alt_r_probit_beta2_beta3" u (logitFact*$1):(logitFact*$2):($3)  w l ls 1

#quit falls Mparam=4



##############################
set out "RC_skript_4alt_r_probit_lnL_beta0_beta4.png"
print "plotting RC_skript_4alt_r_probit_lnL_beta0_beta4.png"
##############################
set xlabel beta0label; set xrange[logitFact*beta0-w_stddev*sig0:logitFact*beta0+w_stddev*sig0]
set ylabel beta4label; set yrange[logitFact*beta4-w_stddev*sig4:logitFact*beta4+w_stddev*sig4]
splot "calibDiscrChoice/RC_skript_4alt_r_probit_beta0_beta4" u (logitFact*$1):(logitFact*$2):($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_probit_lnL_beta1_beta4.png"
print "plotting RC_skript_4alt_r_probit_lnL_beta1_beta4.png"
##############################
set xlabel beta1label; set xrange[logitFact*beta1-w_stddev*sig1:logitFact*beta1+w_stddev*sig1]
set ylabel beta4label; set yrange[logitFact*beta4-w_stddev*sig4:logitFact*beta4+w_stddev*sig4]
splot "calibDiscrChoice/RC_skript_4alt_r_probit_beta1_beta4" u (logitFact*$1):(logitFact*$2):($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_probit_lnL_beta2_beta4.png"
print "plotting RC_skript_4alt_r_probit_lnL_beta2_beta4.png"
##############################
set xlabel beta2label; set xrange[logitFact*beta2-w_stddev*sig2:logitFact*beta2+w_stddev*sig2]
set ylabel beta4label; set yrange[logitFact*beta4-w_stddev*sig4:logitFact*beta4+w_stddev*sig4]
splot "calibDiscrChoice/RC_skript_4alt_r_probit_beta2_beta4" u (logitFact*$1):(logitFact*$2):($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_probit_lnL_beta3_beta4.png"
print "plotting RC_skript_4alt_r_probit_lnL_beta3_beta4.png"
##############################
set xlabel beta3label; set xrange[logitFact*beta3-w_stddev*sig3:logitFact*beta3+w_stddev*sig3]
set ylabel beta4label; set yrange[logitFact*beta4-w_stddev*sig4:logitFact*beta4+w_stddev*sig4]
splot "calibDiscrChoice/RC_skript_4alt_r_probit_beta3_beta4" u (logitFact*$1):(logitFact*$2):($3)  w l ls 1

#quit falls Mparam=5


##############################
set out "RC_skript_4alt_r_probit_lnL_beta0_beta5.png"
print "plotting RC_skript_4alt_r_probit_lnL_beta0_beta5.png"
##############################
set xlabel beta0label; set xrange[logitFact*beta0-w_stddev*sig0:logitFact*beta0+w_stddev*sig0]
set ylabel beta5label; set yrange[logitFact*beta5-w_stddev*sig5:logitFact*beta5+w_stddev*sig5]
splot "calibDiscrChoice/RC_skript_4alt_r_probit_beta0_beta5" u (logitFact*$1):(logitFact*$2):($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_probit_lnL_beta1_beta5.png"
print "plotting RC_skript_4alt_r_probit_lnL_beta1_beta5.png"
##############################
set xlabel beta1label; set xrange[logitFact*beta1-w_stddev*sig1:logitFact*beta1+w_stddev*sig1]
set ylabel beta5label; set yrange[logitFact*beta5-w_stddev*sig5:logitFact*beta5+w_stddev*sig5]
splot "calibDiscrChoice/RC_skript_4alt_r_probit_beta1_beta5" u (logitFact*$1):(logitFact*$2):($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_probit_lnL_beta2_beta5.png"
print "plotting RC_skript_4alt_r_probit_lnL_beta2_beta5.png"
##############################
set xlabel beta2label; set xrange[logitFact*beta2-w_stddev*sig2:logitFact*beta2+w_stddev*sig2]
set ylabel beta5label; set yrange[logitFact*beta5-w_stddev*sig5:logitFact*beta5+w_stddev*sig5]
splot "calibDiscrChoice/RC_skript_4alt_r_probit_beta2_beta5" u (logitFact*$1):(logitFact*$2):($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_probit_lnL_beta3_beta5.png"
print "plotting RC_skript_4alt_r_probit_lnL_beta3_beta5.png"
##############################
set xlabel beta3label; set xrange[logitFact*beta3-w_stddev*sig3:logitFact*beta3+w_stddev*sig3]
set ylabel beta5label; set yrange[logitFact*beta5-w_stddev*sig5:logitFact*beta5+w_stddev*sig5]
splot "calibDiscrChoice/RC_skript_4alt_r_probit_beta3_beta5" u (logitFact*$1):(logitFact*$2):($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_probit_lnL_beta4_beta5.png"
print "plotting RC_skript_4alt_r_probit_lnL_beta4_beta5.png"
##############################
set xlabel beta4label; set xrange[logitFact*beta4-w_stddev*sig4:logitFact*beta4+w_stddev*sig4]
set ylabel beta5label; set yrange[logitFact*beta5-w_stddev*sig5:logitFact*beta5+w_stddev*sig5]
splot "calibDiscrChoice/RC_skript_4alt_r_probit_beta4_beta5" u (logitFact*$1):(logitFact*$2):($3)  w l ls 1

print "get RC_skript_4alt_r_probit_corrMatrix.png from ./calibDiscrChoice/"










