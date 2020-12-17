print "hint: Update lnLmax and beta from RC_skript_4alt_r_eng.outData"

######################################################
# Input aus RC_skript_4alt_r_eng.outData
######################################################

lnLinit=-368.754
lnLmax=-243.318
N=266        #(number of single choices)
Mparam=6
beta0=4.0989; sig0=0.598725
beta1=3.55734; sig1=0.513962
beta2=2.95378; sig2=0.467502
beta3=-1.42553; sig3=0.261516
beta4=-0.484265; sig4=0.0839078
beta5=-0.137775; sig5=0.0488399

w_stddev=3  # number of stddev plotted in contour plots


######################################################
# Formuliere Modell
######################################################

# 0=Ped, 1=bike, 2=PT,3=car

beta0label="{/Symbol b}_1 (AC ped - car)"
beta1label="{/Symbol b}_2 (AC bike - car)"
beta2label="{/Symbol b}_3 (AC PT - car)"
beta3label="{/Symbol b}_4 (Distance sensitivity diff. ped-car)"
beta4label="{/Symbol b}_5 (Distance sensitivity diff. bike-car)"
beta5label="{/Symbol b}_6 (Distance sensitivity diff. PT-car)"


V(i,r)=\
 + r*(beta3*delta(i,0) +  beta4*delta(i,1) +  beta5*delta(i,2))\
 + beta0*delta(i,0) +  beta1*delta(i,1) +  beta2*delta(i,2)

denom(r)=exp(V(0,r))+exp(V(1,r))+exp(V(2,r))+exp(V(3,r))

Prob(i,r)=exp(V(i,r))/denom(r)

logitFact=pi/sqrt(6)

print "AC Alt 0 vs Alt 3 in Ped-km - car-km = -beta0/beta3= ",-beta0/beta3
print "AC Alt 1 vs Alt 3 in Rad-km - car-km = -beta1/beta4= ",-beta1/beta4
print "Diff. Entfernung Ped-car [km] per NE= -1/beta3= ",-1/beta3
print "Diff. Entfernung Rad-car [km] per NE= -1/beta4= ",-1/beta4


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



set style line 11 lt 1 lw 3 pt 7 ps 1.9  lc rgb "#000000" #schwarz,solid,bullet
set style line 12 lt 1 lw 3 pt 2 ps 1.5  lc rgb "#CC0022" #rot, dash, Kreuz
set style line 13 lt 8 lw 3 pt 3 ps 1.2 #blassrot, offener Stern
set style line 14 lt 6 lw 3 pt 4 ps 1.5  lc rgb "#FFAA00"  #gelb, offenes Quadrat
set style line 15 lt 1 lw 3 pt 5 ps 1.5  lc rgb "#00DD22"  #gruen,solid,ClosedBox
set style line 16 lt 5 lw 3 pt 7 ps 1.5  lc rgb "#00CCCC" #offener Kreis
set style line 17 lt 1 lw 3 pt 7 ps 1.5  lc rgb "#1100AA"  #blau,solid,Bullet
set style line 18 lt 4 lw 3 pt 8 ps 1.5 #lila, aufrechtes geschloss. Dreieck
set style line 19 lt 7 lw 3 pt 9 ps 1.5  lc rgb "#999999"  #grau, aufr. gschl. Dreieck

#Sinnvolle point typen (pt)
# 1=Plus,2=Kreuz,4=openQuare,5=closedSquare, 6=openCirc,7=closedCirc,
# 9-11=triangles, 12-13=Rauten





set term pngcairo enhanced color notransparent crop font "Helvetica, 14"


###############################
set out "RC_skript_4alt_r_eng_fProb.png"
print "plotting RC_skript_4alt_r_eng_fProb.png"
###############################

set nolabel
set auto
set xlabel "Person group" offset 0,0.5
set ylabel "Relative frequency" offset 1,0

set key top

set yrange [0:1.3]
set xrange [0.9:]
plot\
  "calibDiscrChoice/RC_skript_4alt_r_eng.outData" u ($1):($2/$10) t "Data ped" w p ls 7,\
  "calibDiscrChoice/RC_skript_4alt_r_eng.outData" u ($1):($3/$10) t "Data Rad" w p ls 5,\
  "calibDiscrChoice/RC_skript_4alt_r_eng.outData" u ($1):($4/$10) t "Data PT" w p ls 2,\
  "calibDiscrChoice/RC_skript_4alt_r_eng.outData" u ($1):($6/$10) t "Model ped" w l ls 7,\
  "calibDiscrChoice/RC_skript_4alt_r_eng.outData" u ($1):($7/$10) t "Model Rad" w l ls 5,\
  "calibDiscrChoice/RC_skript_4alt_r_eng.outData" u ($1):($8/$10) t "Model PT" w l ls 2



###############################
set out "RC_skript_4alt_r_eng_fProbKumDist.png"
print "plotting RC_skript_4alt_r_eng_fProbKumDist.png"
###############################

rmin=0.
rmax=16.
rfunc(n)=(n==1) ? 0.5:(n==2) ? 1.5 : (n==3) ? 3.5 : (n==4) ? 7.5 :15

set xlabel "Distance [km]"

set yrange [0:1.4]

plot[r=rmin:rmax]\
  r, Prob(0,r) t "Ped" w l ls 17,\
  r, Prob(0,r)+Prob(1,r) t "Ped+Rad" w l ls 15,\
  r, 1-Prob(3,r) t "Ped+Rad+OEV" w l ls 12,\
  "calibDiscrChoice/RC_skript_4alt_r_eng.outData" u (rfunc($1)):($2/$10) t "Data ped" w p ls 7,\
  "calibDiscrChoice/RC_skript_4alt_r_eng.outData" u (rfunc($1)):(($2+$3)/$10) t "" w p ls 5,\
  "calibDiscrChoice/RC_skript_4alt_r_eng.outData" u (rfunc($1)):(($2+$3+$4)/$10) t "" w p ls 2,\
  r,1 t "" w l ls 11


########################################################
# Contour plotting
########################################################


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
set out "RC_skript_4alt_r_eng_L_beta0_beta1.png"
print "plotting RC_skript_4alt_r_eng_L_beta0_beta1.png"
##############################

set label 1 "L/L_{max}" at screen 1.3,1.33
set xlabel beta0label offset 0,0.3; 
set xrange[beta0-w_stddev*sig0:beta0+w_stddev*sig0]
set ylabel beta1label offset 0.2,0; 
set yrange[beta1-w_stddev*sig1:beta1+w_stddev*sig1]
set cbrange [0:1.01]
set zrange [0:1.01]
set cntrparam levels incr 0.01,0.05,1

splot "calibDiscrChoice/RC_skript_4alt_r_eng_beta0_beta1" u 1:2:(exp($3-lnLmax))  w l ls 1


set auto x
set auto y

set label 1 "ln L" at screen 1.3,1.33
set cbrange [lnLmax-20:lnLmax+0.1]
set zrange [lnLmax-50:lnLmax+0.1]
set cntrparam levels incr int(lnLmax-12),1,lnLmax


##############################
set out "RC_skript_4alt_r_eng_lnL_beta0_beta1.png"
print "plotting RC_skript_4alt_r_eng_lnL_beta0_beta1.png"
##############################
set xlabel beta0label; set xrange[beta0-w_stddev*sig0:beta0+w_stddev*sig0]
set ylabel beta1label; set yrange[beta1-w_stddev*sig1:beta1+w_stddev*sig1]
splot "calibDiscrChoice/RC_skript_4alt_r_eng_beta0_beta1" u 1:2:($3)  w l ls 1

#quit falls Mparam=2



##############################
set out "RC_skript_4alt_r_eng_lnL_beta0_beta2.png"
print "plotting RC_skript_4alt_r_eng_lnL_beta0_beta2.png"
##############################
set xlabel beta0label; set xrange[beta0-w_stddev*sig0:beta0+w_stddev*sig0]
set ylabel beta2label; set yrange[beta2-w_stddev*sig2:beta2+w_stddev*sig2]
splot "calibDiscrChoice/RC_skript_4alt_r_eng_beta0_beta2" u 1:2:($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_eng_lnL_beta1_beta2.png"
print "plotting RC_skript_4alt_r_eng_lnL_beta1_beta2.png"
##############################
set xlabel beta1label; set xrange[beta1-w_stddev*sig1:beta1+w_stddev*sig1]
set ylabel beta2label; set yrange[beta2-w_stddev*sig2:beta2+w_stddev*sig2]
splot "calibDiscrChoice/RC_skript_4alt_r_eng_beta1_beta2" u 1:2:($3)  w l ls 1

#quit falls Mparam=3


##############################
set out "RC_skript_4alt_r_eng_lnL_beta0_beta3.png"
print "plotting RC_skript_4alt_r_eng_lnL_beta0_beta3.png"
##############################
set xlabel beta0label; set xrange[beta0-w_stddev*sig0:beta0+w_stddev*sig0]
set ylabel beta3label; set yrange[beta3-w_stddev*sig3:beta3+w_stddev*sig3]
splot "calibDiscrChoice/RC_skript_4alt_r_eng_beta0_beta3" u 1:2:($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_eng_lnL_beta1_beta3.png"
print "plotting RC_skript_4alt_r_eng_lnL_beta1_beta3.png"
##############################
set xlabel beta1label; set xrange[beta1-w_stddev*sig1:beta1+w_stddev*sig1]
set ylabel beta3label; set yrange[beta3-w_stddev*sig3:beta3+w_stddev*sig3]
splot "calibDiscrChoice/RC_skript_4alt_r_eng_beta1_beta3" u 1:2:($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_eng_lnL_beta2_beta3.png"
print "plotting RC_skript_4alt_r_eng_lnL_beta2_beta3.png"
##############################
set xlabel beta2label; set xrange[beta2-w_stddev*sig2:beta2+w_stddev*sig2]
set ylabel beta3label; set yrange[beta3-w_stddev*sig3:beta3+w_stddev*sig3]
splot "calibDiscrChoice/RC_skript_4alt_r_eng_beta2_beta3" u 1:2:($3)  w l ls 1

#quit falls Mparam=4



##############################
set out "RC_skript_4alt_r_eng_lnL_beta0_beta4.png"
print "plotting RC_skript_4alt_r_eng_lnL_beta0_beta4.png"
##############################
set xlabel beta0label; set xrange[beta0-w_stddev*sig0:beta0+w_stddev*sig0]
set ylabel beta4label; set yrange[beta4-w_stddev*sig4:beta4+w_stddev*sig4]
splot "calibDiscrChoice/RC_skript_4alt_r_eng_beta0_beta4" u 1:2:($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_eng_lnL_beta1_beta4.png"
print "plotting RC_skript_4alt_r_eng_lnL_beta1_beta4.png"
##############################
set xlabel beta1label; set xrange[beta1-w_stddev*sig1:beta1+w_stddev*sig1]
set ylabel beta4label; set yrange[beta4-w_stddev*sig4:beta4+w_stddev*sig4]
splot "calibDiscrChoice/RC_skript_4alt_r_eng_beta1_beta4" u 1:2:($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_eng_lnL_beta2_beta4.png"
print "plotting RC_skript_4alt_r_eng_lnL_beta2_beta4.png"
##############################
set xlabel beta2label; set xrange[beta2-w_stddev*sig2:beta2+w_stddev*sig2]
set ylabel beta4label; set yrange[beta4-w_stddev*sig4:beta4+w_stddev*sig4]
splot "calibDiscrChoice/RC_skript_4alt_r_eng_beta2_beta4" u 1:2:($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_eng_lnL_beta3_beta4.png"
print "plotting RC_skript_4alt_r_eng_lnL_beta3_beta4.png"
##############################
set xlabel beta3label; set xrange[beta3-w_stddev*sig3:beta3+w_stddev*sig3]
set ylabel beta4label; set yrange[beta4-w_stddev*sig4:beta4+w_stddev*sig4]
splot "calibDiscrChoice/RC_skript_4alt_r_eng_beta3_beta4" u 1:2:($3)  w l ls 1

#quit falls Mparam=5


##############################
set out "RC_skript_4alt_r_eng_lnL_beta0_beta5.png"
print "plotting RC_skript_4alt_r_eng_lnL_beta0_beta5.png"
##############################
set xlabel beta0label; set xrange[beta0-w_stddev*sig0:beta0+w_stddev*sig0]
set ylabel beta5label; set yrange[beta5-w_stddev*sig5:beta5+w_stddev*sig5]
splot "calibDiscrChoice/RC_skript_4alt_r_eng_beta0_beta5" u 1:2:($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_eng_lnL_beta1_beta5.png"
print "plotting RC_skript_4alt_r_eng_lnL_beta1_beta5.png"
##############################
set xlabel beta1label; set xrange[beta1-w_stddev*sig1:beta1+w_stddev*sig1]
set ylabel beta5label; set yrange[beta5-w_stddev*sig5:beta5+w_stddev*sig5]
splot "calibDiscrChoice/RC_skript_4alt_r_eng_beta1_beta5" u 1:2:($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_eng_lnL_beta2_beta5.png"
print "plotting RC_skript_4alt_r_eng_lnL_beta2_beta5.png"
##############################
set xlabel beta2label; set xrange[beta2-w_stddev*sig2:beta2+w_stddev*sig2]
set ylabel beta5label; set yrange[beta5-w_stddev*sig5:beta5+w_stddev*sig5]
splot "calibDiscrChoice/RC_skript_4alt_r_eng_beta2_beta5" u 1:2:($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_eng_lnL_beta3_beta5.png"
print "plotting RC_skript_4alt_r_eng_lnL_beta3_beta5.png"
##############################
set xlabel beta3label; set xrange[beta3-w_stddev*sig3:beta3+w_stddev*sig3]
set ylabel beta5label; set yrange[beta5-w_stddev*sig5:beta5+w_stddev*sig5]
splot "calibDiscrChoice/RC_skript_4alt_r_eng_beta3_beta5" u 1:2:($3)  w l ls 1

##############################
set out "RC_skript_4alt_r_eng_lnL_beta4_beta5.png"
print "plotting RC_skript_4alt_r_eng_lnL_beta4_beta5.png"
##############################
set xlabel beta4label; set xrange[beta4-w_stddev*sig4:beta4+w_stddev*sig4]
set ylabel beta5label; set yrange[beta5-w_stddev*sig5:beta5+w_stddev*sig5]
splot "calibDiscrChoice/RC_skript_4alt_r_eng_beta4_beta5" u 1:2:($3)  w l ls 1

print "get RC_skript_4alt_r_eng_corrMatrix.png from ./calibDiscrChoice/"














