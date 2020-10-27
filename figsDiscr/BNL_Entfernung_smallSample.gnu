



#########################################################
# Minimalbeispiel: 2 Modi, 2 Modellparameter, Studentenbefragung 
# small sample for LR tests etc
# i=1 => Fuss+Rad
# i=2 => \"OPNV und MIV
#########################################################

#########################################################
# User input 
#########################################################

r1=0.5
r2=1.5   
r3=3.5
r4=7.5
r5=15.

y11=7.; y12=2.;
y21=6.; y22=4.;
y31=6.; y32=12.;
y41=1.; y42=10.;
y51=0.; y52=5.;

#########################################################
# ML calibration (simple Newton)
#########################################################




y1=y11+y12   # Zahl aller Befragten in Entf. Klasse 1 
y2=y21+y22 
y3=y31+y32 
y4=y41+y42 
y5=y51+y52 

# Beobachtete Gesamtkilometerzahl in Verkehrsmodus k=1(Fuss+Rad)
R1=y11*r1+y21*r2+y31*r3+y41*r4+y51*r5
R2=y12*r1+y22*r2+y32*r3+y42*r4+y52*r5
R=R1+R2
print ""
print "Merkmalssummen (Alt 1 relevant):\n"
print "R1=",R1
print "R2=",R2

# Beobachtete Gesamtzahl an Entscheidungen fuer k=1
N1=y11+y21+y31+y41+y51
N2=y12+y22+y32+y42+y52
N=N1+N2
print "N1=",N1
print "N2=",N2
print "N=",N

#V_{ni}=V_i(r_n)

V1(r,beta1,beta2)=beta1*r+beta2
V2(r,beta1,beta2)=0.

# logit probabilities

P1(r,beta1,beta2)=1./(1.+exp(-V1(r,beta1,beta2)))
P2(r,beta1,beta2)=1./(1.+exp(V1(r,beta1,beta2)))

#P_{ki}=Pk(r_i)

P11(beta1,beta2)=P1(r1,beta1,beta2)
P21(beta1,beta2)=P1(r2,beta1,beta2)
P31(beta1,beta2)=P1(r3,beta1,beta2)
P41(beta1,beta2)=P1(r4,beta1,beta2)
P51(beta1,beta2)=P1(r5,beta1,beta2)

P12(beta1,beta2)=P2(r1,beta1,beta2)
P22(beta1,beta2)=P2(r2,beta1,beta2)
P32(beta1,beta2)=P2(r3,beta1,beta2)
P42(beta1,beta2)=P2(r4,beta1,beta2)
P52(beta1,beta2)=P2(r5,beta1,beta2)


# Log-Likelihood

lnL(b1,b2)=\
   y11*log(P11(b1,b2))+y12*log(P12(b1,b2))\
 + y21*log(P21(b1,b2))+y22*log(P22(b1,b2))\
 + y31*log(P31(b1,b2))+y32*log(P32(b1,b2))\
 + y41*log(P41(b1,b2))+y42*log(P42(b1,b2))\
 + y51*log(P51(b1,b2))+y52*log(P52(b1,b2))



# Kalibrierungsbedingung 1: F1=<R(k=1)>_theo - R1=0
# Kalibrierungsbedingung 2: F2=<n(k=1)>_theo - N1=0

F1(b1,b2)=y1*r1*P11(b1,b2)+y2*r2*P21(b1,b2)\
   +y3*r3*P31(b1,b2)+y4*r4*P41(b1,b2)+y5*r5*P51(b1,b2) - R1
F2(b1,b2)=y1*P11(b1,b2)+y2*P21(b1,b2)\
   +y3*P31(b1,b2)+y4*P41(b1,b2)+y5*P51(b1,b2) - N1

R1mod(b1,b2)=F1(b1,b2)+R1
N1mod(b1,b2)=F2(b1,b2)+N1

# PPi=P1i*P2i
PP1(b1,b2)=P11(b1,b2)*P12(b1,b2)
PP2(b1,b2)=P21(b1,b2)*P22(b1,b2)
PP3(b1,b2)=P31(b1,b2)*P32(b1,b2)
PP4(b1,b2)=P41(b1,b2)*P42(b1,b2)
PP5(b1,b2)=P51(b1,b2)*P52(b1,b2)

#Functional matrix d(F_k)/d(beta_k')


J11(b1,b2)=y1*r1**2*PP1(b1,b2)+y2*r2**2*PP2(b1,b2)\
                  +y3*r3**2*PP3(b1,b2)+y4*r4**2*PP4(b1,b2)+y5*r5**2*PP5(b1,b2)
J12(b1,b2)=y1*r1*PP1(b1,b2)+y2*r2*PP2(b1,b2)+y3*r3*PP3(b1,b2)\
       +y4*r4*PP4(b1,b2)+y5*r5*PP5(b1,b2)
J21(b1,b2)=J12(b1,b2)
J22(b1,b2)=y1*PP1(b1,b2)+y2*PP2(b1,b2)+y3*PP3(b1,b2)+y4*PP4(b1,b2)+y5*PP5(b1,b2)

# Inverse of functional matrix

detJ(b1,b2)=J11(b1,b2)*J22(b1,b2)-J12(b1,b2)**2

invJ11(b1,b2)=J22(b1,b2)/detJ(b1,b2)
invJ12(b1,b2)=-J12(b1,b2)/detJ(b1,b2)
invJ21(b1,b2)=-J21(b1,b2)/detJ(b1,b2)
invJ22(b1,b2)=J11(b1,b2)/detJ(b1,b2)

# Newton-update step - invJ . F

dbeta1(b1,b2)= -1* (invJ11(b1,b2)*F1(b1,b2) + invJ12(b1,b2)*F2(b1,b2))
dbeta2(b1,b2)= -1* (invJ21(b1,b2)*F1(b1,b2) + invJ22(b1,b2)*F2(b1,b2))


#############################################
# Perform Newton iteration
#############################################

# Zielfunktion zum lotten u. Auswerten; 2*J =Hesse-Matrix davon
errorSQR(b1,b2)=F1(b1,b2)**2+F2(b1,b2)**2

#############################################
b1=0.
b2=0.
it=0
#############################################


print ""
print "Newton-Iteration"
print "============"
print ""
print "Initialwerte: b1=",b1,"  b2=",b2," errorSQR=",errorSQR(b1,b2)
print "R1=",R1," R1mod(b1,b2)=", R1mod(b1,b2), " grad1=F1(b1,b2)=",F1(b1,b2)
print "N1=",N1," N1mod(b1,b2)=", N1mod(b1,b2), " grad2=F2(b1,b2)=",F2(b1,b2)
it=it+1

print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," errorSQR=",errorSQR(b1,b2)
print "R1=",R1," R1mod(b1,b2)=", R1mod(b1,b2), " grad1=F1(b1,b2)=",F1(b1,b2)
print "N1=",N1," N1mod(b1,b2)=", N1mod(b1,b2), " grad2=F2(b1,b2)=",F2(b1,b2)
it=it+1

print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," errorSQR=",errorSQR(b1,b2)
print "R1=",R1," R1mod(b1,b2)=", R1mod(b1,b2), " grad1=F1(b1,b2)=",F1(b1,b2)
print "N1=",N1," N1mod(b1,b2)=", N1mod(b1,b2), " grad2=F2(b1,b2)=",F2(b1,b2)
it=it+1

print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," errorSQR=",errorSQR(b1,b2)
print "R1=",R1," R1mod(b1,b2)=", R1mod(b1,b2), " grad1=F1(b1,b2)=",F1(b1,b2)
print "N1=",N1," N1mod(b1,b2)=", N1mod(b1,b2), " grad2=F2(b1,b2)=",F2(b1,b2)
it=it+1

print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," errorSQR=",errorSQR(b1,b2)
print "R1=",R1," R1mod(b1,b2)=", R1mod(b1,b2), " grad1=F1(b1,b2)=",F1(b1,b2)
print "N1=",N1," N1mod(b1,b2)=", N1mod(b1,b2), " grad2=F2(b1,b2)=",F2(b1,b2)
it=it+1

print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," errorSQR=",errorSQR(b1,b2)
print "R1=",R1," R1mod(b1,b2)=", R1mod(b1,b2), " grad1=F1(b1,b2)=",F1(b1,b2)
print "N1=",N1," N1mod(b1,b2)=", N1mod(b1,b2), " grad2=F2(b1,b2)=",F2(b1,b2)

#######################################################

print "\nFinal calibration result after ",it," steps:"
hatbeta1=b1
hatbeta2=b2
lnLmax=lnL(hatbeta1,hatbeta2)
lnL0=lnL(0,0)

print "hatbeta1=",hatbeta1
print "hatbeta2=",hatbeta2
print "Newton error (data - model property-sum vector)^2=",errorSQR(b1,b2)

print "\nParameter Covariance matrix H^{-1}"
print "V11=invJ11(b1,b2)=",invJ11(b1,b2)
print "V12=invJ12(b1,b2)=",invJ12(b1,b2)
print "V22=invJ22(b1,b2)=",invJ22(b1,b2)

print "\nEstimation error stddev and correlations:"
print "sqrt(V11)=",sqrt(invJ11(b1,b2))
print "sqrt(V22)=",sqrt(invJ22(b1,b2))
print "rbeta12=",invJ12(b1,b2)/(sqrt(invJ11(b1,b2))*sqrt(invJ22(b1,b2)))

print "\nlnLmax=lnL(hatbeta1,hatbeta2)=",lnLmax
print "lnInit=lnL(0,0)=",lnL0



#########################################################
# plotting
#########################################################
set encoding iso_8859_1 

set style line 1 lt 7 lw 2 pt 1 ps 1.2 #schwarz, plus sign
set style line 2 lt 1 lw 4 pt 4 ps 1.2 #rot, open box
set style line 7 lt 3 lw 2 pt 10 ps 1.2 #blau, upside-down open triangle

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x

zmin=lnL0
dz=1 # distance of contour lines
set cbrange [zmin:lnLmax]
set zrange [zmin-20:lnLmax] 

print ""

set term post eps enhanced color solid "Helvetica" 28
set nokey
set noparam

set isosample 60,60

set palette defined ( -30 "#ffffff", 16 "#6666ff",\
       33 "#66ff66", 50 "yellow", 67 "orange", 84 "#ff6666",\
       96 "#ff0055",100 "#aa0000")

unset surface

set pm3d       # bugfrei/am sichersten:  "set pm3d", DANN "set pm3d map" 
set pm3d; set pm3d map 
set contour surface

set cntrparam bspline 
set cntrparam levels 30                          # n equidistant levels from min-max
set cntrparam levels discrete 0.5,1,1.5,2,4,6,8,10,12,14,16,18,20

unset clabel  # dann lauter gleiche Kontourlinien; 
                     # Farbe/Typ mit "w l ls" beim splot-Kommando

set grid   # grid/nogrid: whether 2D grid on xy plane 
set size 1.5,1.5


#################################
# eigentliches Plotten
#################################


##############################
set out "BNL_Entfernung_smallSample_logL.eps"
print "plotting BNL_Entfernung_smallSample_logL.eps"
##############################

set xlabel "{/Symbol b}_1 (Entfernungssensitivitaet Fuss/Rad gegenueber OEV/MIV)"
set xrange [-1:0]
set ylabel "{/Symbol b}_2 (AC Fuss/Rad gegenueber OEV/MIV)"
set yrange [0:2.5]
unset cntrparam
set cntrparam levels incr zmin, dz, 0  #last arg=zmax
set title "Log-Likelihood (Revealed-Choice-Befragungen smallSample)"
splot lnL(x,y) w l lt 6 lw 4 




############################
set out "BNL_Entfernung_smallSample_relHaeuf.eps"
print "plotting BNL_Entfernung_smallSample_relHaeuf.eps"
############################



unset colorbox
#set size 1,1
set key right center
set param
set xlabel "Reiseweite (km)"
set ylabel "Anteilswerte i=1 (Fu\337+Rad)"
xmax=20.
set xrange [0:xmax]
set yrange [0:1]
plot[x=0:xmax]\
 x, P1(x,hatbeta1,hatbeta2) t "Alternative 1 (Fu\337+Rad)" w l ls 7,\
 x, P2(x,hatbeta1,hatbeta2) t "Alternative 2 (\326V+MIV)" w l ls 2,\
 r1,y11/y1 t "Daten (Fu\337+Rad)" w p ls 7,\
 r2,y21/y2 t "" w p ls 7,\
 r3,y31/y3 t "" w p ls 7,\
 r4,y41/y4 t "" w p ls 7,\
 r5,y51/y5 t "" w p ls 7,\
 r1,y12/y1 t "Daten (\326V+MIV)" w p ls 2,\
 r2,y22/y2 t "" w p ls 2,\
 r3,y32/y3 t "" w p ls 2,\
 r4,y42/y4 t "" w p ls 2,\
 r5,y52/y5 t "" w p ls 2

##############################
set out "BNL_Entfernung_smallSample_logL_eng.png"
print "plotting BNL_Entfernung_smallSample_logL_eng.png"
##############################

set term pngcairo enhanced color notransparent crop font "Helvetica,12"
set size 1,1
set colorbox
set noparam
set key right top
set label 1 "ln L" at screen 0.33,0.65 front
set label 2 "distance 1 between" at screen 0.61,0.80 front
set label 3 "contour lines" at screen 0.68,0.76 front
unset title

set xlabel "Differential distance sensitivity {/Symbol b}_1"
set xrange [-1:0]
set ylabel "AC ped/bike vs. PT/car   {/Symbol b}_2"
set yrange [-1:2.5]

unset cntrparam
set cntrparam levels incr zmin, dz, 0  # last arg=zmax
splot lnL(x,y) t "" w l lt 6 lw 1

unset label 1


############################
set out "BNL_Entfernung_smallSample_relHaeuf_eng.png"
print "plotting BNL_Entfernung_smallSample_relHaeuf_eng.png"
############################

set term pngcairo enhanced color notransparent crop font "Helvetica,14"

unset colorbox
#set size 1,1
set key right center
set param
set xlabel "Distance [km]" offset 0,0.5
set ylabel "Modal split i=1 (ped/bike)" offset 1,0
xmax=20.
set xrange [0:xmax]
set yrange [0:1]
plot[x=0:xmax]\
 x, P1(x,hatbeta1,hatbeta2) t "Alternative 1 (ped/bike)" w l ls 7,\
 x, P2(x,hatbeta1,hatbeta2) t "Alternative 2 (PT/car)" w l ls 2,\
 r1,y11/y1 t "Data (ped/bike)" w p ls 7,\
 r2,y21/y2 t "" w p ls 7,\
 r3,y31/y3 t "" w p ls 7,\
 r4,y41/y4 t "" w p ls 7,\
 r5,y51/y5 t "" w p ls 7,\
 r1,y12/y1 t "Data (PT/car)" w p ls 2,\
 r2,y22/y2 t "" w p ls 2,\
 r3,y32/y3 t "" w p ls 2,\
 r4,y42/y4 t "" w p ls 2,\
 r5,y52/y5 t "" w p ls 2





