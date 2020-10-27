



#########################################################
# Minimalbeispiel: 2 Modi, 2 Modellparameter, Studentenbefragung 
# in VL zu Wegeentfernung und Verkehrsmittelwahl 2007-2011
# k=1 => Fuss+Rad
# k=2 => \"OPNV und MIV
#########################################################

set encoding iso_8859_1 

set style line 1 lt 7 lw 6 pt 1 ps 1.2 #schwarz, plus sign
set style line 2 lt 1 lw 8 pt 4 ps 1.2 #rot, open box
set style line 7 lt 3 lw 4 pt 10 ps 1.2 #blau, upside-down open triangle

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x



#########################################################
# Mittl. Entfernung Entf.-klasse n; 

r1=0.5
r2=1.5   
r3=3.5
r4=7.5
r5=15.

#########################################################
# Abs. Haeufigkeiten y_{ni}

y11=33.; y12=10.;
y21=32.; y22=22.;
y31=29.; y32=59.;
y41=7.;  y42=49.;
y51=0.;  y52=25.;


y1=y11+y12   # Zahl aller Befragten in Entf. Klasse 1 
y2=y21+y22 
y3=y31+y32 
y4=y41+y42 
y5=y51+y52 

#V_{ni}=V_i(r_n)

V1(r,beta1,beta2)=beta1*r+beta2
V2(r,beta1,beta2)=0.

#P_{ki}=Pk(r_i)

P1(r,beta1,beta2)=1./(1.+exp(-V1(r,beta1,beta2)))
P2(r,beta1,beta2)=1./(1.+exp(V1(r,beta1,beta2)))

logL(b1,b2)=y11*log(P1(r1,b1,b2))+y12*log(P2(r1,b1,b2))\
 + y21*log(P1(r2,b1,b2))+y22*log(P2(r2,b1,b2))\
 + y31*log(P1(r3,b1,b2))+y32*log(P2(r3,b1,b2))\
 + y41*log(P1(r4,b1,b2))+y42*log(P2(r4,b1,b2))\
 + y51*log(P1(r5,b1,b2))+y52*log(P2(r5,b1,b2))



########################################################
print ""
print "Newton-Iteration, Start"
print "======================="
print ""

beta1=0.
beta2=0.
print "beta1=",beta1,"  beta2=",beta2

P11=P1(r1,beta1,beta2)
P21=P1(r2,beta1,beta2)
P31=P1(r3,beta1,beta2)
P41=P1(r4,beta1,beta2)
P51=P1(r5,beta1,beta2)


# Beobachtete und modellierte Gesamtkilometerzahl in Verkehrsmodus k=1(Fuss+Rad)
R1=y11*r1+y21*r2+y31*r3+y41*r4+y51*r5
R1mod=y1*P11*r1+y2*P21*r2+y3*P31*r3+y4*P41*r4+y5*P51*r5


print "R1=",R1
print "R1mod=",R1mod

# Beobachtete Gesamtzahl an Entscheidungen fuer k=1
N1=y11+y21+y31+y41+y51
N1mod=y1*P11+y2*P21+y3*P31+y4*P41+y5*P51

print "N1=",N1
print "N1mod=",N1mod



#########################################################
# Kalibrierung plotten
#########################################################

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
set out "BNL_Entfernung_2007-2011_logL.eps"
print "plotting BNL_Entfernung_2007-2011_logL.eps"
##############################

set xlabel "{/Symbol b}_1 (Entfernungssensitivitaet Fuss/Rad gegenueber OEV/MIV)"
set xrange [-1:0]
set ylabel "{/Symbol b}_2 (AC Fuss/Rad gegenueber OEV/MIV)"
set yrange [0:2.5]
zmax=-26
set cbrange [-160:]
set zrange [-180:] 
unset cntrparam
set cntrparam levels incr -150,1,zmax
set title "Log-Likelihood (Revealed-Choice-Befragungen 2007-2011)"
splot logL(x,y) w l lt 6 lw 4 




############################
set out "BNL_Entfernung_2007-2011_relHaeuf.eps"
print "plotting BNL_Entfernung_2007-2011_relHaeuf.eps"
############################

hatbeta1=-0.456565
hatbeta2=1.09616
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


############################
print "Merkmalssummen nach Kalibrierung:"
############################

P11=P1(r1,hatbeta1,hatbeta2)
P21=P1(r2,hatbeta1,hatbeta2)
P31=P1(r3,hatbeta1,hatbeta2)
P41=P1(r4,hatbeta1,hatbeta2)
P51=P1(r5,hatbeta1,hatbeta2)


# Beobachtete und modellierte Gesamtkilometerzahl in Verkehrsmodus k=1(Fuss+Rad)
R1mod=y1*P11*r1+y2*P21*r2+y3*P31*r3+y4*P41*r4+y5*P51*r5
print "R1=",R1
print "R1mod=",R1mod

# Beobachtete Gesamtzahl an Entscheidungen fuer k=1
N1=y11+y21+y31+y41+y51
N1mod=y1*P11+y2*P21+y3*P31+y4*P41+y5*P51
print "N1=",N1
print "N1mod=",N1mod



