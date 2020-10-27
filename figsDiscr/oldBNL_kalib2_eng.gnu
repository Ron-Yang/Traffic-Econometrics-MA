

min(x,y)=(x<y) ? x : y


#########################################################
# Minimalbeispiel: 2 Modi, 2 Modellparameter, Studentenbefragung 
# in VL zu Wegeentfernung und Verkehrsmittelwahl 2008
# k=1 => Fuss+Rad
# k=2 => \"OPNV und MIV
#########################################################

#########################################################
# Mittl. Entfernung Entf.-klasse l; 
# BNL_kalib2.gnu unterscheidet sich nur in anderen Klassenmittelwerten!
r1=1.   
r2=3.
r3=7.5
r4=15.
#########################################################

h11=8.  # Beobachtete absolute Hauef. Alternative k, Klasse i
h12=8.
h13=2.
h14=0.

h21=5.
h22=9.
h23=11.
h24=6.

n1=h11+h21   # Zahl aller Befragten in Entf. Klasse 1 
n2=h12+h22 
n3=h13+h23 
n4=h14+h24 

# Beobachtete Gesamtkilometerzahl in Verkehrsmodus k=1(Fuss+Rad)
R1=h11*r1+h12*r2+h13*r3+h14*r4
R2=h21*r1+h22*r2+h23*r3+h24*r4
R=R1+R2
print "R1=",R1
print "R2=",R2

# Beobachtete Gesamtzahl an Entscheidungen fuer k=1
N1=h11+h12+h13+h14
N2=h21+h22+h23+h24
n=N1+N2
print "N1=",N1
print "N2=",N2
print "n=",n
#scalR=R/n
scalR=1.

print ""
print "Skalierung auf <r>: scalR=(N1+N2)/(R1+R2)=",scalR

r1=r1/scalR
r2=r2/scalR
r3=r3/scalR
r4=r4/scalR
R1=R1/scalR
R2=R2/scalR

#V_{ki}=Vk(r_i)

V1(r,beta1,beta2)=beta1*r+beta2
V2(r,beta1,beta2)=0.

#P_{ki}=Pk(r_i)

P1(r,beta1,beta2)=1./(1.+exp(-V1(r,beta1,beta2)))
P2(r,beta1,beta2)=1./(1.+exp(V1(r,beta1,beta2)))

P11(beta1,beta2)=P1(r1,beta1,beta2)
P12(beta1,beta2)=P1(r2,beta1,beta2)
P13(beta1,beta2)=P1(r3,beta1,beta2)
P14(beta1,beta2)=P1(r4,beta1,beta2)

P21(beta1,beta2)=P2(r1,beta1,beta2)
P22(beta1,beta2)=P2(r2,beta1,beta2)
P23(beta1,beta2)=P2(r3,beta1,beta2)
P24(beta1,beta2)=P2(r4,beta1,beta2)

# Kalibrierungsbedingung 1: F1=<R(k=1)>_theo - R1=0
# Kalibrierungsbedingung 2: F2=<n(k=1)>_theo - N1=0

F1(b1,b2)=n1*r1*P11(b1,b2)+n2*r2*P12(b1,b2)\
   +n3*r3*P13(b1,b2)+n4*r4*P14(b1,b2) - R1
F2(b1,b2)=n1*P11(b1,b2)+n2*P12(b1,b2)\
   +n3*P13(b1,b2)+n4*P14(b1,b2) -N1

# PPi=P1i*P2i
PP1(b1,b2)=P11(b1,b2)*P21(b1,b2)
PP2(b1,b2)=P12(b1,b2)*P22(b1,b2)
PP3(b1,b2)=P13(b1,b2)*P23(b1,b2)
PP4(b1,b2)=P14(b1,b2)*P24(b1,b2)

#Functional matrix d(F_k)/d(beta_k')


J11(b1,b2)=n1*r1**2*PP1(b1,b2)+n2*r2**2*PP2(b1,b2)\
                  +n3*r3**2*PP3(b1,b2)+n4*r4**2*PP4(b1,b2)
J12(b1,b2)=n1*r1*PP1(b1,b2)+n2*r2*PP2(b1,b2)+n3*r3*PP3(b1,b2)+n4*r4*PP4(b1,b2)
J21(b1,b2)=J12(b1,b2)
J22(b1,b2)=n1*PP1(b1,b2)+n2*PP2(b1,b2)+n3*PP3(b1,b2)+n4*PP4(b1,b2)

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
b1=-0.
b2=0.
it=1
#############################################


print ""
print "Newton-Iteration"
print "============"
print ""
print "Initialwerte: b1=",b1,"  b2=",b2," errorSQR=",errorSQR(b1,b2)

print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," errorSQR=",errorSQR(b1,b2)
it=it+1

print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," errorSQR=",errorSQR(b1,b2)
it=it+1

print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," errorSQR=",errorSQR(b1,b2)
it=it+1

print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," errorSQR=",errorSQR(b1,b2)
it=it+1

print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," errorSQR=",errorSQR(b1,b2)
it=it+1

#########################################################
# Kalibrierung plotten
#########################################################

print ""

set term post eps enhanced color solid "Helvetica" 24
set nokey
set noparam

set isosample 60,60

set palette defined ( 0 "#dd00ff", 5 "blue", \
      15 "green", 25 "yellow", 40 "orange", 60 "red", 100 "#aa0066")

set pm3d       # bugfrei/am sichersten:  "set pm3d", DANN "set pm3d map" 
set pm3d; set pm3d map 
unset surface
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

set xlabel "{/Symbol b}_1"
set xrange [-0.4*scalR:-0.25*scalR]
#set xrange [-2.5:-1]
set ylabel "{/Symbol b}_2"
set yrange [0.5:1.35]
set label 1 "sqrt(F_1^2+F_2^2)" at screen 1.2,1.35

errorSQRrestr(b1,b2)=min(110, errorSQR(b1,b2))
errorRMSrestr(b1,b2)=min(22., sqrt(errorSQR(b1,b2)))

set out "BNL_kalib2_errorRMS.eps"
print "plotting BNL_kalib2_errorRMS.eps"

splot errorRMSrestr(x,y) w l lt 6 lw 4 

set nolabel 1

##############################

set out "BNL_kalib2_errorSQR_b1.eps"
print "plotting  BNL_kalib2_errorSQR_b1.eps"
set yrange [0:10]
set key
plot\
 errorSQRrestr(x,0.7) w l ls 8,\
 errorSQRrestr(x,0.8) w l ls 1,\
 errorSQRrestr(x,0.9) w l ls 2,\
 errorSQRrestr(x,1.0) w l ls 3,\
 errorSQRrestr(x,1.1) w l ls 4,\
 errorSQRrestr(x,1.2) w l ls 5,\
 errorSQRrestr(x,1.3) w l ls 6,\
 errorSQRrestr(x,1.4) w l ls 7


#########################################################
# Ergebnis plotten
#########################################################

print ""
print "Ergebnisse"
print "========="
print ""
print "mittlere Reiseweite scalR =<r>= ",scalR," (km)"
print "beta1 (Vorteilsaend. Fu3+Rad pro km)=", b1/scalR, "1/km"
print "beta2 (Vorteil Fu3+Rad bei r=0) =",b2
set out "BNL_result2.eps"
print "plotting BNL_result2.eps"

set size 1,1
set key at screen 0.93,0.6
set param
set xlabel "Reiseweite (km)"
set ylabel "Anteilswerte k=1 (Fu3+Rad)"
xmax=20.
set xrange [0:xmax]
set yrange [0:1]
unset colorbox
plot[x=0:xmax/scalR]\
 scalR*x, P1(x,b1,b2) t "Alternative 1 (Fu3+Rad)" w l ls 7,\
 scalR*x, P2(x,b1,b2) t "Alternative 2 (OV+MIV)" w l ls 2,\
 scalR*r1,h11/n1 t "Daten (Fu3+Rad)" w p ls 7,\
 scalR*r2,h12/n2 t "" w p ls 7,\
 scalR*r3,h13/n3 t "" w p ls 7,\
 scalR*r4,h14/n4 t "" w p ls 7,\
 scalR*r1,h21/n1 t "Daten (OV+MIV)" w p ls 2,\
 scalR*r2,h22/n2 t "" w p ls 2,\
 scalR*r3,h23/n3 t "" w p ls 2,\
 scalR*r4,h24/n4 t "" w p ls 2










quit
############################################################


set out "BNL_F1.eps"
print "plotting BNL_F1.eps"
splot F1(x,y) w l lt 6 lw 4 

set out "BNL_F2.eps"
print "plotting BNL_F2.eps"
splot F2(x,y) w l lt 6 lw 4 

