


# siehe ~/info/gnuplot, ~/info/gnuTemplate.gnu

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



#kind of symbols of point style ps file (aug05) 
#BUG: Keine oen circles,squares,triangles!!
#     obskure Abhaengigkeit der Symbole von Reihenfolge bzw linestyles
#--------------------------------------

#p 1 = vertical dash
#p 2 = douple diamond
#p 3 = star
#p 4 = closed square 
#p 5 = another closed square
#p 6 = closed circle
#p 7 = closed hexagon
#p 8 = upwards closed triangle
#p 9 = small upwards closed triangle
#p 10= upside-down closed triangle
#p 11= small upside-down closed triangle
#p 12= closed diamond
#p 13= small closed diamond
#p 14= closed pentagon
#p 15/16= more hexagons
#p 17= semi-open circle
#p 18=another closed circle
#p 19/20=similar to 17/18

############### Beispiele fuer Funktionen ####################

gauss(mu,sigma,x) = exp(-(x-mu)**2/(2* sigma**2)) / (sigma * sqrt(2*pi))

xlimited(x)       =(x<-5) ? -5 : ( (x>5) ? 5 : x)
phi(x)            =0.5*(erf(xlimited(x/sqrt(2)))+1) 
                  # integral der Standardnormalverteilung

lorenz(x)         = 1/(pi*(1+x**2))
expo(lambda,x)    = lambda*exp(-lambda*x)
gleich(a,b,x)     = ((x>=a)&&(x<=b)) ? 1/(b-a) : 0

studnorm(nu)  = gamma(0.5*(nu+1)) / (sqrt(nu*pi)*gamma(0.5*nu))
student(x,nu) =  studnorm(nu) / (1.+x**2/nu)**(0.5*(nu+1))

chinorm(nu)   = 1./(2**(0.5*nu)*gamma(0.5*nu))
chi2(x,nu)    =  chinorm(nu) *exp(-0.5*x) * x**(0.5*nu-1)

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x

eulerKonst=0.577216

#gamma(n+1)=n! for integer n>=0

#erlang(1,x)=standard Poisson distribution

b_erlang(n)=gamma(n+1) /gamma(n)
a_erlang(n)=(b_erlang(n))**(n)/gamma(n)
erlang(n,x)=a_erlang(n)*x**(n-1)*exp(-b_erlang(n)*x)

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
r2=3.5
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


#Log-Likelihood 

logL(b1,b2)=h11*log(P11(b1,b2))+h12*log(P12(b1,b2))\
  +h13*log(P13(b1,b2))+h14*log(P14(b1,b2))\
  +h21*log(P21(b1,b2))+h22*log(P22(b1,b2))\
  +h23*log(P23(b1,b2))+h24*log(P24(b1,b2))

#Functional matrix d(F_k)/d(beta_k')=-d^2 (log L)/(d beta_k d beta_l)

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
print "Initialwerte: b1=",b1,"  b2=",b2," F1^2+F2^2=",errorSQR(b1,b2)
print "X1data=R1=",R1," X1sim=R1sim=",R1+F1(b1,b2)," \tF1=R1sim+R1=",F1(b1,b2)
print "X2data=N1=",N1," X2sim=N1sim=",N1+F2(b1,b2)," \tF2=N1sim+N1=",F2(b1,b2)
print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," F1^2+F2^2=",errorSQR(b1,b2)
print "X1data=R1=",R1," X1sim=R1sim=",R1+F1(b1,b2)," \tF1=R1sim+R1=",F1(b1,b2)
print "X2data=N1=",N1," X2sim=N1sim=",N1+F2(b1,b2)," \tF2=N1sim+N1=",F2(b1,b2)
it=it+1

print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," F1^2+F2^2=",errorSQR(b1,b2)
print "X1data=R1=",R1," X1sim=R1sim=",R1+F1(b1,b2)," \tF1=R1sim+R1=",F1(b1,b2)
print "X2data=N1=",N1," X2sim=N1sim=",N1+F2(b1,b2)," \tF2=N1sim+N1=",F2(b1,b2)
it=it+1

print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," F1^2+F2^2=",errorSQR(b1,b2)
print "X1data=R1=",R1," X1sim=R1sim=",R1+F1(b1,b2)," \tF1=R1sim+R1=",F1(b1,b2)
print "X2data=N1=",N1," X2sim=N1sim=",N1+F2(b1,b2)," \tF2=N1sim+N1=",F2(b1,b2)
it=it+1

print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," F1^2+F2^2=",errorSQR(b1,b2)
print "X1data=R1=",R1," X1sim=R1sim=",R1+F1(b1,b2)," \tF1=R1sim+R1=",F1(b1,b2)
print "X2data=N1=",N1," X2sim=N1sim=",N1+F2(b1,b2)," \tF2=N1sim+N1=",F2(b1,b2)
it=it+1

print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," F1^2+F2^2=",errorSQR(b1,b2)
print "X1data=R1=",R1," X1sim=R1sim=",R1+F1(b1,b2)," \tF1=R1sim+R1=",F1(b1,b2)
print "X2data=N1=",N1," X2sim=N1sim=",N1+F2(b1,b2)," \tF2=N1sim+N1=",F2(b1,b2)
it=it+1


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

set xlabel "{/Symbol b}_1"
set xrange [-0.4*scalR:-0.25*scalR]
#set xrange [-2.5:-1]
set ylabel "{/Symbol b}_2"
set yrange [0.5:1.35]
set label 1 "sqrt(F_1^2+F_2^2)" at screen 1.34,1.25

errorSQRrestr(b1,b2)=min(110, errorSQR(b1,b2))
errorRMSrestr(b1,b2)=min(22., sqrt(errorSQR(b1,b2)))

#################################

set out "BNL_kalib_errorRMS.eps"
print "plotting BNL_kalib_errorRMS.eps"

splot errorRMSrestr(x,y) w l lt 6 lw 4 


##############################
# Plotte L(beta1,beta2) und ln L(beta1,beta2)
##############################

set xrange [-0.8*scalR:0*scalR]
set yrange [-0.5:2.5] 


##############################
c=exp(-3)
set out "BNL_kalib_L.eps"
print "plotting BNL_kalib_L.eps"
unset cntrparam
set cntrparam levels discrete 0.1,0.2,0.3,0.4,0.5,c-0.002,c-0.001,c,c+0.001,c+0.002,0.7,0.8,0.9    
set label 1 "Likelihood"  at screen 0.9,1.1 front 

set style line 99 lt 7 lw 1        # Farbe+Dicke fuer optionales Gitternetz (schwarz)
set pm3d hidden3d 99


set nogrid
set view 25,25

set style line 99 lt 7 lw 1        # Farbe+Dicke fuer optionales Gitternetz
#set pm3d  hidden3d 99 
#set surface
#set contour surface
set isosample 40,40
set ztics 0.2
unset pm3d; set pm3d 
logLmax=logL(b1,b2)

splot exp(logL(x,y)-logLmax) w l lt 6 lw 2 

unset surface
set isosample 60,60
unset pm3d
set pm3d
set pm3d; set pm3d map
set grid

##############################
set out "BNL_kalib_logL.eps"
print "plotting BNL_kalib_logL.eps"
##############################
set xrange [-1*scalR:0.2*scalR]
set yrange [-1:2.5] 

zmax=-26
set cbrange [-50:zmax] 
unset cntrparam
set cntrparam levels incr -40,1,zmax
set label 1 "ln (L)"  at b1-0.05,b2 front
set title "Log-Likelihood"
splot logL(x,y) w l lt 6 lw 4 

##############################
set out "BNL_kalib_logLdurchLmax.eps"
print "plotting BNL_kalib_logLdurchLmax.eps"
##############################

zmax=0
set cbrange [-25:zmax] 
set cntrparam levels incr -15,1,zmax
set label 1 "ln (L/L_{max})"  at b1-0.05,b2 front
splot logL(x,y)-logLmax w l lt 6 lw 4 



##############################
# Multipliziere alle rel. Haeuf. mit 10!
##############################

print ""
print "!! Multipliziere alle hki mit 10 und plotte L und ln L erneut!"
print ""
h11=80.  # Beobachtete absolute Hauef. Alternative k, Klasse i
h12=80.
h13=20.
h14=0.

h21=50.
h22=90.
h23=110.
h24=60.



##############################
print "plotting BNL_kalib_logL10.eps"
set out "BNL_kalib_logL10.eps"
##############################

set title "Log-Likelihood, 10 facher Stichprobenumfang"
set label 1 "ln (L/L_{max})"  at b1+0.1,b2+0.5 front
logLmax=logL(b1,b2)
splot logL(x,y)-logLmax w l lt 6 lw 4 


#########################################################
# Ergebnis plotten
#########################################################

set term post eps enhanced color solid "Helvetica" 24
set nolabel 1

print ""
print "!! Setze alle hki wieder auf Originalwerte zurueck!"

h11=8.  # Beobachtete absolute Hauef. Alternative k, Klasse i
h12=8.
h13=2.
h14=0.

h21=5.
h22=9.
h23=11.
h24=6.

print ""
print "Ergebnisse"
print "========="
print ""
print "mittlere Reiseweite <r>=R/n= ",R/n," (km)"
print "beta1 (Vorteilsaend. Fu\337+Rad pro km)=", b1/scalR, "1/km"
print "beta2 (Vorteil Fu\337+Rad bei r=0) =",b2

v1_kmh=10.
v2_kmh=30.
invlambda_min=60*(1./v2_kmh-1./v1_kmh)/b1  # 1/\lambda
stdev_min=pi/sqrt(6.)*invlambda_min
T20=b2*invlambda_min                      # beta_2/\lambda

print "Annahme v1_kmh=",v1_kmh," v2_kmh=",v2_kmh
print "Gumbel-Skalierungsfaktor: lambda=1/(",invlambda_min," min)"
print "Standardabw. des Zufallsnutzens pi/sqrt(6.)*invlambda_min=",stdev_min
print "Ruestzeit OEPNV/Kfz: T20=",T20," min"

print ""
print "Indukt. Statistik: Unsicherheiten/Konfidenzintervalle direkt durch ML-Ansatz!"
print "S_ij=Kovarianzmatrix der (b1,b2)-Unsicherheit=J^{-1}/2"
print ""
s11=0.5*invJ11(b1,b2)
s12=0.5*invJ12(b1,b2)
s21=s12
s22=0.5*invJ22(b1,b2)
print "s11=",s11,"  sqrt(s11)=",sqrt(s11)
print "s12=s21=",s12
print "s22=",s22,"  sqrt(s22)=",sqrt(s22)
print "r12=s12/sqrt(s11*s22)=",s12/sqrt(s11*s22)
print "ln(Lmax)=",logL(b1,b2)

############################
set out "BNL_result1.eps"
print "plotting BNL_result1.eps"
############################

unset colorbox
set size 1,1
set key at screen 0.83,0.6
set param
set xlabel "Reiseweite (km)"
set ylabel "Anteilswerte k=1 (Fu\337+Rad)"
xmax=20.
set xrange [0:xmax]
set yrange [0:1] reverse
plot[x=0:xmax/scalR]\
 scalR*x, P1(x,b1,b2) t "Alternative 1 (Fu\337+Rad)" w l ls 7,\
 scalR*x, P2(x,b1,b2) t "Alternative 2 (\326V+MIV)" w l ls 2,\
 scalR*r1,h11/n1 t "Daten (Fu\337+Rad)" w p ls 7,\
 scalR*r2,h12/n2 t "" w p ls 7,\
 scalR*r3,h13/n3 t "" w p ls 7,\
 scalR*r4,h14/n4 t "" w p ls 7,\
 scalR*r1,h21/n1 t "Daten (\326V+MIV)" w p ls 2,\
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

