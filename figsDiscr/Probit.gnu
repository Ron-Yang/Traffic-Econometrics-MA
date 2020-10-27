


# siehe ~/info/gnuplot, ~/info/gnuTemplate.gnu

##########################################################
#geordnet nach hue (ps)
# set style line: linetype lt, point type pt

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

set style line 99 lt 1 lw 1 pt 11 ps 0.8 lc rgb "#000000"


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
gaussNorm(x)=exp(-0.5*x**2) /sqrt(2*pi)
xlimited(x)       =(x<-5) ? -5 : ( (x>5) ? 5 : x)
phi(x)            =0.5*(erf(xlimited(x/sqrt(2)))+1) 
                  # integral der Standardnormalverteilung
pi(x)=norm(x)
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


sqrt2=sqrt(2.)
P1binom(dV)=phi(dV/sqrt2)
P2binom(dV)=1.-phi(dV/sqrt2)

set key
set term post eps enhanced color solid "Helvetica" 16
set out "test.eps"
print "plotting test.eps"
plot[x=-3:3]\
 x, P1binom(x) t "Probit" w l ls 1,\
 x, 1./(1+exp(-x)) t "Logit" w l ls 2
unset key


#########################################################
# Probit-Wahrsch. bei drei Alternativen als f(V1, V2) ; V3=0 und Var(eps)=1
#########################################################

dx=0.2
integrand1(x,V1,V2)=gaussNorm(x)*phi(x+V1-V2)*phi(x+V1)

P1trinom(V1,V2)=dx*( \
  integrand1(0,V1,V2) \
 +  integrand1(1*dx,V1,V2)+ integrand1(-1*dx,V1,V2)\
 +  integrand1(2*dx,V1,V2)+ integrand1(-2*dx,V1,V2)\
 +  integrand1(3*dx,V1,V2)+ integrand1(-3*dx,V1,V2)\
 +  integrand1(4*dx,V1,V2)+ integrand1(-4*dx,V1,V2)\
 +  integrand1(5*dx,V1,V2)+ integrand1(-5*dx,V1,V2)\
 +  integrand1(6*dx,V1,V2)+ integrand1(-6*dx,V1,V2)\
 +  integrand1(7*dx,V1,V2)+ integrand1(-7*dx,V1,V2)\
 +  integrand1(8*dx,V1,V2)+ integrand1(-8*dx,V1,V2)\
 +  integrand1(9*dx,V1,V2)+ integrand1(-9*dx,V1,V2)\
 +  integrand1(10*dx,V1,V2)+ integrand1(-10*dx,V1,V2)\
 +  integrand1(11*dx,V1,V2)+ integrand1(-11*dx,V1,V2)\
 +  integrand1(12*dx,V1,V2)+ integrand1(-12*dx,V1,V2)\
 +  integrand1(13*dx,V1,V2)+ integrand1(-13*dx,V1,V2)\
 +  integrand1(14*dx,V1,V2)+ integrand1(-14*dx,V1,V2)\
 +  integrand1(15*dx,V1,V2)+ integrand1(-15*dx,V1,V2)\
 +  integrand1(16*dx,V1,V2)+ integrand1(-16*dx,V1,V2)\
 +  integrand1(17*dx,V1,V2)+ integrand1(-17*dx,V1,V2)\
 +  integrand1(18*dx,V1,V2)+ integrand1(-18*dx,V1,V2)\
 +  integrand1(19*dx,V1,V2)+ integrand1(-19*dx,V1,V2)\
 +  integrand1(20*dx,V1,V2)+ integrand1(-20*dx,V1,V2)\
)

P2trinom(V1,V2)=P1trinom(V2,V1)
P3trinom(V1,V2)=1. - P1trinom(V2,V1) - P2trinom(V2,V1)

print "P1trinom(0,0)=",P1trinom(0,0)
print "P1trinom(0,-10)=",P1trinom(0,-10)
print "P1trinom(0,10)=",P1trinom(0,10)
print "P1trinom(1,1)=",P1trinom(1,1)
print "P1trinom(10,10)=",P1trinom(10,10)


#########################################################
# Plotten
#########################################################

print ""

set nokey
set noparam

set isosample 60,60

set palette defined ( 0 "#aa88ff", 16 "#6666ff", \
       33 "#66ff66", 50 "yellow", 67 "orange", 84 "#ff6666", 100 "#ff0055")
unset surface

set pm3d       # bugfrei/am sichersten:  "set pm3d", DANN "set pm3d map" 
set pm3d; set pm3d map 
set contour surface

set cntrparam bspline 

unset clabel  # dann lauter gleiche Kontourlinien; 
                     # Farbe/Typ mit "w l ls" beim splot-Kommando

set grid   # grid/nogrid: whether 2D grid on xy plane 


#################################
# eigentliches Plotten
#################################

#################################
set out "pProbitTrinom_V1V2.eps"
print "plotting pProbitTrinom_V1V2.eps"
#################################

V1min=-3
V1max=4
V2min=-4
V2max=4

set size 1.0,1.1
yrel=0.97

set label 1 "i.i.d. Probit" at screen 0.4,yrel front
set label 2 "P_1" at screen 0.86,yrel front
set xlabel "V_1-V_3 (Units {/Symbol s_e})"
set xrange [V1min:V1max]

set ylabel "V_2-V_3 (Units {/Symbol s_e})" offset -1,0
set yrange [V2min:V2max]
unset cntrparam
set cntrparam levels discrete\
 0.02,0.02,0.06,0.08,0.099,\
 0.10,0.101,0.12,0.14,0.16,0.18,0.198,\
 0.20,0.202,0.22,0.24,0.26,0.28,0.298,\
 0.30,0.302,0.32,0.34,0.36,0.38,0.398,\
 0.40,0.402,0.42,0.44,0.46,0.48,0.498,\
 0.50,0.502,0.52,0.54,0.56,0.58,0.598,\
 0.60,0.602,0.62,0.64,0.66,0.68,0.698,\
 0.70,0.702,0.72,0.74,0.76,0.78,0.798,\
 0.80,0.802,0.82,0.84,0.86,0.88,0.899,\
 0.90,0.901,0.92,0.94,0.96,0.98

splot P1trinom(x,y) w l ls 99

#unset label 1
unset label 2

#################################
set out "kalProbitBinom_beta1beta2.eps"
print "plotting kalProbitBinom_beta1beta2.eps"
#################################

#======== Input  ==============
ni=5.  # 5 Wege pro befrage Person
I=6. # befragte Personen
# (T1i-T2i)    y_{1i}		 y_{2i}=ni-y_{1i}
dT1=-15.;  y11=3.;	y21=ni-y11
dT2=-5.;  y12=2.;	y22=ni-y12
dT3=0. ;  y13=1.;	y23=ni-y13
dT4=5. ;  y14=1.;	y24=ni-y14
dT5=10.;  y15=0.;	y25=ni-y15
dT6=30.;  y16=0.;	y26=ni-y16
#======== End input ============

# Zahl der insgesamt befragten Entscheidungen

n=ni*I  

# Nutzenfunktion Vki=(beta1+beta2*(T1i-T2i)) * delta(k,1)

dV(dT,beta1,beta2)=beta1+beta2*dT  # dV=V1-V2

# Probit-Wahrscheinlichkeit, falls Nutzenfunktion in Einheiten von 
# sigeps

# Log-Likelihood dieser Daten ohne konstanten Binomialkoeffizient-Term

logL(beta1,beta2)=\
     y11*log(P1binom(dV(dT1,beta1,beta2)))\
  + y21*log(P2binom(dV(dT1,beta1,beta2)))\
  + y12*log(P1binom(dV(dT2,beta1,beta2)))\
  + y22*log(P2binom(dV(dT2,beta1,beta2)))\
  + y13*log(P1binom(dV(dT3,beta1,beta2)))\
  + y23*log(P2binom(dV(dT3,beta1,beta2)))\
  + y14*log(P1binom(dV(dT4,beta1,beta2)))\
  + y24*log(P2binom(dV(dT4,beta1,beta2)))\
  + y15*log(P1binom(dV(dT5,beta1,beta2)))\
  + y25*log(P2binom(dV(dT5,beta1,beta2)))\
  + y16*log(P1binom(dV(dT6,beta1,beta2)))\
  + y26*log(P2binom(dV(dT6,beta1,beta2)))


beta1min=-2
beta1max=0
beta2min=-0.25
beta2max=0
zmin=-22.
zmax=-12.2

set label 1 "i.i.d. Probit" at screen 0.48,yrel front
set label 2 "ln L" at screen 0.90,yrel front

set xlabel "{/Symbol b}_1"
set xrange [beta1min:beta1max]

set ylabel "{/Symbol b}_2"
set yrange [beta2min:beta2max]

unset cntrparam
set cntrparam levels incr zmin,0.1*(zmax-zmin),zmax
#set cntrparam levels 10


set zrange [zmin:zmax]
set cbrange [zmin:zmax]

set palette defined ( 0 "#ffffff", 16 "#6666ff", \
       33 "#66ff66", 50 "yellow", 67 "orange", 84 "#ff6666",\
      99 "#ff0055", 100 "#770022")


set param
splot[x=beta1min:beta1max][y=beta2min:beta2max]\
  x,y,logL(x,y) w l ls 99


#################################
set out "kalProbitTrinom_beta1beta2.eps"
print "plotting kalProbitTrinom_beta1beta2.eps\n"
#################################

#======== Input  ==============
ni=5.  # 5 Wege pro befrage Person
I=6. # befragte Personen
# (T1i-T3i)  (T2i-T3i)      y_{1i}       y_{2i}   y_{3i}=ni-y_{1i}
dT11=-20.;  dT21=-20;  y11=3.;   y21=2;	y31=ni-y11-y21
dT12=-5.;   dT22=-10;  y12=1.;   y22=1; 	y32=ni-y12-y22
dT13=0. ;  dT23=1e3;   y13=1.;   y23=0; 	y33=ni-y13-y23
dT14=5. ;  dT24=-0.0;   y14=1.;   y24=0;	y34=ni-y14-y24
dT15=10.;  dT25=-20;   y15=0.;   y25=4; 	y35=ni-y15-y25
dT16=30.;  dT26=1e3;   y16=0.;   y26=0; 	y36=ni-y16-y26
#======== End input ============

# Zahl der insgesamt befragten Entscheidungen

n=ni*I  

# Nutzenfunktion Vki-V3i=(beta1+beta2*(Tki-T3i)) * (1-delta(k,3))

dV1(dT,beta1,beta2)=beta1+beta2*dT  # dV1=V1-V3
dV2(dT,beta1,beta2)=beta1+beta2*dT  # dV2=V2-V3 # identisch, da T nur gen. Var

# Probit-Wahrscheinlichkeit, falls Nutzenfunktion in Einheiten von sqrt(2*sigeps)

#P1trinom(V1,V2), P2trinom(V1,V2), P3trinom(V1,V2) von oben


# Log-Likelihood dieser Daten ohne konstanten Binomialkoeffizient-Term

logLtri(beta1,beta2)=\
  + y11*log(P1trinom(dV1(dT11,beta1,beta2), dV2(dT21,beta1,beta2)))\
  + y21*log(P2trinom(dV1(dT11,beta1,beta2), dV2(dT21,beta1,beta2)))\
  + y31*log(P3trinom(dV1(dT11,beta1,beta2), dV2(dT21,beta1,beta2)))\
  + y12*log(P1trinom(dV1(dT12,beta1,beta2), dV2(dT22,beta1,beta2)))\
  + y22*log(P2trinom(dV1(dT12,beta1,beta2), dV2(dT22,beta1,beta2)))\
  + y32*log(P3trinom(dV1(dT12,beta1,beta2), dV2(dT22,beta1,beta2)))\
  + y13*log(P1trinom(dV1(dT13,beta1,beta2), dV2(dT23,beta1,beta2)))\
  + y23*log(P2trinom(dV1(dT13,beta1,beta2), dV2(dT23,beta1,beta2)))\
  + y33*log(P3trinom(dV1(dT13,beta1,beta2), dV2(dT23,beta1,beta2)))\
  + y14*log(P1trinom(dV1(dT14,beta1,beta2), dV2(dT24,beta1,beta2)))\
  + y24*log(P2trinom(dV1(dT14,beta1,beta2), dV2(dT24,beta1,beta2)))\
  + y34*log(P3trinom(dV1(dT14,beta1,beta2), dV2(dT24,beta1,beta2)))\
  + y15*log(P1trinom(dV1(dT15,beta1,beta2), dV2(dT25,beta1,beta2)))\
  + y25*log(P2trinom(dV1(dT15,beta1,beta2), dV2(dT25,beta1,beta2)))\
  + y35*log(P3trinom(dV1(dT15,beta1,beta2), dV2(dT25,beta1,beta2)))\
  + y16*log(P1trinom(dV1(dT16,beta1,beta2), dV2(dT26,beta1,beta2)))\
  + y26*log(P2trinom(dV1(dT16,beta1,beta2), dV2(dT26,beta1,beta2)))\
  + y36*log(P3trinom(dV1(dT16,beta1,beta2), dV2(dT26,beta1,beta2)))



beta1min=-3
beta1max=0
beta2min=-0.3
beta2max=0
zmin=-27
zmax=-17

set xrange [beta1min:beta1max]
set yrange [beta2min:beta2max]
set zrange [zmin:zmax]
set cbrange [zmin:zmax]

set isosamples 20,20
unset cntrparam
set cntrparam levels incr zmin,0.1*(zmax-zmin),zmax

splot[x=beta1min:beta1max][y=beta2min:beta2max]\
  x,y,logLtri(x,y) w l ls 99


quit


print "Vergleich Probitmodell-Daten mit geschaetzten Parametern:"
b1=-0.7
b2=-0.08
print "dV1=",dV(dT1,b1,b2)," P11=",P1binom(dV(dT1,b1,b2))," f11=",y11/ni
print "dV2=",dV(dT2,b1,b2)," P12=",P1binom(dV(dT2,b1,b2))," f12=",y12/ni
print "dV3=",dV(dT3,b1,b2)," P13=",P1binom(dV(dT3,b1,b2))," f13=",y13/ni
print "dV4=",dV(dT4,b1,b2)," P14=",P1binom(dV(dT4,b1,b2))," f14=",y14/ni
print "dV5=",dV(dT5,b1,b2)," P15=",P1binom(dV(dT5,b1,b2))," f15=",y15/ni
print "dV6=",dV(dT6,b1,b2)," P16=",P1binom(dV(dT6,b1,b2))," f16=",y16/ni


print "Vergleich Probitmodell-Daten mit geschaetzten Parametern:"
b1=-1.4
b2=-0.13
dV11=dV1(dT11,b1,b2); dV21=dV2(dT21,b1,b2)
dV12=dV1(dT12,b1,b2); dV22=dV2(dT22,b1,b2)
dV13=dV1(dT13,b1,b2); dV23=dV2(dT23,b1,b2)
dV14=dV1(dT14,b1,b2); dV24=dV2(dT24,b1,b2)
dV15=dV1(dT15,b1,b2); dV25=dV2(dT25,b1,b2)
dV16=dV1(dT16,b1,b2); dV26=dV2(dT26,b1,b2)

P11=P1trinom(dV11,dV21); P21=P2trinom(dV11,dV21); P31=1-P11-P21
P12=P1trinom(dV12,dV22); P22=P2trinom(dV12,dV22); P32=1-P12-P22
P13=P1trinom(dV13,dV23); P23=P2trinom(dV13,dV23); P33=1-P13-P23
P14=P1trinom(dV14,dV24); P24=P2trinom(dV14,dV24); P34=1-P14-P24
P15=P1trinom(dV15,dV25); P25=P2trinom(dV15,dV25); P35=1-P15-P25
P16=P1trinom(dV16,dV26); P26=P2trinom(dV16,dV26); P36=1-P16-P26

print "dV11=",dV11," dV21=",dV21
print "  P11=",P11," P21=",P21," f11=",y11/ni," f21=",y21/ni
print "dV12=",dV12," dV22=",dV22
print "  P12=",P12," P22=",P22," f12=",y12/ni," f22=",y22/ni
print "dV13=",dV13," dV23=",dV23
print "  P13=",P13," P23=",P23," f13=",y13/ni," f23=",y23/ni
print "dV14=",dV14," dV24=",dV24
print "  P14=",P14," P24=",P24," f14=",y14/ni," f24=",y24/ni
print "dV15=",dV15," dV25=",dV25
print "  P15=",P15," P25=",P25," f15=",y15/ni," f25=",y25/ni
print "dV16=",dV16," dV26=",dV26
print "  P16=",P16," P26=",P26," f16=",y16/ni," f26=",y26/ni

