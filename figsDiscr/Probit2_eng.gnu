
# Trinomiales in Probit.gnu

#======== Input  ==============
ni=5.  # 5 Wege pro befrage Person
I=6. # befragte Personen
# (T1i-T2i)    y_{1i}		 y_{2i}=ni-y_{1i}
dT1=-15.;  y11=4.;	y21=ni-y11
dT2=-5.;  y12=1.;	y22=ni-y12
dT3=0. ;  y13=3.;	y23=ni-y13
dT4=5. ;  y14=2.;	y24=ni-y14
dT5=10.;  y15=1.;	y25=ni-y15
dT6=30.;  y16=0.;	y26=ni-y16
#======== End input ============

# Zahl der insgesamt befragten Entscheidungen

n=ni*I  


##########################################################
#geordnet nach hue (ps)
# set style line: linetype lt, point type pt

set encoding iso_8859_1 # dann -bäöüßÄÖÜ durch woertl. Eingabe korrekt-A

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



#########################################################
# Plotten
#########################################################

print ""

set term post eps enhanced color solid "Helvetica" 30
set nokey
set noparam

set isosample 60,60


 # bugfrei/am sichersten:  "set pm3d", DANN "set pm3d map" 
set pm3d; set pm3d map 
unset surface
set contour surface

set cntrparam bspline 

unset clabel  # dann lauter gleiche Kontourlinien; 
                     # Farbe/Typ mit "w l ls" beim splot-Kommando

set grid   # grid/nogrid: whether 2D grid on xy plane 
set size 1.5,1.5


#################################
# eigentliches Plotten
#################################

#################################
set out "kalProbitBinom2_beta1beta2.eps"
print "plotting kalProbitBinom2_beta1beta2.eps\n"
#################################


# Nutzenfunktion Vki=(beta1+beta2*(T1i-T2i)) * delta(k,1)

dV(dT,beta1,beta2)=beta1+beta2*dT  # dV=V1-V2

# Probit-Wahrscheinlichkeit, falls Nutzenfunktion in Einheiten von sqrt(2*sigeps)

P1binom(dV)=phi(dV)
P2binom(dV)=1.-phi(dV)

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

beta1min=-1.5
beta1max=0.2
beta2min=-0.15
beta2max=0.02
zmin=-22.02
zmax=-16.
unset cntrparam
set cntrparam levels incr zmin, 0.5, zmax
#set cntrparam levels incr zmin,0.1*(zmax-zmin),zmax
#set cntrparam levels 10


set xlabel "{/Symbol b}_1"
set xrange [beta1min:beta1max]

set ylabel "{/Symbol b}_2"
set yrange [beta2min:beta2max]



set zrange [zmin:zmax]
set cbrange [zmin:zmax]

set palette defined ( 0 "#ffffff", 16 "#6666ff", \
       33 "#66ff66", 50 "yellow", 67 "orange", 84 "#ff6666",\
       99 "#ff0055", 100 "#770022")

scalefactBi=sqrt(2.) #!! im Plot andersrum wie bei Logit, da hier Skalierung
# nicht in Nutzenfunktion, sondern beim Ausrechnen der Differenz eps_1=eps_2

set param
splot[x=beta1min:beta1max][y=beta2min:beta2max]\
  x*scalefactBi,y*scalefactBi,logL(x,y) w l lt 6 lw 4 

