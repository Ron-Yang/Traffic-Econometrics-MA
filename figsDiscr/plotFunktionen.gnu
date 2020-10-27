


# siehe ~/info/gnuplot, ~/info/gnuTemplate.gnu

##########################################################
#geordnet nach hue (ps)
# set style line: linetype lt, point type pt

set encoding iso_8859_1 # dann -bäöüßÄÖÜ durch woertl. Eingabe korrekt-A

set style line 1 lt 7 lw 6 pt 1 ps 0.9 #schwarz, plus sign
set style line 2 lt 1 lw 8 pt 4 ps 0.5 #rot, open box
set style line 3 lt 8 lw 4 pt 5 ps 0.7 #blassrot, closed square
set style line 4 lt 6 lw 4 pt 6 ps 0.7 #gelb, open circle
set style line 5 lt 2 lw 4 pt 8  ps 0.7 #gruen, open triangle
set style line 6 lt 5 lw 4 pt 9  ps 0.7 #blasstuerkisblau, closed triangle
set style line 7 lt 3 lw 4 pt 10 ps 0.7 #blau, upside-down open triangle
set style line 8 lt 4 lw 4 pt 11 ps 1.4 #lila, upside-down closed triangle

set style line 11 lt 7 lw 1 pt 4 ps 1.2 #rot, open box
set style line 12 lt 1 lw 1 pt 4 ps 1.2 #rot, open box
set style line 13 lt 8 lw 1 pt 5 ps 1.5 #blassrot, closed square
set style line 14 lt 6 lw 1 pt 6 ps 1.5 #gelb, open circle
set style line 15 lt 2 lw 1 pt 8  ps 1.5 #gruen, open triangle
set style line 16 lt 5 lw 1 pt 9  ps 1.8 #blasstuerkisblau, closed triangle
set style line 17 lt 3 lw 1 pt 10 ps 1.5 #blau, upside-down open triangle
set style line 18 lt 4 lw 1 pt 11 ps 2.8 #lila, upside-down closed triangle



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
chi2Cum(x,nu)=igamma(0.5*nu,0.5*x)

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x

eulerKonst=0.577216

#gamma(n+1)=n! for integer n>=0

#erlang(1,x)=standard Poisson distribution

b_erlang(n)=gamma(n+1) /gamma(n)
a_erlang(n)=(b_erlang(n))**(n)/gamma(n)
erlang(n,x)=a_erlang(n)*x**(n-1)*exp(-b_erlang(n)*x)


set term post eps enhanced color solid "Helvetica" 24
set size 1,1
set nokey

#########################################################
set out "fGauss.eps"
print "plotting fGauss.eps"
#########################################################

set size 1,1
set key at screen 0.83,0.92
xmin=-3.1
xmax=3.1
set xrange [xmin:xmax]
set xlabel "Zufallsnutzen {/Symbol e}"
set ylabel "f_{Gauss}" offset 1,0
set yrange [0:0.49]
#set auto y
plot [x=xmin:xmax]\
 x, gauss(0,1,x) t "Dichtefunktion f_1({/Symbol e})=f_2({/Symbol e})" w l ls 7,\
 x, gauss(0,2**0.5,x) t "Faltung (f_1*f_2)({/Symbol e})" w l ls 2


#########################################################
set out "binProbit_P1.eps"
print "plotting binProbit_P1.eps"
#########################################################
set key at screen 0.92,0.35
set yrange [0:1.05]
set xlabel "Mittlere Nutzendifferenz V_1-V_2 (Einheiten {/Symbol s_e})"
set ylabel "Auswahlwahrscheinlichkeit P_1"
plot [x=xmin:xmax]\
 x, phi(x/2**0.5) t "P_1" w l ls 2,\
 x, phi(x) t "Verteilungsfunktion F_1" w l ls 7

#########################################################
set out "multiProbit_P1.eps"
print "plotting multiProbit_P1.eps"
#########################################################

# alle deterministischen Nuzten ausser V1 sind gleich Null !

integrand(x,V1,K)=gauss(0,1,x)*phi(x+V1)**(K-1)
integral(x1,x2,V1,K)=0.05*(x2-x1)*(integrand(x1,V1,K)+integrand(x2,V1,K))\
+0.1*(x2-x1)*(integrand(0.9*x1+0.1*x2,V1,K)\
        + integrand(0.8*x1+0.2*x2,V1,K)\
        + integrand(0.7*x1+0.3*x2,V1,K)\
        + integrand(0.6*x1+0.4*x2,V1,K)\
        + integrand(0.5*x1+0.5*x2,V1,K)\
        + integrand(0.4*x1+0.6*x2,V1,K)\
        + integrand(0.3*x1+0.7*x2,V1,K)\
        + integrand(0.2*x1+0.8*x2,V1,K)\
        + integrand(0.1*x1+0.9*x2,V1,K)\
)
print "integral(-2,2,0,3)=",integral(-2,2,0,3)
print "integral(-4,4,9,3)=",integral(-4,4,9,3)
print "integral(-3,3,9,3)=",integral(-3,3,9,3)

set key at screen 0.42,0.9
set xlabel "V_1-V_i (Einheiten {/Symbol s})"
xmax=4.
set xrange [xmin:xmax]
plot [x=xmin:xmax]\
 x, integral(-4,4,x,2) t "I=2" w l ls 2,\
 x, integral(-4,4,x,3) t  "I=3" w l ls 3,\
 x, integral(-4,4,x,5) t  "I=5" w l ls 5,\
 x, integral(-4,4,x,10) t  "I=10" w l ls 7,\
 x, integral(-4,4,x,20) t  "I=20" w l ls 1


#########################################################
set out "multiLogit_P1.eps"
print "plotting multiLogit_P1.eps"
#########################################################

# alle deterministischen Nuzten ausser V1 sind gleich Null, Varianz = 1 !
# \sigma=\frac{\pi}{\sqrt{6}\mu}.
set key at screen 0.51,0.92

mu=pi/sqrt(6.)
P1Logit(V1,K)=exp(mu*V1)/(exp(mu*V1)+K-1)
xmax=4.5
set xrange [xmin:xmax]
plot [x=xmin:xmax]\
 x, P1Logit(x,2) t "Logit, I=2" w l ls 2,\
 x, integral(-4,4,x,2) t "Probit, I=2" w l ls 12,\
 x, P1Logit(x,4) t "Logit, I=4" w l ls 3,\
 x, integral(-4,4,x,4) t "Probit, I=4" w l ls 13,\
 x, P1Logit(x,10) t "Logit, I=10" w l ls 7,\
 x, integral(-4,4,x,10) t  "Probit, I=10" w l ls 17,\
 x, P1Logit(x,50) t "Logit, I=50" w l ls 1,\
 x, integral(-4,4,x,50) t  "Probit, I=50" w l ls 11


#########################################################
set out "gumbel_f.eps"
print "plotting gumbel_f.eps"
#########################################################

set noparam

lambda=1.
theta(x)=(x>0) ? 1 : 0
FE(x,k)=theta(x)*(1.-exp(-lambda*x))**k
fE(x,k)=theta(x)*k*(1.-exp(-lambda*x))**(k-1) * lambda*exp(-lambda*x)

FGumbel(x,eta,lambda)=exp(-exp(-lambda*(x-eta)))
fGumbel(x,eta,lambda)=lambda*exp(-lambda*(x-eta))*FGumbel(x,eta,lambda)

etamax(k)=log(k+0.)/lambda

set key at screen 0.94,0.9

set ylabel "f_{Gu}(x)"
set yrange [0:0.8]
set xlabel "x"
set xrange [-2.2:5]
plot\
 fGumbel(x,0,1) w l ls 1,\
 fGumbel(x,0,2) w l ls 2,\
 fGumbel(x,1,1) w l ls 3

#########################################################
set out "chi2_f.eps"
print "plotting chi2_f.eps"
#########################################################
xmin=0
xmax=10.
set xrange [xmin:xmax]
set xlabel "x"
set ylabel "f(x)"
set yrange [0:0.6]
plot\
 chi2(x,1) t "{/Symbol c}^2(1)" w l ls 1,\
 chi2(x,2) t "{/Symbol c}^2(2)"  w l ls 2,\
 chi2(x,3) t "{/Symbol c}^2(3)"  w l ls 3,\
 chi2(x,4) t "{/Symbol c}^2(4)"  w l ls 7

#########################################################
set out "chi2_F.eps"
print "plotting chi2_F.eps"
#########################################################

set ylabel "F(x)"
set yrange [0:1.001]
#set ytics 0.1
set key at screen 0.9,0.4
plot\
 chi2Cum(x,1) t "{/Symbol c}^2(1)" w l ls 1,\
 chi2Cum(x,2) t "{/Symbol c}^2(2)"  w l ls 2,\
 chi2Cum(x,3) t "{/Symbol c}^2(3)"  w l ls 3,\
 chi2Cum(x,4) t "{/Symbol c}^2(4)"  w l ls 7,\
 0.95 t "F=0.95" w l ls 11
 
#########################################################
set out "chi2_F_detail.eps"
print "plotting chi2_F_detail.eps"
#########################################################
set xrange [3:20]
set yrange [0.94:1.0001]
#set ytics 0.005
plot\
 chi2Cum(x,1) t "{/Symbol c}^2(1)" w l ls 1,\
 chi2Cum(x,2) t "{/Symbol c}^2(2)"  w l ls 2,\
 chi2Cum(x,3) t "{/Symbol c}^2(3)"  w l ls 3,\
 chi2Cum(x,4) t "{/Symbol c}^2(4)"  w l ls 7,\
 0.95 t "F=0.95" w l ls 11
 
#########################################################
set out "chi2_F_detail2.eps"
print "plotting chi2_F_detail2.eps"
#########################################################
set xrange [2:16]
set yrange [0.80:1.0001]
#set ytics 0.02
plot\
 chi2Cum(x,1) t "{/Symbol c}^2(1)" w l ls 1,\
 chi2Cum(x,2) t "{/Symbol c}^2(2)"  w l ls 2,\
 chi2Cum(x,3) t "{/Symbol c}^2(3)"  w l ls 3,\
 chi2Cum(x,4) t "{/Symbol c}^2(4)"  w l ls 7,\
 0.95 t "F=0.95" w l ls 11
 
#########################################################
set out "gumbelGrenz_F.eps"
print "plotting gumbelGrenz_F.eps"
#########################################################

set samples 100


set key at screen 0.96,0.40
set noparam
xmin=-1
xmax=9.
set xrange [xmin:xmax]
set xlabel "x"
set ylabel "F(x)"
set yrange [0:1.05]

plot\
 FE(x,1) t "X_1" w p ls 1,\
 FGumbel(x,etamax(1),lambda) t "" w l ls 11,\
 FE(x,2) t "max(X_1,X_2)" w p ls 2,\
 FGumbel(x,etamax(2),lambda) t "" w l ls 12,\
 FE(x,20) t "max(X_1, ..., X_{20})" w p ls 7,\
 FGumbel(x,etamax(20),lambda) t "" w l ls 17,\
 FE(x,100) t "max(X_1, ..., X_{100})" w p ls 8,\
 FGumbel(x,etamax(100),lambda) t "Gu(ln I, 1)" w l ls 18


#########################################################
set out "gumbelGrenz_f.eps"
print "plotting gumbelGrenz_f.eps"
#########################################################

set key at screen 0.94,0.93
set ylabel "f(x)"
set yrange [0:0.7]
set samples 300

plot\
 fE(x,1) t "X_1" w p ls 1,\
 fGumbel(x,etamax(1),lambda) t "Gu(ln 1, 1" w l ls 11,\
 fE(x,2) t "max(X_1,X_2)" w p ls 2,\
 fGumbel(x,etamax(2),lambda) t "Gu(ln 2, 1" w l ls 12,\
 fE(x,20) t "max(X_1, ..., X_{20})" w p ls 7,\
 fGumbel(x,etamax(20),lambda) t "Gu(ln 20, 1" w l ls 17,\
 fE(x,100) t "max(X_1, ..., X_{100})" w p ls 8,\
 fGumbel(x,etamax(100),lambda) t "" w l ls 18
set samples 100

#########################################################
set out "gumbelGrenz2_F.eps"
print "plotting gumbelGrenz2_F.eps"
#########################################################

set samples 80

a=1.5
#const2=(1.-exp(-a*lambda))/a
fconst=lambda*exp(-lambda*a)
Fconst=1.-exp(-lambda*a)
dxconst=Fconst/fconst
xmin=a-dxconst
print "xmin=",xmin
print ""

FE2(x,k)=theta(x-xmin) * (  (x<a) ? fconst*(x-xmin) : 1.-exp(-lambda*x))**k
fE2(x,k)= theta(x-xmin)\
 * ((x<a) ? k*fconst*(fconst*(x-xmin))**(k-1) :  k*(1.-exp(-lambda*x))**(k-1) * lambda*exp(-lambda*x))
set key at screen 0.96,0.40
set xrange [xmin-0.5:12]
set ylabel "F(x)"
set yrange [0:1.05]

plot\
 FE2(x,1) t "X_1" w p ls 1,\
 FGumbel(x,etamax(1),lambda) t "" w l ls 11,\
 FE2(x,2) t "max(X_1,X_2)" w p ls 2,\
 FGumbel(x,etamax(2),lambda) t "" w l ls 12,\
 FE2(x,20) t "max(X_1, ..., X_{20})" w p ls 7,\
 FGumbel(x,etamax(20),lambda) t "" w l ls 17,\
 FE2(x,100) t "max(X_1, .., X_{100})" w p ls 8,\
 FGumbel(x,etamax(100),lambda) t "Gu(ln I, 1)" w l ls 18

#########################################################
set out "gumbelGrenz2_f.eps"
print "plotting gumbelGrenz2_f.eps"
#########################################################

set key at screen 0.96,0.92
set ylabel "f(x)"
set yrange [0:0.4]
plot\
 fE2(x,1) t "X_1" w p ls 1,\
 fGumbel(x,etamax(1),lambda) t "" w l ls 11,\
 fE2(x,2) t "max(X_1,X_2)" w p ls 2,\
 fGumbel(x,etamax(2),lambda) t "" w l ls 12,\
 fE2(x,20) t "max(X_1, ..., X_{20})" w p ls 7,\
 fGumbel(x,etamax(20),lambda) t "" w l ls 17,\
 fE2(x,100) t "max(X_1,..,X_{100})" w p ls 8,\
 fGumbel(x,etamax(100),lambda) t "Gu(ln I, 1)" w l ls 18


set samples 100

quit

####!! ab hier obsolet !!!#################

#########################################################
set out "fGumbel.eps"
#########################################################

#set size 0.8,1
set size 1,0.6
set auto y

set param
set key at screen 0.27,0.5
xmin=-10.
xmax=10.
set xlabel "z_1 bzw. z_2"
set ylabel "f_G"
#FGumbel(x,beta)=exp(-exp(-(beta*x+eulerKonst)))
#fGumbel(x,beta)=beta*exp(-(beta*x+eulerKonst))*FGumbel(x,beta)


set xrange [xmin:xmax]
set auto x
plot [x=xmin:xmax]\
 x, fGumbel(x,1) t "{/Symbol b}=1" w lp ls 2,\
 x, fGumbel(x,0.5) t "{/Symbol b}=1/2" w l ls 5,\
 x, fGumbel(x,0.25) t "{/Symbol b}=1/4" w l ls 7

#########################################################
set out "FGumbel.eps"
#########################################################

set ylabel "F_G"
plot [x=xmin:xmax]\
 x, FGumbel(x,1) t "{/Symbol b}=1" w lp ls 2,\
 x, FGumbel(x,0.5) t "{/Symbol b}=1/2" w l ls 5,\
 x, FGumbel(x,0.25) t "{/Symbol b}=1/4" w l ls 7

#########################################################
set out "logit_P1.eps"
#########################################################

#set size 1,1
set xlabel "Mittlere Nutzendifferenz V_1-V_2"
set ylabel "Auswahlwahrscheinlichkeit A_1"

P1_logit(u12,beta)=exp(beta*u12)/(1+exp(beta*u12))
plot [x=xmin:xmax]\
 x, P1_logit(x,1) t "{/Symbol b}=1" w lp ls 2,\
 x, P1_logit(x,0.5) t "{/Symbol b}=1/2" w l ls 5,\
 x, P1_logit(x,0.25) t "{/Symbol b}=1/4" w l ls 7

#########################################################
set out "AnteileUni.eps"
#########################################################

set size 1,1
set samples 500
set xlabel "Weglänge x (km)"
set ylabel "Modal-Split" offset 1,0

# Parameter des Multinom-Logit Modells:
# u*=det. Anteil der Utilities (neg. Zeiten in min.)
# beta=Parameter d. stoch. Anteils (1/beta approx Stdev in min)

betaConst(x)=0.2
beta(x)=betaConst(x) 

vFuss_kmh=5.
vRad_kmh=15.
vOEV_kmh=35.
vMIV_kmh=30.

vOEV2_kmh=30. # index 2 realistic for non-studis!
vMIV2_kmh=33.

# Ruestzeiten, Zugangszeiten, Wartezeiten etc

cFuss=0. 
cRad=5.
cOEV=10.
cMIV=16.
cOEV2=10.
cMIV2=10.

uFuss(x)=-(cFuss+60*x/vFuss_kmh)
uRad(x) =-(cRad +60*x/vRad_kmh)
uOEV(x) =-(cOEV +60*x/vOEV_kmh)
uMIV(x) =-(cMIV +60*x/vMIV_kmh)
uOEV2(x) =-(cOEV2 +60*x/vOEV2_kmh)
uMIV2(x) =-(cMIV2 +60*x/vMIV2_kmh)

sum(x) =exp(beta(x)*uFuss(x)) + exp(beta(x)*uRad(x))\
      + exp(beta(x)*uOEV(x)) + exp(beta(x)*uMIV(x))
sum2(x) =exp(beta(x)*uFuss(x)) + exp(beta(x)*uRad(x))\
      + exp(beta(x)*uOEV2(x)) + exp(beta(x)*uMIV2(x))

A_Fuss(x)=exp(beta(x)*uFuss(x))/sum(x)
A_Rad(x)=exp(beta(x)*uRad(x))/sum(x)
A_OEV(x)=exp(beta(x)*uOEV(x))/sum(x)
A_MIV(x)=exp(beta(x)*uMIV(x))/sum(x)

A2_Fuss(x)=exp(beta(x)*uFuss(x))/sum2(x)
A2_Rad(x)=exp(beta(x)*uRad(x))/sum2(x)
A2_OEV(x)=exp(beta(x)*uOEV2(x))/sum2(x)
A2_MIV(x)=exp(beta(x)*uMIV2(x))/sum2(x)

xmin=0
xmax=20
set xrange [xmin:xmax]
set key at screen 0.82,0.78

plot [x=xmin:xmax]\
 "AnteileUni.data" u 1:($2/($2+$3+$4+$5)) t "Data, Fuß" w p ls 17,\
 "AnteileUni.data" u 1:(($2+$3)/($2+$3+$4+$5)) t "Data, Fuß+Rad"  w p ls 15,\
 "AnteileUni.data" u 1:(($2+$3+$4)/($2+$3+$4+$5)) t "Data, Fuß+Rad+ÖV"  w p ls 13,\
 0,0 t "." w d 0,\
 0,0 t "Multinomial-Logit-Ansatz:" w d 0,\
 0,0 t "." w d 0,\
 x, A_Fuss(x) t "Fuß" w l ls 7,\
 x, A_Fuss(x)+A_Rad(x)  t "Fuß+Rad"  w l ls 5,\
 x, A_Fuss(x)+A_Rad(x)+A_OEV(x)  t "Fuß+Rad+ÖV"  w l ls 3
# x, A_MIV(x)  t "MIV"  w l ls 2


#########################################################
set out "AnteileLogitKum.eps"
#########################################################

xmaxKum=0.5*(xmin+xmax)
set nokey
set xrange [xmin:xmaxKum]
set label "Fuß" at xmin+0.06*(xmaxKum-xmin), 0.05
set label "Rad" at xmin+0.25*(xmaxKum-xmin), 0.15
set label "ÖV"  at xmin+0.45*(xmaxKum-xmin), 0.45
set label "MIV" at xmin+0.45*(xmaxKum-xmin), 0.95
plot [x=xmin:xmaxKum]\
 x, A_Fuss(x) t "Fuß" w l ls 7,\
 x, A_Fuss(x)+A_Rad(x)  t "F+R"  w l ls 5,\
 x, A_Fuss(x)+A_Rad(x)+A_OEV(x)  t "f+R+ÖV"  w l ls 3
set key

set nolabel

#########################################################
set out "AnteileUniLogitVar.eps"
#########################################################

#betaVar(x)=1./(4*x**0.1)
betaVar(x)=1./(4*x**0.2)
beta(x)=betaVar(x) ##!!! nun beta variabel !!!

## neue Def. sum(x), A_i(x) nicht noetig !!!


plot [x=xmin:xmax]\
 "AnteileUni.data" u 1:($2/($2+$3+$4+$5)) t "Data, Fuß" w p ls 17,\
 "AnteileUni.data" u 1:($3/($2+$3+$4+$5)) t "Data, Rad"  w p ls 15,\
 "AnteileUni.data" u 1:($4/($2+$3+$4+$5)) t "Data, ÖV"  w p ls 13,\
 "AnteileUni.data" u 1:($5/($2+$3+$4+$5)) t "Data, MIV"  w p ls 12,\
 0,0 t "." w d 0,\
 0,0 t "Multinomial-Logit-Ansatz:" w d 0,\
 0,0 t "." w d 0,\
 x, A_Fuss(x) t "Fuß" w l ls 7,\
 x, A_Rad(x)  t "Rad"  w l ls 5,\
 x, A_OEV(x)  t "ÖV"  w l ls 3,\
 x, A_MIV(x)  t "MIV"  w l ls 2
#########################################################
set out "AnteileUniKirchhoff.eps"
#########################################################

set key at screen 0.8,0.65

#Geaenderte kalibrierte Parameter notwendig!

vFuss_kmh=2.
vRad_kmh=5.
vOEV_kmh=50.
vMIV_kmh=20.

# Ruestzeiten, Zugangszeiten, Wartezeiten etc

cFuss=0. 
cRad=5.
cOEV=10.
cMIV=50.

uFuss(x)=-(cFuss+60*x/vFuss_kmh)
uRad(x) =-(cRad +60*x/vRad_kmh)
uOEV(x) =-(cOEV +60*x/vOEV_kmh)
uMIV(x) =-(cMIV +60*x/vMIV_kmh)

exp=1.
sum(x)=-(1./uFuss(x)**exp+1./uRad(x)**exp + 1./uOEV(x)**exp + 1./uMIV(x)**exp)
A_Fuss(x)=(-1./uFuss(x)**exp)/sum(x)
A_Rad(x) =(-1./uRad(x)**exp) /sum(x)
A_OEV(x) =(-1./uOEV(x)**exp) /sum(x)
A_MIV(x) =(-1./uMIV(x)**exp) /sum(x)


plot [x=xmin:xmax]\
 "AnteileUni.data" u 1:($2/($2+$3+$4+$5)) t "" w p ls 17,\
 "AnteileUni.data" u 1:(($2+$3)/($2+$3+$4+$5)) t ""  w p ls 15,\
 "AnteileUni.data" u 1:(($2+$3+$4)/($2+$3+$4+$5)) t ""  w p ls 13,\
 0,0 t "." w d 0,\
 0,0 t "Kirchhoffsche Regel linear:" w d 0,\
 0,0 t "." w d 0,\
 x, A_Fuss(x) t "Fuß" w l ls 7,\
 x, A_Fuss(x) +A_Rad(x)  t "Rad+Fuß"  w l ls 5,\
 x, A_Fuss(x) +A_Rad(x) + A_OEV(x)  t "ÖV+Rad+Fuß"  w l ls 3
# x, A_MIV(x)  t "MIV"  w l ls 2


#########################################################
set out "AnteileUniKirchhoff2.eps"
#########################################################

#Geaenderte kalibrierte Parameter notwendig!

vFuss_kmh=4.
vRad_kmh=15.
vOEV_kmh=40.
vMIV_kmh=20.

# Ruestzeiten, Zugangszeiten, Wartezeiten etc

cFuss=0. 
cRad=7.
cOEV=10.
cMIV=20.

uFuss(x)=-(cFuss+60*x/vFuss_kmh)
uRad(x) =-(cRad +60*x/vRad_kmh)
uOEV(x) =-(cOEV +60*x/vOEV_kmh)
uMIV(x) =-(cMIV +60*x/vMIV_kmh)

exp=3.
sum(x)=-(1./uFuss(x)**exp+1./uRad(x)**exp + 1./uOEV(x)**exp + 1./uMIV(x)**exp)
A_Fuss(x)=(-1./uFuss(x)**exp)/sum(x)
A_Rad(x) =(-1./uRad(x)**exp) /sum(x)
A_OEV(x) =(-1./uOEV(x)**exp) /sum(x)
A_MIV(x) =(-1./uMIV(x)**exp) /sum(x)

plot [x=xmin:xmax]\
 "AnteileUni.data" u 1:($2/($2+$3+$4+$5)) t "" w p ls 17,\
  "AnteileUni.data" u 1:(($2+$3)/($2+$3+$4+$5)) t ""  w p ls 15,\
   "AnteileUni.data" u 1:(($2+$3+$4)/($2+$3+$4+$5)) t ""  w p ls 13,\
    0,0 t "." w d 0,\
     0,0 t "Kirchhoffsche Regel nichtlinear:" w d 0,\
      0,0 t "." w d 0,\
       x, A_Fuss(x) t "Fuß" w l ls 7,\
        x, A_Fuss(x) +A_Rad(x)  t "Rad+Fuß"  w l ls 5,\
	 x, A_Fuss(x) +A_Rad(x) + A_OEV(x)  t "ÖV+Rad+Fuß"  w l ls 3
	 # x, A_MIV(x)  t "MIV"  w l ls 2
	 #

#########################################################
#########################################################
# mai2006
#########################################################


set samples 200

vFuss2_kmh=5.
vRad2_kmh=15.
vOEV2_kmh=30. # index 2 realistic for non-studis!
vMIV2_kmh=33.

# Ruestzeiten, Zugangszeiten, Wartezeiten etc

cFuss2=0. 
cRad2=5.
cOEV2=10.
cMIV2=10.

uFuss2(x)=-(cFuss2+60*x/vFuss2_kmh)
uRad2(x) =-(cRad2 +60*x/vRad2_kmh)
uOEV2(x) =-(cOEV2 +60*x/vOEV2_kmh)
uMIV2(x) =-(cMIV2 +60*x/vMIV2_kmh)

#########################################################
set out "AnteileLogitVarKum2.eps"
#########################################################

sum2(x) =exp(beta(x)*uFuss2(x)) + exp(beta(x)*uRad2(x))\
      + exp(beta(x)*uOEV2(x)) + exp(beta(x)*uMIV2(x))

A2_Fuss(x)=exp(beta(x)*uFuss2(x))/sum2(x)
A2_Rad(x)=exp(beta(x)*uRad2(x))/sum2(x)
A2_OEV(x)=exp(beta(x)*uOEV2(x))/sum2(x)
A2_MIV(x)=exp(beta(x)*uMIV2(x))/sum2(x)

beta(x)=1./(4*x**0.3)
xmax=10
xmaxKum=0.5*(xmin+xmax)
set nokey
set xrange [xmin:xmaxKum]
set label "Fuß" at xmin+0.06*(xmaxKum-xmin), 0.05
set label "Rad" at xmin+0.25*(xmaxKum-xmin), 0.25
set label "ÖV"  at xmin+0.45*(xmaxKum-xmin), 0.50
set label "MIV" at xmin+0.45*(xmaxKum-xmin), 0.90
plot [x=xmin:xmaxKum]\
 x, A2_Fuss(x) t "Fuß" w l ls 7,\
 x, A2_Fuss(x)+A2_Rad(x)  t "F+R"  w l ls 5,\
 x, A2_Fuss(x)+A2_Rad(x)+A2_OEV(x)  t "f+R+ÖV"  w l ls 3
set nolabel

#########################################################
set out "AnteileLogitVarKum2_large.eps"
#########################################################

xmax=50
xmaxKum=0.5*(xmin+xmax)
xmaxLabel=7.
set nokey
set xrange [xmin:xmaxKum]
set label "Fuß" at xmin+0.06*(xmaxLabel-xmin), 0.05
set label "Rad" at xmin+0.25*(xmaxLabel-xmin), 0.25
set label "ÖV"  at xmin+0.45*(xmaxLabel-xmin), 0.50
set label "MIV" at xmin+0.45*(xmaxKum-xmin), 0.90
plot [x=xmin:xmaxKum]\
 x, A2_Fuss(x) t "Fuß" w l ls 7,\
 x, A2_Fuss(x)+A2_Rad(x)  t "F+R"  w l ls 5,\
 x, A2_Fuss(x)+A2_Rad(x)+A2_OEV(x)  t "f+R+ÖV"  w l ls 3
set nolabel

#########################################################
set out "AnteileLogitKum2.eps"
#########################################################

beta(x)=0.2
xmax=10
xmaxKum=0.5*(xmin+xmax)
set nokey
set nolabel
set xrange [xmin:xmaxKum]
set label "Fuß" at xmin+0.06*(xmaxKum-xmin), 0.05
set label "Rad" at xmin+0.25*(xmaxKum-xmin), 0.25
set label "ÖV"  at xmin+0.45*(xmaxKum-xmin), 0.50
set label "MIV" at xmin+0.45*(xmaxKum-xmin), 0.90
plot [x=xmin:xmaxKum]\
 x, A2_Fuss(x) t "Fuß" w l ls 7,\
 x, A2_Fuss(x)+A2_Rad(x)  t "F+R"  w l ls 5,\
 x, A2_Fuss(x)+A2_Rad(x)+A2_OEV(x)  t "f+R+ÖV"  w l ls 3

set nolabel


#########################################################
set out "AnteileLogitKum2_large.eps"
#########################################################

xmax=50
xmaxKum=0.5*(xmin+xmax)
set nokey
set nolabel
set xrange [xmin:xmaxKum]
set label "Fuß" at xmin+0.06*(xmaxLabel-xmin), 0.05
set label "Rad" at xmin+0.25*(xmaxLabel-xmin), 0.25
set label "ÖV"  at xmin+0.45*(xmaxLabel-xmin), 0.50
set label "MIV" at xmin+0.45*(xmaxKum-xmin), 0.90
plot [x=xmin:xmaxKum]\
 x, A2_Fuss(x) t "Fuß" w l ls 7,\
 x, A2_Fuss(x)+A2_Rad(x)  t "F+R"  w l ls 5,\
 x, A2_Fuss(x)+A2_Rad(x)+A2_OEV(x)  t "f+R+ÖV"  w l ls 3

set nolabel

#########################################################
set out "AnteileLogitBetarelKum2.eps"
#########################################################

beta(x)=1./(2*x**1.)
xmax=10
xmaxKum=0.5*(xmin+xmax)
set nokey
set nolabel
set xrange [xmin:xmaxKum]
set label "Fuß" at xmin+0.06*(xmaxKum-xmin), 0.05
set label "Rad" at xmin+0.25*(xmaxKum-xmin), 0.25
set label "ÖV"  at xmin+0.45*(xmaxKum-xmin), 0.50
set label "MIV" at xmin+0.45*(xmaxKum-xmin), 0.90
plot [x=xmin:xmaxKum]\
 x, A2_Fuss(x) t "Fuß" w l ls 7,\
 x, A2_Fuss(x)+A2_Rad(x)  t "F+R"  w l ls 5,\
 x, A2_Fuss(x)+A2_Rad(x)+A2_OEV(x)  t "f+R+ÖV"  w l ls 3
set nolabel


#########################################################
set out "AnteileLogitBetarelKum2_large.eps"
#########################################################

xmax=50
xmaxKum=0.5*(xmin+xmax)
set nokey
set xrange [xmin:xmaxKum]
set label "Fuß" at xmin+0.06*(xmaxLabel-xmin), 0.05
set label "Rad" at xmin+0.25*(xmaxLabel-xmin), 0.25
set label "ÖV"  at xmin+0.45*(xmaxLabel-xmin), 0.50
set label "MIV" at xmin+0.45*(xmaxKum-xmin), 0.90
plot [x=xmin:xmaxKum]\
 x, A2_Fuss(x) t "Fuß" w l ls 7,\
 x, A2_Fuss(x)+A2_Rad(x)  t "F+R"  w l ls 5,\
 x, A2_Fuss(x)+A2_Rad(x)+A2_OEV(x)  t "f+R+ÖV"  w l ls 3

set key
set nolabel
xmax=20

#########################################################
set out "AnteileKirchhoffKum2.eps"
#########################################################

xmax=10
xmaxKum=0.5*(xmin+xmax)
set nokey
set xrange [xmin:xmaxKum]

exp=1.
sum2(x)=-(1./uFuss2(x)**exp+1./uRad2(x)**exp\
          + 1./uOEV2(x)**exp + 1./uMIV2(x)**exp)
A2_Fuss(x)=(-1./uFuss2(x)**exp)/sum2(x)
A2_Rad(x) =(-1./uRad2(x)**exp) /sum2(x)
A2_OEV(x) =(-1./uOEV2(x)**exp) /sum2(x)
A2_MIV(x) =(-1./uMIV2(x)**exp) /sum2(x)
set label "Fuß" at xmin+0.06*(xmaxKum-xmin), 0.05
set label "Rad" at xmin+0.25*(xmaxKum-xmin), 0.25
set label "ÖV"  at xmin+0.45*(xmaxKum-xmin), 0.50
set label "MIV" at xmin+0.45*(xmaxKum-xmin), 0.90
plot [x=xmin:xmaxKum]\
 x, A2_Fuss(x) t "Fuß" w l ls 7,\
 x, A2_Fuss(x)+A2_Rad(x)  t "F+R"  w l ls 5,\
 x, A2_Fuss(x)+A2_Rad(x)+A2_OEV(x)  t "f+R+ÖV"  w l ls 3
set nolabel

#########################################################
set out "AnteileKirchhoffKum2_large.eps"
#########################################################

xmax=50
xmaxKum=0.5*(xmin+xmax)
set xrange [xmin:xmaxKum]
set nokey
set label "Fuß" at xmin+0.06*(xmaxLabel-xmin), 0.05
set label "Rad" at xmin+0.25*(xmaxLabel-xmin), 0.25
set label "ÖV"  at xmin+0.45*(xmaxLabel-xmin), 0.50
set label "MIV" at xmin+0.45*(xmaxKum-xmin), 0.90
plot [x=xmin:xmaxKum]\
 x, A2_Fuss(x) t "Fuß" w l ls 7,\
 x, A2_Fuss(x)+A2_Rad(x)  t "F+R"  w l ls 5,\
 x, A2_Fuss(x)+A2_Rad(x)+A2_OEV(x)  t "f+R+ÖV"  w l ls 3

set key
set nolabel
xmax=20

#########################################################
set out "AnteileEVAKum2.eps"
#########################################################

xmax=10
xmaxKum=0.5*(xmin+xmax)
set nokey
set xrange [xmin:xmaxKum]
set yrange [0:1]

vFuss2_kmh=5.
vRad2_kmh=15.
vOEV2_kmh=30. # index 2 realistic for non-studis!
vMIV2_kmh=33.

cFuss2=0.
cRad2=5. # EVA enthaelt Ruestzeit ??? doch nicht!
cOEV2=10.
cMIV2=10.

uFuss2(x)=-(cFuss2+60*x/vFuss2_kmh)
uRad2(x) =-(cRad2 +60*x/vRad2_kmh)
uOEV2(x) =-(cOEV2 +60*x/vOEV2_kmh)
uMIV2(x) =-(cMIV2 +60*x/vMIV2_kmh)

E=3.
F=5.
w0=20.
phi(w)=E/(1.+exp(F*(1.-w/w0)))
eva(w)=1./((1.+w)**phi(w))

sum2(x)=eva(-uFuss2(x))+eva(-uRad2(x))+eva(-uOEV2(x))+eva(-uMIV2(x))
A2_Fuss(x)=eva(-uFuss2(x))/sum2(x)
A2_Rad(x)  =eva(-uRad2(x))/sum2(x)
A2_OEV(x)  =eva(-uOEV2(x))/sum2(x)
A2_MIV(x)  =eva(-uMIV2(x))/sum2(x)

print "x=20: uRad2(20)=",uRad2(20)," phi(-uRad2(20))=",phi(-uRad2(20))
print "x=20: eva(-uRad2(20))=",eva(-uRad2(20))," A2_Rad(20)=",A2_Rad(20)
print "x=40: uRad2(40)=",uRad2(40)," phi(-uRad2(40))=",phi(-uRad2(40))
print "x=40: eva(-uRad2(40))=",eva(-uRad2(40))," A2_Rad(40)=",A2_Rad(40)
print "x=20: uOEV2(20)=",uOEV2(20)," phi(-uOEV2(20))=",phi(-uOEV2(20))
print "x=20: eva(-uOEV2(20))=",eva(-uOEV2(20))," A2_OEV(20)=",A2_OEV(20)
print "x=40: uOEV2(40)=",uOEV2(40)," phi(-uOEV2(40))=",phi(-uOEV2(40))
print "x=40: eva(-uOEV2(40))=",eva(-uOEV2(40))," A2_OEV(40)=",A2_OEV(40)
set label "Fuß" at xmin+0.06*(xmaxKum-xmin), 0.05
set label "Rad" at xmin+0.25*(xmaxKum-xmin), 0.25
set label "ÖV"  at xmin+0.45*(xmaxKum-xmin), 0.50
set label "MIV" at xmin+0.45*(xmaxKum-xmin), 0.90
plot [x=xmin:xmaxKum]\
 x, A2_Fuss(x) t "Fuß" w l ls 7,\
 x, A2_Fuss(x)+A2_Rad(x)  t "F+R"  w l ls 5,\
 x, A2_Fuss(x)+A2_Rad(x)+A2_OEV(x)  t "f+R+ÖV"  w l ls 3
set nolabel

#########################################################
set out "AnteileEVAKum2_large.eps"
#########################################################

xmax=50
xmaxKum=0.5*(xmin+xmax)
set xrange [xmin:xmaxKum]
set nokey
set label "Fuß" at xmin+0.06*(xmaxLabel-xmin), 0.05
set label "Rad" at xmin+0.25*(xmaxLabel-xmin), 0.25
set label "ÖV"  at xmin+0.45*(xmaxLabel-xmin), 0.50
set label "MIV" at xmin+0.45*(xmaxKum-xmin), 0.90
plot [x=xmin:xmaxKum]\
 x, A2_Fuss(x) t "Fuß" w l ls 7,\
 x, A2_Fuss(x)+A2_Rad(x)  t "F+R"  w l ls 5,\
 x, A2_Fuss(x)+A2_Rad(x)+A2_OEV(x)  t "f+R+ÖV"  w l ls 3

set key
set nolabel
xmax=20

#########################################################
#########################################################

print "Kalibrierung Vorlesungsbeispiel Logit:"

vFuss=5.
vRad=15.
x1=0.5
x2=2.
x3=4.
x4=7.5
T1Fuss=60.*x1/vFuss
T2Fuss=60.*x2/vFuss
T3Fuss=60.*x3/vFuss
T4Fuss=60.*x4/vFuss
T1Rad=60.*x1/vRad
T2Rad=60.*x2/vRad
T3Rad=60.*x3/vRad
T4Rad=60.*x4/vRad
afar1=3.5
afar2=0.8
afar3=0.2
afar4=0.02

y1=log(afar1)
y2=log(afar2)
y3=log(afar3)
y4=log(afar4)

print "Klasse1: x1=",x1,"  Delta T-Delta T0=",(T1Fuss-T1Rad),"  log(afar1)=",y1
print "Klasse2: x2=",x2,"  Delta T-Delta T0=",(T2Fuss-T2Rad),"  log(afar2)=",y2
print "Klasse3: x3=",x3,"  Delta T-Delta T0=",(T3Fuss-T3Rad),"  log(afar3)=",y3
print "Klasse4: x4=",x4,"  Delta T-Delta T0=",(T4Fuss-T4Rad),"  log(afar4)=",y4

print "\nKalibrierung Vorlesungsbeispiel Kirchhoff:"
T0Rad=10.
z1=T1Fuss/(T0Rad+T1Rad)
z2=T2Fuss/(T0Rad+T2Rad)
z3=T3Fuss/(T0Rad+T3Rad)
z4=T4Fuss/(T0Rad+T4Rad)
y1=1./(afar1)
y2=1./(afar2)
y3=1./(afar3)
y4=1./(afar4)

print "Klasse1: x1=",x1,"  (T-T0)_Fuss/(T-T0)_Rad=",z1,"  1/(afar1)=",y1
print "Klasse2: x2=",x2,"  (T-T0)_Fuss/(T-T0)_Rad=",z2,"  1/(afar2)=",y2
print "Klasse3: x3=",x3,"  (T-T0)_Fuss/(T-T0)_Rad=",z3,"  1/(afar3)=",y3
print "Klasse4: x4=",x4,"  (T-T0)_Fuss/(T-T0)_Rad=",z4,"  1/(afar4)=",y4
