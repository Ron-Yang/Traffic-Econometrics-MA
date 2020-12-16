


# siehe ~/info/gnuplot, ~/info/gnuTemplate.gnu



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
set xlabel "Random utility {/Symbol e}"
set ylabel "f_{Gauss}" offset 1,0
set yrange [0:0.49]
#set auto y
plot [x=xmin:xmax]\
 x, gauss(0,1,x) t "Density f_1({/Symbol e})=f_2({/Symbol e})" w l ls 7,\
 x, gauss(0,2**0.5,x) t "Convolution (f_1*f_2)({/Symbol e})" w l ls 2


#########################################################
set out "binProbit_P1.eps"
print "plotting binProbit_P1.eps"
#########################################################
set key at screen 0.92,0.35
set yrange [0:1.05]
set xlabel "Avg. utility difference V_1-V_2 (units {/Symbol s_e})"
set ylabel "Choice probability P_1"
plot [x=xmin:xmax]\
 x, phi(x/2**0.5) t "P_1" w l ls 2,\
 x, phi(x) t "Distribution function F_1" w l ls 7

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
set xlabel "V_1-V_i (units {/Symbol s})"
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
 x, integral(-4,4,x,10) t  "Probit, I=10" w l ls 16,\
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

set key at screen 0.96,0.40
set noparam
xmin=-1
xmax=9.
set xrange [xmin:xmax]
set xlabel "x"
set ylabel ""
set label 1 "F(x)" at screen 0.8,0.65
set yrange [0:1.05]

plot\
 FE(x,1) t "X_1" w l ls 1,\
 FGumbel(x,etamax(1),lambda) t "" w l ls 11,\
 FE(x,2) t "max(X_1,X_2)" w l ls 3,\
 FGumbel(x,etamax(2),lambda) t "" w l ls 13,\
 FE(x,10) t "max(X_1, ..., X_{10})" w l ls 5,\
 FGumbel(x,etamax(10),lambda) t "" w l ls 15,\
 FE(x,50) t "max(X_1, ..., X_{50})" w l ls 7,\
 FGumbel(x,etamax(50),lambda) t "Gumbel(ln I, 1)" w l ls 17

set out "gumbelGrenz_F1.eps"
print "plotting gumbelGrenz_F1.eps"
plot\
 FE(x,1) t "X_1" w l ls 1,\
 FGumbel(x,etamax(1),lambda) t "Gumbel(ln 1, 1)" w l ls 11

set out "gumbelGrenz_F2.eps"
print "plotting gumbelGrenz_F2.eps"
plot\
 FE(x,1) t "X_1" w l ls 1,\
 FGumbel(x,etamax(1),lambda) t "Gumbel(ln 1, 1)" w l ls 11,\
 FE(x,2) t "max(X_1,X_2)" w l ls 3,\
 FGumbel(x,etamax(2),lambda) t "Gumbel(ln 2, 1)" w l ls 13


set out "gumbelGrenz_F3.eps"
print "plotting gumbelGrenz_F3.eps"
plot\
 FE(x,1) t "X_1" w l ls 1,\
 FGumbel(x,etamax(1),lambda) t "Gumbel(ln 1, 1)" w l ls 11,\
 FE(x,2) t "max(X_1,X_2)" w l ls 3,\
 FGumbel(x,etamax(2),lambda) t "Gumbel(ln 2, 1)" w l ls 13,\
 FE(x,10) t "max(X_1, ..., X_{10})" w l ls 5,\
 FGumbel(x,etamax(10),lambda) t "Gumbel(ln 10, 1)" w l ls 15



#########################################################
set out "gumbelGrenz_f.eps"
print "plotting gumbelGrenz_f.eps"
#########################################################

set key at screen 0.94,0.93
set ylabel ""
set label 1 ""
set label 2 "f(x)" at screen 0.45,0.88
set yrange [0:0.7]


plot\
 fE(x,1) t "X_1" w l ls 1,\
 fGumbel(x,etamax(1),lambda) t "Gu(ln 1, 1)" w l ls 11,\
 fE(x,2) t "max(X_1,X_2)" w l ls 3,\
 fGumbel(x,etamax(2),lambda) t "Gu(ln 2, 1)" w l ls 13,\
 fE(x,10) t "max(X_1, ..., X_{10})" w l ls 5,\
 fGumbel(x,etamax(10),lambda) t "Gu(ln 10, 1)" w l ls 15,\
 fE(x,50) t "max(X_1, ..., X_{50})" w l ls 7,\
 fGumbel(x,etamax(50),lambda) t "Gu(ln 50, 1)" w l ls 17

set out "gumbelGrenz_f1.eps"
print "plotting gumbelGrenz_f1.eps"
plot\
 fE(x,1) t "X_1" w l ls 1,\
 fGumbel(x,etamax(1),lambda) t "Gu(ln 1, 1)" w l ls 11

set out "gumbelGrenz_f2.eps"
print "plotting gumbelGrenz_f2.eps"
plot\
 fE(x,1) t "X_1" w l ls 1,\
 fGumbel(x,etamax(1),lambda) t "Gu(ln 1, 1)" w l ls 11,\
 fE(x,2) t "max(X_1,X_2)" w l ls 3,\
 fGumbel(x,etamax(2),lambda) t "Gu(ln 2, 1)" w l ls 13

set out "gumbelGrenz_f3.eps"
print "plotting gumbelGrenz_f3.eps"
plot\
 fE(x,1) t "X_1" w l ls 1,\
 fGumbel(x,etamax(1),lambda) t "Gu(ln 1, 1)" w l ls 11,\
 fE(x,2) t "max(X_1,X_2)" w l ls 3,\
 fGumbel(x,etamax(2),lambda) t "Gu(ln 2, 1)" w l ls 13,\
 fE(x,10) t "max(X_1, ..., X_{10})" w l ls 5,\
 fGumbel(x,etamax(10),lambda) t "Gu(ln 10, 1)" w l ls 15


#########################################################
set out "gumbelGrenz2_F.eps"
print "plotting gumbelGrenz2_F.eps"
#########################################################

set key at screen 0.92,0.56
a=1.5
#const2=(1.-exp(-a*lambda))/a
fconst=lambda*exp(-lambda*a)
Fconst=1.-exp(-lambda*a)
dxconst=Fconst/fconst
xmin=a-dxconst

FE2(x,k)=theta(x-xmin) * (  (x<a) ? fconst*(x-xmin) : 1.-exp(-lambda*x))**k
fE2(x,k)= theta(x-xmin)\
 * ((x<a) ? k*fconst*(fconst*(x-xmin))**(k-1) :  k*(1.-exp(-lambda*x))**(k-1) * lambda*exp(-lambda*x))
set xrange [xmin-0.5:12]
set ylabel ""
set label 2 ""
set label 1 "F(x)"
set yrange [0:1.05]

plot\
 FE2(x,1) t "X_1" w l ls 1,\
 FGumbel(x,etamax(1),lambda) t "Gu(ln 1, 1)" w l ls 11,\
 FE2(x,2) t "max(X_1,X_2)" w l ls 3,\
 FGumbel(x,etamax(2),lambda) t "Gu(ln 2, 1)" w l ls 13,\
 FE2(x,10) t "max(X_1, ..., X_{10})" w l ls 5,\
 FGumbel(x,etamax(10),lambda) t "Gu(ln 10, 1)" w l ls 15,\
 FE2(x,50) t "max(X_1, .., X_{50})" w l ls 7,\
 FGumbel(x,etamax(50),lambda) t "Gu(ln 50, 1)" w l ls 17


set out "gumbelGrenz2_F1.eps"
print "plotting gumbelGrenz2_F1.eps"
plot\
 FE2(x,1) t "X_1" w l ls 1,\
 FGumbel(x,etamax(1),lambda) t "Gu(ln 1, 1)" w l ls 11

set out "gumbelGrenz2_F2.eps"
print "plotting gumbelGrenz2_F2.eps"
plot\
 FE2(x,1) t "X_1" w l ls 1,\
 FGumbel(x,etamax(1),lambda) t "Gu(ln 1, 1)" w l ls 11,\
 FE2(x,2) t "max(X_1,X_2)" w l ls 3,\
 FGumbel(x,etamax(2),lambda) t "Gu(ln 2, 1)" w l ls 13

set out "gumbelGrenz2_F3.eps"
print "plotting gumbelGrenz2_F3.eps"
plot\
 FE2(x,1) t "X_1" w l ls 1,\
 FGumbel(x,etamax(1),lambda) t "Gu(ln 1, 1)" w l ls 11,\
 FE2(x,2) t "max(X_1,X_2)" w l ls 3,\
 FGumbel(x,etamax(2),lambda) t "Gu(ln 2, 1)" w l ls 13,\
 FE2(x,10) t "max(X_1, ..., X_{10})" w l ls 5,\
 FGumbel(x,etamax(10),lambda) t "Gu(ln 10, 1)" w l ls 15



#########################################################
set out "gumbelGrenz2_f.eps"
print "plotting gumbelGrenz2_f.eps"
#########################################################

set key at screen 0.96,0.92
set ylabel ""
set label 1 ""
set label 2 "f(x)"
set yrange [0:0.4]
plot\
 fE2(x,1) t "X_1" w l ls 1,\
 fGumbel(x,etamax(1),lambda) t "" w l ls 11,\
 fE2(x,2) t "max(X_1,X_2)" w l ls 3,\
 fGumbel(x,etamax(2),lambda) t "" w l ls 13,\
 fE2(x,10) t "max(X_1, ..., X_{10})" w l ls 5,\
 fGumbel(x,etamax(10),lambda) t "" w l ls 15,\
 fE2(x,50) t "max(X_1,..,X_{50})" w l ls 7,\
 fGumbel(x,etamax(50),lambda) t "Gu(ln I, 1)" w l ls 17


set out "gumbelGrenz2_f1.eps"
print "plotting gumbelGrenz2_f1.eps"
plot\
 fE2(x,1) t "X_1" w l ls 1,\
 fGumbel(x,etamax(1),lambda) t "Gu(ln 1, 1)" w l ls 11

set out "gumbelGrenz2_f2.eps"
print "plotting gumbelGrenz2_f2.eps"
plot\
 fE2(x,1) t "X_1" w l ls 1,\
 fGumbel(x,etamax(1),lambda) t "Gu(ln 1, 1)" w l ls 11,\
 fE2(x,2) t "max(X_1,X_2)" w l ls 3,\
 fGumbel(x,etamax(2),lambda) t "Gu(ln 2, 1)" w l ls 13

set out "gumbelGrenz2_f3.eps"
print "plotting gumbelGrenz2_f3.eps"
plot\
 fE2(x,1) t "X_1" w l ls 1,\
 fGumbel(x,etamax(1),lambda) t "Gu(ln 1, 1)" w l ls 11,\
 fE2(x,2) t "max(X_1,X_2)" w l ls 3,\
 fGumbel(x,etamax(2),lambda) t "Gu(ln 2, 1)" w l ls 13,\
 fE2(x,10) t "max(X_1, ..., X_{10})" w l ls 5,\
 fGumbel(x,etamax(10),lambda) t "Gu(ln 3, 1)" w l ls 15


