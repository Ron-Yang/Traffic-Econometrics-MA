


# siehe ~/info/gnuplot, ~/info/gnuTemplate.gnu

##########################################################
#geordnet nach hue (ps)
# set style line: linetype lt, point type pt

set encoding iso_8859_1   # dann äöüßÄÖÜ durch woertl. Eingabe korrekt


set style line 1 lt 7 lw 4 pt 1 ps 1.5 #schwarz, plus sign
                               # (screen: tausche mit lt4=braun); 
set style line 2 lt 1 lw 4 pt 4 ps 1.5 #rot, open box
set style line 3 lt 8 lw 4 pt 5 ps 1.5 #blassrot, closed square
set style line 4 lt 6 lw 4 pt 6 ps 1.5 #gelb, open circle
                               #  (screen: tausche mit lt 1=orange-ocker)
set style line 5 lt 2 lw 4 pt 8 ps 1.5 #gruen, open triangle
set style line 6 lt 5 lw 4 pt 9 ps 1.5 #blasstuerkisblau, closed triangle
set style line 7 lt 3 lw 4 pt 10 ps 1.5 #blau, upside-down open triangle
set style line 8 lt 4 lw 4 pt 11 ps 1.5 #lila, upside-down closed triangle
set style line 9 lt 7 lw 4 pt 12 ps 1.5 #blau, upside-down open diamond
set style line 10 lt 1 lw 4 pt 13 ps 1.5 #lila, upside-down closed diamond


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

#gamma(n+1)=n! for integer n>=0

studnorm(nu)  = gamma(0.5*(nu+1)) / (sqrt(nu*pi)*gamma(0.5*nu))
student(x,nu) =  studnorm(nu) / (1.+x**2/nu)**(0.5*(nu+1))

chinorm(nu)   = 1./(2**(0.5*nu)*gamma(0.5*nu))
chi2(x,nu)    =  chinorm(nu) *exp(-0.5*x) * x**(0.5*nu-1)

fac(x)=gamma(x+1)

#erlang(1,x)=standard Poisson distribution

b_erlang(n)=gamma(n+1) /gamma(n)
a_erlang(n)=(b_erlang(n))**(n)/gamma(n)
erlang(n,x)=a_erlang(n)*x**(n-1)*exp(-b_erlang(n)*x)

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x

############### Beispiele fuer Output ####################


set samples 200

#############################################
# Input
#############################################

n=2500.
theta1=0.4       # Anteil Schicht 1 "Studenten"
theta2=1.-theta1 # Anteil Schicht 2 "Nicht-Studenten"
mu1=0.125        # Erwartungswert Anteil MIV-Bevorzugung Studis
mu2=0.75

sxx1=mu1*(1.-mu1)
sxx2=mu2*(1.-mu2)
sx1=sqrt(sxx1)
sx2=sqrt(sxx2)
print "mu1=Anteil MIV Bevorzugung Studis=", mu1, " sxx1=",sxx1
print "mu2=Anteil MIV Bevorzugung Sonst=", mu2, " sxx2=",sxx2


#############################################
print "\nX=MIV-Bevorzugung: Zufallsauswahl ohne Schichtstruktur"
#############################################

mu=theta1*mu1+theta2*mu2
sxxquerZufall=1./n *(theta1*(sxx1+(mu1-mu)**2)\
                     + theta2*(sxx2+(mu2-mu)**2))
sxquerZufall=sqrt(sxxquerZufall)
print "mu=",mu
print "sxxquerZufall=",sxxquerZufall
print "sxquerZufall=",sxquerZufall

print "Test direkt: sxxquerZufall=",(mu*(1.-mu)/n)


#############################################
print "\nX=MIV-Bevorzugung: Proportionale Schichtauswahl"
#############################################


sxxquerSchicht=1./n *(theta1*sxx1+theta2*sxx2)
sxquerSchicht=sqrt(sxxquerSchicht)
print "sxxquerSchicht=",sxxquerSchicht
print "sxquerSchicht=",sxquerSchicht



#############################################
print "\nX=MIV-Bevorzugung: Zufallsauswahl mit Entzerrung gemaess Schichtstruktur"
#############################################

n1=1050.   # statt n*theta1
n2=n-n1    # statt n*theta2

E1=theta1/(n1/n)
E2=theta2/(n2/n)

sxxquerEntzerr=1./n *(E1*theta1*sxx1+E2*theta2*sxx2)
sxquerEntzerr=sqrt(sxxquerEntzerr)
print "Zufaellig gefundene Schichtanteile in SP: n1=",n1," n2=",n2
print "Entzerrungsfaktoren: E1=",E1," E2=",E2
print "sxxquerEntzerr=",sxxquerEntzerr
print "sxquerEntzerr=",sxquerEntzerr


#############################################
print "\nX=MIV-Bevorzugung: Optimale Schichtauswahl"
#############################################

denom=theta1*sx1+theta2*sx2
n1opt=n*theta1*sx1/denom
n2opt=n*theta2*sx2/denom
E1=theta1/(n1opt/n)
E2=theta2/(n2opt/n)

sxxquerOpt=1./n *(E1*theta1*sxx1+E2*theta2*sxx2)
sxquerOpt=sqrt(sxxquerOpt)
print "Optimale Schichtanteile in SP:"
print "  n1opt=",n1opt," n2opt=",n2opt
print "Entzerrungsfaktoren der optimalen Anteile:"
print "  E1=",E1," E2=",E2
print "sxxquerOpt=",sxxquerOpt
print "sxquerOpt=",sxquerOpt
#print "Probe: sxxquerOpt=",(1./n*denom**2)






#############################################
# Input2
#############################################

n=2000.
theta1=0.5       # Anteil Schicht 1 "Frauen"
theta2=1.-theta1 # Anteil Schicht 2
mu1=2.5           # Erwartungswert Mobilitaet (Wege/Tag)
mu2=3.7 #3.7

sx1=3.  #3
sx2=1.  #1
sxx1=sx1**2
sxx2=sx2**2
print "mu1=Mobilitaet (Wege/Tag) Frauen=", mu1, " sxx1=",sxx1
print "mu2=", mu2, " sxx2=",sxx2


#############################################
print "\nX=Mobilitaet: Zufallsauswahl ohne Schichtstruktur"
#############################################

mu=theta1*mu1+theta2*mu2
sxxquerZufall=1./n *(theta1*(sxx1+(mu1-mu)**2)\
                     + theta2*(sxx2+(mu2-mu)**2))
sxquerZufall=sqrt(sxxquerZufall)
print "mu=",mu
print "sxxquerZufall=",sxxquerZufall
print "sxquerZufall=",sxquerZufall


#############################################
print "\nX=Mobilitaet: Proportionale Schichtauswahl"
#############################################


sxxquerSchicht=1./n *(theta1*sxx1+theta2*sxx2)
sxquerSchicht=sqrt(sxxquerSchicht)
print "sxxquerSchicht=",sxxquerSchicht
print "sxquerSchicht=",sxquerSchicht



#############################################
print "\nX=Mobilitaet: Zufallsauswahl mit Entzerrung gemaess Schichtstruktur"
#############################################

n1=1050.   # statt n*theta1
n2=n-n1    # statt n*theta2

E1=theta1/(n1/n)
E2=theta2/(n2/n)

sxxquerEntzerr=1./n *(E1*theta1*sxx1+E2*theta2*sxx2)
sxquerEntzerr=sqrt(sxxquerEntzerr)
print "Zufaellig gefundene Schichtanteile in SP: n1=",n1," n2=",n2
print "Entzerrungsfaktoren: E1=",E1," E2=",E2
print "sxxquerEntzerr=",sxxquerEntzerr
print "sxquerEntzerr=",sxquerEntzerr


#############################################
print "\nX=Mobilitaet: Optimale Schichtauswahl"
#############################################

denom=theta1*sx1+theta2*sx2
n1opt=n*theta1*sx1/denom
n2opt=n*theta2*sx2/denom
E1=theta1/(n1opt/n)
E2=theta2/(n2opt/n)

sxxquerOpt=1./n *(E1*theta1*sxx1+E2*theta2*sxx2)
sxquerOpt=sqrt(sxxquerOpt)
print "Optimale Schichtanteile in SP:"
print "  n1opt=",n1opt," n2opt=",n2opt
print "Entzerrungsfaktoren der optimalen Anteile:"
print "  E1=",E1," E2=",E2
print "sxxquerOpt=",sxxquerOpt
print "sxquerOpt=",sxquerOpt


#############################################
print "\nX=Mobilitaet: Zufallsauswahl mit optimalen Gewichten
#############################################

n1=1050. # 1050  # statt n*theta1
n2=n-n1    # statt n*theta2

f1=n1/n
f2=n2/n
denom=f1*mu1**2/sxx1+f2*mu2**2/sxx2

E1=mu*mu1/sxx1/ denom
E2=mu*mu2/sxx2/ denom

#E1=theta1/(n1/n)
#E2=theta2/(n2/n)
#E1=0.3

sxxquerZufOpt=1./n *(E1**2*f1*sxx1+E2**2*f2*sxx2)
sxquerZufOpt=sqrt(sxxquerZufOpt)
print "Zufaellig gefundene Schichtanteile in SP: n1=",n1," n2=",n2
print "f1=",f1," f2=",f2
print "Optimale Gewichtungsfaktoren: E1=",E1," E2=",E2
print "sxxquerZufOpt=",sxxquerZufOpt
print "sxquerZufOpt=",sxquerZufOpt




#############################################
#############################################
print "\nErsparnis repraesentativ und optimal geschichtet gegenueber Zufasll"
#############################################

print "2 Schichten; Referenz: sigma1=1, mu1=0"
sxx1=1. #immer
mu1=0.  #immer


theta1_init=0.5     #einstellen!
sxx2_init  =9     #einstellen!
mu2_init   =1        #einstellen!



#############################################



theta1 = theta1_init;
sxx2   = sxx2_init;
mu2    = mu2_init;


muGes(theta1,mu1,mu2)=theta1*mu1+(1-theta1)*mu2

#### Zufall

ns2xquerZuf(theta1, mu1,mu2, sxx1,sxx2)\
  =theta1*(sxx1+(mu1-muGes(theta1,mu1,mu2))**2)\
 + (1.-theta1)*(sxx2+(mu2-muGes(theta1,mu1,mu2))**2)

### Proport. Schichtung

ns2xquerProp(theta1, sxx1, sxx2)\
  =theta1*sxx1+(1.-theta1)*sxx2

### Opt. Schichtung

denom(theta1, sxx1, sxx2)=theta1*sqrt(sxx1)+(1.-theta1)*sqrt(sxx2)

f1opt(theta1, sxx1, sxx2)=theta1*sqrt(sxx1)/denom(theta1, sxx1, sxx2)
f2opt(theta1, sxx1, sxx2)=(1.-theta1)*sqrt(sxx2)/denom(theta1, sxx1, sxx2)
ns2xquerOpt(theta1, sxx1, sxx2) = denom(theta1, sxx1, sxx2)**2


#############################################

set term post eps enhanced color dashed "Helvetica" 20
set out "geschichtet1.eps"
print " plotting geschichtet1.eps ..."

set xlabel "{/Symbol m}_2"
set ylabel "n {/Symbol s}^2_{xquer}"
set key screen 0.8,0.3


plot[mu2=0:4]\
 mu2, ns2xquerZuf(theta1, mu1, mu2, sxx1, sxx2) t "Zufall" w l ls 1,\
 mu2, ns2xquerProp(theta1, sxx1, sxx2) t "Prop. Schichtung" w l ls 2,\
 mu2, ns2xquerOpt(theta1, sxx1, sxx2) t "Opt. Schichtung" w l ls 3

#############################################

set out "geschichtet2.eps"
print " plotting geschichtet2.eps ..."

sxx2=sxx2_init
mu2=mu2_init

set xlabel "{/Symbol s}_2^2 / {/Symbol s}_1^2 "
set ylabel "{/Symbol s}^2_{xquer}/{/Symbol s}^2_{xquer,Zufall}"

plot[sxx2=0.05:3]\
 sxx2, 1 t "Zufall=1" w l ls 1,\
 sxx2, ns2xquerProp(theta1, sxx1, sxx2)/ns2xquerZuf(theta1, mu1, mu2, sxx1, sxx2)\
 t "Prop. Schichtung" w l ls 2,\
 sxx2, ns2xquerOpt(theta1, sxx1, sxx2)/ns2xquerZuf(theta1, mu1, mu2, sxx1, sxx2)\
 t "Opt. Schichtung" w l ls 3



#############################################

set out "geschichtet3.eps"
print " plotting geschichtet3.eps ..."

set key screen 0.8,0.9

#sxx2=sxx2_init;
mu2=mu2_init;
set ylabel "Optimaler Stichprobenanteil f_1 ({/Symbol q}_1=0.5)"

plot[sxx2=0.05:3]\
 sxx2, 0.5 t "Anteil in GG=Prop. Schichtung" w l ls 1,\
 sxx2, f1opt(theta1, sxx1, sxx2) t "Opt. Schichtung" w l ls 7


#############################################

set out "geschichtet4.eps"
print " plotting geschichtet4.eps ..."
set xlabel "Anteil f_1 in Stichprobe"
set ylabel "{/Symbol s}^2_{xquer}/{/Symbol s}^2_{xquer,opt}"
set yrange [1:5]

ns2xquerAllg(f1, theta1, sxx1, sxx2)=theta1**2*sxx1/f1+(1.-theta1)**2*sxx2/(1.-f1)
set key screen 0.7,0.9

sxx2=sxx2_init;
mu2=mu2_init;

var2a=1.
var2b=9.
plot[f1=0.001:0.999]\
 f1opt(theta1, sxx1, var2a), 1+2*f1\
 t "Optimale Schichtung {/Symbol s}_2={/Symbol s}_1" w l ls 6,\
 f1opt(theta1, sxx1, var2b), 1+2*f1\
 t "Optimale Schichtung {/Symbol s}_2=3 {/Symbol s}_1" w l ls 3,\
 theta1, 1+2*f1 t "Proportionale Schichtung" w l ls 1,\
 f1, ns2xquerAllg(f1, theta1, sxx1, var2a)/ns2xquerOpt(theta1, sxx1, var2a)\
 t "{/Symbol s}_2={/Symbol s}_1" w l ls 7,\
 f1, ns2xquerAllg(f1, theta1, sxx1, var2b)/ns2xquerOpt(theta1, sxx1, var2b)\
 t "{/Symbol s}_2=3{/Symbol s}_1" w l ls 2


##########################################################
# Spezialfall: Merkmale= Anteilswerte:
##########################################################

theta1 = 0.5; #Anteil der SCHICHT 1

A1=0.05   # Anteil fuer Alternative 1 in GG Schicht 1
A2=0.5   # Anteil fuer Alternative 1 in GG Schicht 2

B1=0.2   # andere Alterbnativen-Anteile
B2=0.8   # Anteil fuer Alternative 1 in GG Schicht 2


mu(A)=A
sxx(A) =A*(1.-A)


#############################

set out "geschichtet1A.eps"
print " plotting geschichtet1A.eps ..."

set xlabel "{/Symbol q}_1"
set ylabel "n V ({/Symbol m}Dach)"
set label 1 "{/Symbol m}_1=5%, {/Symbol m}_2=50%" at screen 0.25,0.5
set key screen 0.55,0.40
set auto y


plot[theta1=0:1]\
 theta1, ns2xquerZuf(theta1, mu(A1), mu(A2), sxx(A1),sxx(A2))\
 t "Zufall " w l ls 1,\
 theta1, ns2xquerProp(theta1, sxx(A1), sxx(A2))\
 t "Prop. Schichtung " w l ls 2,\
 theta1, ns2xquerOpt(theta1, sxx(A1), sxx(A2))\
 t "Opt. Schichtung " w l ls 3


set out "geschichtet1B.eps"
print " plotting geschichtet1B.eps ..."

set auto y

set label 1 "{/Symbol m}_1=20%, {/Symbol m}_2=80%" at screen 0.4,0.5
set key screen 0.55,0.40
plot[theta1=0:1]\
 theta1, ns2xquerZuf(theta1, mu(B1), mu(B2), sxx(B1),sxx(B2))\
 t "Zufall symm" w l ls 7,\
 theta1, ns2xquerProp(theta1, sxx(B1), sxx(B2))\
 t "Prop. Schichtung symm" w l ls 6,\
 theta1, ns2xquerOpt(theta1, sxx(B1), sxx(B2))\
 t "Opt. Schichtung symm" w l ls 5

#############################

set out "geschichtet2A.eps"
print " plotting geschichtet2A.eps ..."

set xlabel "{/Symbol q}_1"
set ylabel "Stichpoben-Schichtanteil f_1"
set label 1 "{/Symbol m}_1=5%, {/Symbol m}_2=50%" at screen 0.7,0.4
set key screen 0.8,0.3
set auto y

plot[theta1=0:1]\
 theta1, theta1\
 t "Prop. Schichtung" w l ls 1,\
 theta1, f1opt(theta1, sxx(A1), sxx(A2))\
 t "Opt. Schichtung" w l ls 3

#############################

set out "geschichtet2B.eps"
print " plotting geschichtet2B.eps ..."

set xlabel "{/Symbol q}_1"
set ylabel "Stichpoben-Schichtanteil f_1"
set label 1 "{/Symbol m}_1=20%, {/Symbol m}_2=80%" at screen 0.8,0.5
set key screen 0.8,0.3
set auto y


plot[theta1=0:1]\
 theta1, theta1\
 t "Prop. Schichtung" w l ls 1,\
 theta1, f1opt(theta1, sxx(B1), sxx(B2))\
 t "Opt. Schichtung symm" w l ls 7





#############################

set out "geschichtet4A.eps"
print " plotting geschichtet4A.eps ..."

set xlabel "Stichpoben-Schichtanteil f_1"
set ylabel "n V ({/Symbol m}Dach)"

set label 1 "{/Symbol m}_1=5%, {/Symbol m}_2=50%" at screen 0.5,0.9
set key screen 0.65,0.82
set yrange [0:0.5]

plot[f1=0:1]\
 0.5, 0.2*f1 t "Prop. Aufteilung {/Symbol q}_1=0.5" w l ls 1,\
 f1opt(0.5, sxx(A1), sxx(A2)), 0.2*f1\
 t "Optimale Aufteilung" w l ls 3,\
 f1, ns2xquerAllg(f1, 0.5, sxx(A1), sxx(A2))\
 t "Allgemeine Schichtung" w l ls 2
# 0.8, 0.2*f1 t "Prop. Aufteilung {/Symbol q}_1=0.8" w l ls 1,\
# f1, ns2xquerAllg(f1, 0.8, sxx(A1), sxx(A2))\
# t "Allgemeine Schichtung, \"\"" w l ls 7,\
# f1opt(0.8, sxx(A1), sxx(A2)), 0.2*f1\
# t "Optimale Aufteilung \"\"" w l ls 6







