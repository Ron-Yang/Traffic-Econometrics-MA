
set encoding iso_8859_1   # dann %G�������%@ durch woertl. Eingabe korrekt
#\304     "A
#\326     "O
#\334     "U
#\344     "a
#\366     "o
#\374     "u
#\337     &szlig;

# Siehe auch
# figsOekonometrie/geschichtet_jun06.gnu 


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

#######################################

set term post eps enhanced color solid "Helvetica" 24
set out "IOMdyn1.eps"
print "plotting IOMdyn1.eps"
a=0.2
ya=100.
ye=120.
tau=1.
b=1./(1-a)
y(t)=(t>=0) ? ye : ya
x(t)=(t>=0) ? b*(ye+(ya-ye)*exp(-t/(b*tau))) : b*ya
set noparam
set xlabel "Zeit (Einheiten von {/Symbol t})"
set ylabel "Angebot und Nachfrage (Geldeinheiten)"
set samples 100
set key at screen 0.8,0.7

plot [t=-1:4]\
 y(t) w l ls 7 t "Nachfrage y(t)",\
 x(t)/b w l ls 2 t "Angebot x(t)/b",\
 x(t) w l ls 3 t "Gesamtproduktion x(t)

#######################################

set term post eps enhanced color solid "Helvetica" 22
set out "IOMdyn2.eps"
print "plotting IOMdyn2.eps"
set key at screen 0.84,0.62

set xlabel "Zeit (Jahre)"
set ylabel "Angebot/Nachfrage bezogen auf x_2(0)

y2(t)=(t>=0) ? 2: 1

x1supply(t)=(t>=0) ? 1.-0.5*exp(-0.5*t) : 0.5
x2supply(t)=(t>=0) ? 2.  - 2./3.*exp(-0.5*t) - 1./3.*exp(-2*t) : 1
x2demand(t)=(t>=0) ? 2.  - exp(-2*t)  : 1
plot [t=-1:4]\
 y2(t) w l ls 7 t "Nachfrage Solarzellen y_2(t)",\
 x2supply(t) w l ls 2 t "Solarzellen x_2(t)",\
 x1supply(t) w l ls 3 t "Silizium x_1(t)",\
 x2demand(t) w l ls 5 t "x_2(t) ohne Siliziummangel"


#######################################

set term post eps enhanced color solid "Helvetica" 22
set out "supplyChain1.eps"
print "plotting supplyChain1.eps"
set key at screen 0.84,0.92

set xlabel "Zeit (Jahre)"
set ylabel "Materialfluss, Lagerbestand"
set yrange [0:3.5]
y2(t)=(t>=0) ? 2: 1
omega=1.
x2(t)=(t>=0) ? -sin(omega*t) : 0
x1(t)=(t>=0) ? 2-cos(omega*t) : 1

plot [t=-1:8.5]\
 y2(t) w l ls 7 t "Nachfrage  y_2(t)",\
 1+x2(t) w l ls 2 t "Lagerbestand/Sollbestand 1+x_2(t)",\
 x1(t) w l ls 3 t "Produktion x_1(t)"

##############################
set out "bullwhip.eps"
print "plotting bullwhip.eps"

omega=2*pi
tdelay2=9./7
tdelay3=5./7
beta2=2.*pi/tdelay2
beta3=2.*pi/tdelay3
y0=1
delta_y=0.1
delta_x23=delta_y/(1-(omega/beta3)**2)
delta_x12=delta_x23/(1-(omega/beta2)**2)

y(t)=y0+delta_y*sin(omega*t)
x23(t)=y0+delta_x23*sin(omega*t)
x12(t)=y0+delta_x12*sin(omega*t)
set xlabel "Zeit (Wochen)"
set ylabel "Bierumsatz (a.u.)
tmin=0
tmax=3*2*pi/omega
set yrange [y0-abs(1.1*delta_x12):y0+abs(1.55*delta_x12)]
plot [t=tmin:tmax]\
 y(t) t "Nachfrage" w l ls 2,\
 x23(t) t "Bierfluss Gro\337handel-Kneipe" w l ls 3,\
 x12(t) t "Bierfluss Brauerei-Gro\337handel" w l ls 7
















