
# gnuplot4-Syntax!

# von  ~/info/gnuTemplate.gnu
# siehe auch
# ~/info/gnuplot,   
# ~/info/gnuTemplate42.gnu,  
# ~/info/gnuColoredContour/*.gnu


##########################################################
#geordnet nach hue (ps)
# set style line: linetype lt, point type pt

set encoding iso_8859_1   # dann äöüßÄÖÜ durch woertl. Eingabe korrekt

# (Punkttypen Feb. 2007)
# BUG: Manchmal
#     obskure Abhaengigkeit der Symbole von Reihenfolge bzw linestyles

set style line 1 lt 7 lw 10 pt 1 ps 1.5 #schwarz, plus sign
set style line 2 lt 1 lw 6 pt 2 ps 1.5 #rot, Kreuz
set style line 3 lt 8 lw 2 pt 3 ps 1.5 #blassrot, offener Stern
set style line 4 lt 6 lw 2 pt 4 ps 1.5 #gelb, offenes Quadrat
set style line 5 lt 2 lw 2 pt 5 ps 1.5 #gruen, geschl. Quadrat
set style line 6 lt 5 lw 2 pt 6 ps 1.5 #blasstuerkisblau, offener Kreis
set style line 7 lt 3 lw 2 pt 7 ps 1.5 #blau, Bullet=geschloss. Kreis!
set style line 8 lt 4 lw 2 pt 8 ps 1.5 #lila, aufrechtes geschloss. Dreieck
set style line 9 lt 7 lw 2 pt 9 ps 1.5 #schwarz, aufrechtes geschl. Dreieck
set style line 10 lt 1 lw 2 pt 10 ps 1.5 #rot, upside-down offenes Dreieck
set style line 11 lt 7 lw 10 pt 11 ps 1.5 #schwarz, upside-down geschl. Dreieck
set style line 12 lt 1 lw 2 pt 12 ps 1.5 #rot, offene Raute
set style line 13 lt 8 lw 2 pt 13 ps 1.5 #blassrot, geschl. Raute
set style line 14 lt 6 lw 2 pt 14 ps 1.5 #gelb, "Sonderzeichen" ...
set style line 15 lt 2 lw 2 pt 15 ps 1.5 #gruen, geschl. Fuenfeck
set style line 16 lt 5 lw 2 pt 16 ps 1.5 #blasstuerkisblau, 
set style line 17 lt 3 lw 2 pt 17 ps 1.5 #blau, 
set style line 18 lt 4 lw 2 pt 18 ps 1.5 #lila,
set style line 19 lt 7 lw 2 pt 19 ps 1.5
set style line 20 lt 1 lw 2 pt 20 ps 1.5



############### Beispiele fuer Funktionen ####################

gauss(mu,sigma,x) = exp(-(x-mu)**2/(2* sigma**2)) / (sigma * sqrt(2*pi))

xlimited(x)       =(x<-5) ? -5 : ( (x>5) ? 5 : x)

  # integral der Standardnormalverteilung
phi(x)            =0.5*(erf(xlimited(x/sqrt(2)))+1) 
                
 # Quantilder Standardnormalverteilung
phiQuantil(q)=sqrt(2.)*inverf(-1.+2*q)
                


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

set term post eps enhanced color solid "Helvetica" 24
set out "limitedGrowth_eng.eps"
print "plotting limitedGrowth_eng.eps ..."

set xlabel "Year"
set ylabel "Car Penetration (%)"
y0=3. # Prozent in 1950
ys=60. # Saettigung in Prozent
t0=1950.
tau=10.
y(t,y0,t0,ys,tau)=ys/(1. + (ys/y0-1)*exp(-(t-t0)/tau))
set noparam
set xrange [1950:2020]
plot y(x,y0,t0,ys,tau) w l ls 2
