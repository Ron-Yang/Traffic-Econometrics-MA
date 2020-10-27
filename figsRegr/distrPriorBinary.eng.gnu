

set style line 99 lt 1 lw 3 pt 4 ps 1.5 linecolor rgb "#1100EE" #blau, solid, open box

set style line 1 lt 1 lw 3 pt 7 ps 1.9  lc rgb "#000000" #schwarz,solid,bullet
set style line 2 lt 1 lw 3 pt 2 ps 1.5  lc rgb "#CC0022" #rot, solid, Kreuz
set style line 3 lt 8 lw 3 pt 4 ps 1.2  lc rgb "#FF3300"#orange, offenes Quadrat
set style line 4 lt 6 lw 3 pt 4 ps 1.5  lc rgb "#FFAA00"  #gelb, offenes Quadrat
set style line 5 lt 1 lw 3 pt 5 ps 1.5  lc rgb "#00DD22"  #gruen,solid,ClosedBox
set style line 6 lt 5 lw 3 pt 4 ps 1.5  lc rgb "#00AAAA" #offenes Quadrat
set style line 7 lt 3 lw 3 pt 4 ps 2.0  lc rgb "#1100FF"  #blau,gepunktet,offenes Quadrat
set style line 8 lt 4 lw 3 pt 8 ps 1.5  lc rgb "#220088"
set style line 9 lt 7 lw 3 pt 9 ps 1.5  lc rgb "#999999"  #grau, aufr. gschl. Dreieck

set style line 11 lt 1 lw 6 pt 7 ps 1.9  lc rgb "#000000" #schwarz,solid,bullet
set style line 12 lt 1 lw 6 pt 2 ps 1.5  lc rgb "#CC0022" #rot, dash, Kreuz
set style line 13 lt 8 lw 6 pt 4 ps 1.2  lc rgb "#FF3300"#orange, offenes Quadrat
set style line 14 lt 6 lw 6 pt 4 ps 1.5  lc rgb "#FFAA00"  #gelb, offenes Quadrat
set style line 15 lt 1 lw 6 pt 5 ps 1.5  lc rgb "#00DD22"  #gruen,solid,ClosedBox
set style line 16 lt 5 lw 6 pt 7 ps 1.5  lc rgb "#00AAAA" #offener Kreis
set style line 17 lt 1 lw 6 pt 7 ps 1.5  lc rgb "#1100FF"  #blau,solid,Bullet
set style line 18 lt 4 lw 6 pt 8 ps 1.5  lc rgb "#220088"
set style line 19 lt 7 lw 6 pt 9 ps 1.5  lc rgb "#999999"  #grau, aufr. gschl. Dreieck

#Sinnvolle point typen (pt)
# 1=Plus,2=Kreuz,4=openQuare,5=closedSquare, 6=openCirc,7=closedCirc,
# 9-11=triangles, 12-13=Rauten



############### Beispiele fuer Funktionen ####################

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x
mod(x,interval)=x-(floor(x/interval)*interval) # x%interval for reals

fac(x)=gamma(x+1)
lnfac(x)=(x<50) ? log(gamma(x+1)) : x*log(x)-x


 ############################
 # Binomialverteilung und Derivate
 ############################

f_binom(n,theta,x)=fac(n)/(fac(x)*fac(n-x)) * theta**x*(1.-theta)**(n-x)

 ############################
 # Normierte Gauss-Verteilungen und Derivate
 ############################


gauss(mu,sigma,x) = exp(-(x-mu)**2/(2.* sigma**2)) / (sigma * sqrt(2*pi))
gaussDens(mu,sig2,x) = exp(-(x-mu)**2/(2.* sig2)) / sqrt(2*pi*sig2) #Alt.

gauss2d(x,y,sx,sy,r) \
 = exp(-(0.5/(1-r*r))*((x/sx)**2+(y/sy)**2-2*r*x*y/(sx*sy)))\
  / (2*pi*sx*sy*sqrt(1-r*r))
gauss2d_invCovMatrix(x1,x2,invs11, invs12,invs22)\
 = sqrt(invs11*invs22-invs12**2) /  (2*pi)\
 * exp(-0.5*(x1*invs11*x1+2*x1*invs12*x2+x2*invs22*x2))


  # integral der Standardnormalverteilung
phi(x)            =norm(x)
                
 # Quantil der Standardnormalverteilung
phiQuantil(q)=invnorm(q)


 ############################
 # Student-t-Verteilung  und Derivate, nu=Zahl der FG
 ############################

# gamma(n+1)=n! for integer n>=0
# nu=Zahl der FG

studnorm(nu)  = gamma(0.5*(nu+1)) / (sqrt(nu*pi)*gamma(0.5*nu))
student(x,nu) =  studnorm(nu) / (1.+x**2/nu)**(0.5*(nu+1))
studentCum(x,nu) = ibeta(0.5*nu,0.5*nu,0.5*(x+sqrt(x**2+nu))/sqrt(x**2+nu))

#Iteration; Startwert: phiQuantil(p)
tQuantil_it(t,p,nu)=t-(studentCum(t,nu)-p)/student(t,nu)
tQuantil(p,nu)=tQuantil_it(tQuantil_it(tQuantil_it(tQuantil_it(tQuantil_it(\
  tQuantil_it(tQuantil_it(tQuantil_it(tQuantil_it(tQuantil_it(\
  phiQuantil(p),p,nu), p,nu), p,nu), p,nu), p,nu), p,nu),\
  p,nu), p,nu), p,nu), p,nu)


 

##################################################################
# Bayes Approach:
# beta    = {0 with prob PH0 and d with prob 1-PH0} diskret binaer verteilt
# hatbeta = beta+U, U=sigb*Z, Z sim N(0,1) gaussvert Schaetzer (Dichte h(u))
# b       = konkrete Realisierung  von hatbeta ("Beobachtung" bzw. "Messung")
# B       = [b-delta/2, b+delta/2] Beobachtungsintervall bzw. "Beobachtung"
# H0      = {beta=0} oder {beta in [-infty,0]}, so dass Prior-H0-Wahrsch. PH0
# Grenzwert der H0-Menge: beta0=0
# Extremmenge E_b={hatbeta: hatbeta \ge b}
# p-Wert p_b=P(E_b|beta=beta0)=1-Phi(b/sigb)
# bzw. b als Funktion von p: b(p)=sigb Phi^{-1}(1-p)   (1)

# Bayes:  P(H0|B)=P(B|H0)*P(H0)/P(B)
# P(H0)=PH0
# P(B|H0)=delta*h(b)  (geht direkt, da H0 Punktmenge, sonst nur P(H0
#                     \cap B) ausrechenbar)
# P(B)=PH0*P(B|H0)+(1-PH0)*P(B|\bar{H_0})
# =delta*(PH0*h(b)+(1-PH0)*h(b-d))
 

##################################################################

b(p,sigb)=sigb*phiQuantil(1-p)
PH0post(b, sigb, d, PH0)=PH0*gaussDens(0,sigb**2,b)\
 / (PH0*gaussDens(0,sigb**2,b)+(1-PH0)*gaussDens(d,sigb**2,b))

print "PH0post(20,10,50,0.8)=",PH0post(20,10,50,0.8)


##################################################################
#set term post eps enhanced color solid "Helvetica" 16
set term pngcairo enhanced color notransparent crop font "Helvetica, 14"

set noparam
##################################################################

##################################################################
set out "PH0_PriorBinaryProb05_b.png"
print   "plotting PH0_PriorBinaryProb05_b.png"
##################################################################

set key left bottom box
set label 1 at screen 0.14,0.44
set label 2 at screen 0.14,0.40
set label 3 at screen 0.77,0.90

PH0=0.5
d=1.

strTitle(PH0)=sprintf("A-priori probability\
 P(H_0)=P(on freeway)=%1.2f", PH0)
set xlabel "GPS measurement hat(y)/d (0: freeway, 1: parallel road)" offset 0,0.5
set ylabel "A posteriori probability P(H_0|hat(y))" offset 1,0
#set label 1 "Binary variable  {/Symbol b}=0 (H_0) und {/Symbol b}=d"
set label 2 "Stddev  of estimation error : {/Symbol s}"
set label 3 sprintf("P(H_0)=%1.2f", PH0)
#set title strTitle(PH0)


plot[r=0:1]\
 PH0post(d*r, 0.2, d, PH0) t "{/Symbol s}/d=0.2" w l ls 1,\
 PH0post(d*r, 0.3, d, PH0) t "{/Symbol s}/d=0.3" w l ls 2,\
 PH0post(d*r, 0.5, d, PH0) t "{/Symbol s}/d=0.5" w l ls 3,\
 PH0post(d*r, 1.0, d, PH0) t "{/Symbol s}/d=1.0" w l ls 5

##################################################################
set out "PH0_PriorBinaryProb08_b.png"
print   "plotting PH0_PriorBinaryProb08_b.png"
##################################################################
PH0=0.8
set label 3 sprintf("P(H_0)=%1.2f", PH0)
#set title strTitle(PH0)
plot[r=0:1]\
 PH0post(d*r, 0.2, d, PH0) t "{/Symbol s}/d=0.2" w l ls 1,\
 PH0post(d*r, 0.3, d, PH0) t "{/Symbol s}/d=0.3" w l ls 2,\
 PH0post(d*r, 0.5, d, PH0) t "{/Symbol s}/d=0.5" w l ls 3,\
 PH0post(d*r, 1.0, d, PH0) t "{/Symbol s}/d=1.0" w l ls 5


##################################################################
set out "PH0_PriorBinaryProb05_p.png"
print   "plotting PH0_PriorBinaryProb05_p.png"
##################################################################


set noparam
set key right bottom box

set xlabel "p value for test of H_0"
#set ylabel "A-Posteriori Wahrscheinlichkeit von H0: \"Autobahn\""
set label 1 at screen 0.50,0.40
set label 2 at screen 0.50,0.35
#set label 1 "Moegliche wahre Werte: {/Symbol b}=0 (H_0) und {/Symbol b}=d"
#set label 2 "Standardabweichung der Messfehler: {/Symbol s}_b"

PH0=0.5
set label 3 sprintf("P(H_0)=%1.2f", PH0)
sigb=10.
d=50. # variabel im Plot

plot[p=0:0.5]\
 PH0post(b(p,sigb), sigb, 50, PH0) t "{/Symbol s}/d=0.2" w l ls 1,\
 PH0post(b(p,sigb), sigb, 33, PH0) t "{/Symbol s}/d=0.3" w l ls 2,\
 PH0post(b(p,sigb), sigb, 20, PH0) t "{/Symbol s}/d=0.5" w l ls 3,\
 PH0post(b(p,sigb), sigb, 10, PH0) t "{/Symbol s}/d=1.0" w l ls 5,\
 PH0post(b(p,sigb), sigb, 5, PH0)  t "{/Symbol s}/d=2.0" w l ls 6,\
 PH0post(b(p,sigb), sigb, 2, PH0)  t "{/Symbol s}/d=5.0" w l ls 7

##################################################################
set out "PH0_PriorBinaryProb08_p.png"
print   "plotting PH0_PriorBinaryProb08_p.png"
##################################################################
PH0=0.8
set label 3 sprintf("P(H_0)=%1.2f", PH0)
sigb=10.
d=50.

plot[p=0:0.5]\
 PH0post(b(p,sigb), sigb, 50, PH0) t "{/Symbol s}/d=0.2" w l ls 1,\
 PH0post(b(p,sigb), sigb, 33, PH0) t "{/Symbol s}/d=0.3" w l ls 2,\
 PH0post(b(p,sigb), sigb, 20, PH0) t "{/Symbol s}/d=0.5" w l ls 3,\
 PH0post(b(p,sigb), sigb, 10, PH0) t "{/Symbol s}/d=1.0" w l ls 5,\
 PH0post(b(p,sigb), sigb, 5, PH0)  t "{/Symbol s}/d=2.0" w l ls 6,\
 PH0post(b(p,sigb), sigb, 2, PH0)  t "{/Symbol s}/d=5.0" w l ls 7

print "\nIn-Text-Aufgabe Map-Matching:"
print "1-phi(2)=",1-phi(2)
pH0=0.8
pBH0=gaussDens(0,100,20)
pBH0bar=gaussDens(0,100,-30)
pB=pH0*pBH0+(1-pH0)*pBH0bar
pH0B=pH0*pBH0/pB
print "pBH0=delta*",pBH0
print "pBH0bar=delta*",pBH0bar
print "pB=delta*",pB
print "pH0*pBH0/pB=",pH0*pBH0/pB
