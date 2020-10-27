

set style line 99 lt 1 lw 3 pt 4 ps 1.5 linecolor rgb "#1100EE" #blau, solid, open box

set style line 1 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000" #schwarz,solid,bullet
set style line 2 lt 1 lw 2 pt 2 ps 1.5  lc rgb "#CC0022" #rot, solid, Kreuz
set style line 3 lt 8 lw 2 pt 4 ps 1.2  lc rgb "#FF3300"#orange, offenes Quadrat
set style line 4 lt 6 lw 2 pt 4 ps 1.5  lc rgb "#FFAA00"  #gelb, offenes Quadrat
set style line 5 lt 1 lw 2 pt 5 ps 1.5  lc rgb "#00DD22"  #gruen,solid,ClosedBox
set style line 6 lt 5 lw 2 pt 4 ps 1.5  lc rgb "#00AAAA" #offenes Quadrat
set style line 7 lt 3 lw 2 pt 4 ps 2.0  lc rgb "#1100FF"  #blau,gepunktet,offenes Quadrat
set style line 8 lt 4 lw 2 pt 8 ps 1.5  lc rgb "#220088"
set style line 9 lt 7 lw 2 pt 9 ps 1.5  lc rgb "#999999"  #grau, aufr. gschl. Dreieck

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


gauss(mu,sigma,x) = exp(-(x-mu)**2/(2* sigma**2)) / (sigma * sqrt(2*pi))
gaussDens(mu,sig2,x) = exp(-(x-mu)**2/(2* sig2)) / sqrt(2*pi*sig2) #Alt.

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
# beta \sim N(0, sigbeta^2) = zu schaetzende Groesse mit A-Priori
#      Gaussverteilung der Dichte g(beta) (o.E.d.A Erw-Wert=0)
# hatbeta=beta+U, U=sigb*Z, Z sim N(0,1) Schaetzer von beta
# U=sigb*Z unverzerrt-gaussverteilte Unschaerfe/Messfehler (allgemein Dichte h(u))
# b=konkrete Realisierung  von hatbeta ("Beobachtung" bzw. "Messung"
#      von beta) aus den Daten 
# Trick, um Bayes auf endl. Wahrsch. anzuwenden: Kleines eps
# B=[b-delta/2, b+delta/2] Beobachtungsintervall bzw. "Beobachtung"

# H0: beta \le beta0, so dass Prior-H0-Wahrsch. P(H_0)=p0=Phi(beta0/sigbeta) 
# Grenzwert der H0-Menge: beta=beta0
# Extremmenge E_b={hatbeta: hatbeta \ge b}
# p-Wert p_b=P(E_b|beta=beta0)=1-Phi((b-beta0)/sigb)
# bzw. b als Funktion von p: b(p)=sigb Phi^{-1}(1-p)+beta0    (1)

# Ziehung des wahren Wertes aus A-Priori-Verteilung und Messprozess U
# unabhaengig von beta mit (keine grosse Eischraenkung) bekannter Varianz sigb

# 1. Unabhaengigkeit => hatbeta=beta+U  mit Dichten g(beta) bzw. h(u)
# hat Dichte aus Faltungssatz: f(b)=\int(dbeta)g(beta)*h(b-beta)
# => P(B_b)=delta * f(b)   = delta * int(-infty,infty,dbeta) g(beta)*h(b-beta)   (2)

# 2. Wahrsch. dafuer, dass H0 UND  Z in D:

# P(H0 \cap B)=P(B|H0)*P(H0)=delta * int(-infty,beta0,dbeta) g(beta)*h(b-beta)    (3)

# 3. Bayes:  P(H0|B)=P(B|H0)*P(H0)/P(B) mit (2), (3) und Ausrechnen:

##################################################################
# P(H0|B)=P(H0|d)=P(p,beta0,sigbeta,sigb)= Phi((beta0-mu)/sigma), 
# mu=b * sigbeta**2/(sigbeta**2+sigb**2)
# sigma=sigbeta*sigb/sqrt(sigbeta**2+sigb**2)
# b=sigb Phi^{-1}(1-p)+beta0 
##################################################################
#
# mit (1): mu= (beta0+sigb Phi^{-1}(1-p)) * sigbeta**2/(sigbeta**2+sigb**2)
# 
# => Vereinfachung zu
#
# P(H0|p)=P(H0|B|B liefert p-Werte um p)
# =Phi(beta0/sigma-(beta0+sigb*sigbeta/sqrt(sigbeta**2+sigb**2)*Phi^{-1}(1-p)))
# special case for beta0=0:
# P(H0|p)=Phi(-sigbeta/sqrt(sigbeta**2+sigb**2)*Phi^{-1}(1-p))

##################################################################

sigma(sigbeta,sigb)=sigbeta*sigb/sqrt(sigbeta**2+sigb**2)
P_H0(p, beta0, sigbeta, sigb)=phi(beta0/sigma(sigbeta,sigb) \
 -(beta0/sigb+phiQuantil(1-p))*sigbeta/sqrt(sigbeta**2+sigb**2))

P_H0test(p,beta0,sigbeta,sigb)=phi((beta0-mu(p,beta0,sigbeta,sigb))/sigma(sigbeta,sigb))
mu(p,beta0,sigbeta,sigb)=b(p,beta0,sigb) * sigbeta**2/(sigbeta**2+sigb**2)
b(p,beta0,sigb)=beta0+sigb*phiQuantil(1-p)




##################################################################
#set term post eps enhanced color solid "Helvetica" 14
#set term post eps enhanced color dashed "Helvetica" 14
set term pngcairo enhanced color notransparent crop font "Helvetica, 12"

set noparam
set yrange [0:1]
##################################################################

##################################################################
set out "PH0_PriorGauss_beta0eq0.png"
print "plotting PH0_PriorGauss_beta0eq0.png"
##################################################################

beta0=0.
set key right top box
set label 1 at screen 0.65,0.68
set label 2 at screen 0.65,0.64
set label 3 at screen 0.50,0.90

strTitle(beta0)=sprintf("Gaussian a-priori distribution,\
 P(H_0)=%1.4f", phi(beta0))
set xlabel "p value" offset 0,0.4
set ylabel "Bayes a-posteriori probability P(H_0|data)"
set label 1 "stddev of true {/Symbol b} values: {/Symbol s}_{/Symbol b}"
set label 2 "stddev of the estimator: {/Symbol s}_b"
set label 3 sprintf("P(H_0)=%1.2f", phi(beta0))
#print "beta0=",beta0," H0 probability phi(beta0)=",phi(beta0)
set title strTitle(beta0)
unset title

plot[p=0:0.5]\
 P_H0(p, beta0, 1, 0.1) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=0.1" w l ls 1,\
 P_H0(p, beta0, 1, 0.5) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=0.5" w l ls 2,\
 P_H0(p, beta0, 1, 1) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=1" w l ls 3,\
 P_H0(p, beta0, 1, 2) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=2" w l ls 5,\
 P_H0(p, beta0, 1, 10) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=10" w l ls 7

##################################################################
set out "PH0_PriorGauss_beta0eq1sigbeta.png"
print "plotting PH0_PriorGauss_beta0eq1sigbeta.png"
##################################################################
beta0=1.
set key right bottom box
set label 1 at screen 0.65,0.45
set label 2 at screen 0.65,0.41
set label 3 at screen 0.80,0.51

#print "beta0=",beta0," H0 probability phi(beta0)=",phi(beta0)
#set title strTitle(beta0)
set label 3 sprintf("P(H_0)=%1.2f", phi(beta0))


plot[p=0:0.5]\
 P_H0(p, beta0, 1, 0.1) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=0.1" w l ls 1,\
 P_H0(p, beta0, 1, 0.5) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=0.5" w l ls 2,\
 P_H0(p, beta0, 1, 1) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=1" w l ls 3,\
 P_H0(p, beta0, 1, 2) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=2" w l ls 5,\
 P_H0(p, beta0, 1, 10) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=10" w l ls 7


##################################################################
set out "PH0_PriorGauss_beta0eq3sigbeta.png"
print "plotting PH0_PriorGauss_beta0eq3sigbeta.png"
##################################################################
beta0=3.
#set title strTitle(beta0)
set label 3 sprintf("P(H_0)=%1.4f", phi(beta0))

plot[p=0:0.5]\
 P_H0test(p, beta0, 1, 0.1) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=0.1" w l ls 1,\
 P_H0test(p, beta0, 1, 0.5) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=0.5" w l ls 2,\
 P_H0test(p, beta0, 1, 1) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=1" w l ls 3,\
 P_H0test(p, beta0, 1, 2) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=2" w l ls 5,\
 P_H0test(p, beta0, 1, 10) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=10" w l ls 7


##################################################################
set out "PH0_PriorGauss_beta0eqm1sigbeta.png"
print "plotting PH0_PriorGauss_beta0eqm1sigbeta.png"
##################################################################
beta0=-1.
set key right top box
set label 1 at screen 0.65,0.68
set label 2 at screen 0.65,0.64
set label 3 at screen 0.50,0.90

#set title strTitle(beta0)
set label 3 sprintf("P(H_0)=%1.2f", phi(beta0))

plot[p=0:0.5]\
 P_H0test(p, beta0, 1, 0.1) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=0.1" w l ls 1,\
 P_H0test(p, beta0, 1, 0.5) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=0.5" w l ls 2,\
 P_H0test(p, beta0, 1, 1) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=1" w l ls 3,\
 P_H0test(p, beta0, 1, 2) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=2" w l ls 5,\
 P_H0test(p, beta0, 1, 10) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=10" w l ls 7

##################################################################
set out "PH0_PriorGauss_beta0eqm0.5sigbeta.png"
print "plotting PH0_PriorGauss_beta0eqm0.5sigbeta.png"
##################################################################
beta0=-0.5
#set title strTitle(beta0)
set label 3 sprintf("P(H_0)=%1.2f", phi(beta0))

plot[p=0:0.5]\
 P_H0test(p, beta0, 1, 0.1) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=0.1" w l ls 1,\
 P_H0test(p, beta0, 1, 0.5) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=0.5" w l ls 2,\
 P_H0test(p, beta0, 1, 1) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=1" w l ls 3,\
 P_H0test(p, beta0, 1, 2) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=2" w l ls 5,\
 P_H0test(p, beta0, 1, 10) t "{/Symbol s}_b/{/Symbol s}_{/Symbol b}=10" w l ls 7

quit

print "\nIn-Text-Aufgabe MIV-Anteil:"

sigbeta=3.
sigb=3.
b_aufgabe=49.
beta0=55.
PH0=0.5
p=phi( (b_aufgabe-beta0)/sigb)

beta0=sigb*phiQuantil(PH0)

print "PH0=",PH0
print "p=",p
print "sigma=",sigma(sigbeta,sigb)
print "mu=",mu(p,beta0,sigbeta,sigb)
print "b=",b(p,beta0,sigb)
print "P_H0post_test=",P_H0test(p,beta0,sigbeta,sigb)
print "P_H0post=",P_H0(p,beta0,sigbeta,sigb)
print "p=",p


sigb=0.01; print "nun sigbeta=",sigbeta," sigb=",sigb,":"
print "P_H0post_test=",P_H0test(p,beta0,sigbeta,sigb)
print "P_H0post=",P_H0(p,beta0,sigbeta,sigb)
print "p=",p
