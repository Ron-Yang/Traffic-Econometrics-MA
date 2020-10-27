
###########################
# Template:
# von ~/vorlesungen/Verkehrsoekonometrie_Ma/skript/figsRegr/Standardverteilungen.gnu
############################

set encoding iso_8859_1   # dann äöüßÄÖÜ durch woertl. Eingabe korrekt
#\304     "A
#\326     "O
#\334     "U
#\344     "a
#\366     "o
#\374     "u
#\337     &szlig;

##########################################################
# (Line- and Punkttypen Oct. 2009)
#geordnet nach hue (ps)
# BUG: Manchmal
#     obskure Abhaengigkeit of the Symbole von Reihenfolge bzw linestyles
##########################################################

# set style line: linetype lt, point type pt

# if independently color, dash-style, point-style gnuplot 4.2 and later
set style line 98 lt 7 lw 0 linecolor rgb "#000000" # beliebige Farben:Schwarz
set style line 99 lt 1 lw 5 linecolor rgb "#000000" # beliebige Farben:Schwarz


set style line 1 lt 1 lw 1 pt 7 ps 1.9  lc rgb "#000000" #schwarz,solid,bullet
set style line 2 lt 1 lw 1 pt 2 ps 1.5  lc rgb "#CC0022" #rot, dash, Kreuz
set style line 3 lt 8 lw 1 pt 4 ps 1.2 #blassrot, offenes Quadrat
set style line 4 lt 6 lw 1 pt 4 ps 1.5  lc rgb "#FFAA00"  #gelb, offenes Quadrat
set style line 5 lt 1 lw 1 pt 5 ps 1.5  lc rgb "#00DD22"  #gruen,solid,ClosedBox
set style line 6 lt 5 lw 1 pt 4 ps 1.5  lc rgb "#00AAAA" #offenes Quadrat
set style line 7 lt 1 lw 1 pt 4 ps 2.0  lc rgb "#1100AA"  #blau,offenes Quadrat
set style line 8 lt 4 lw 1 pt 8 ps 1.5 #lila, aufrechtes geschloss. Dreieck
set style line 9 lt 7 lw 1 pt 9 ps 1.5  lc rgb "#999999"  #grau, aufr. gschl. Dreieck

set style line 11 lt 1 lw 6 pt 8 ps 1.9  lc rgb "#000000" #schwarz,dreieck
set style line 12 lt 1 lw 4 pt 2 ps 1.5  lc rgb "#CC0022" #rot, dash, Kreuz
set style line 13 lt 8 lw 4 pt 3 ps 1.2 #blassrot, offener star
set style line 14 lt 6 lw 4 pt 4 ps 1.5  lc rgb "#FFAA00"  #gelb, offenes Quadrat
set style line 15 lt 1 lw 4 pt 5 ps 1.5  lc rgb "#00DD22"  #gruen,solid,ClosedBox
set style line 16 lt 5 lw 4 pt 7 ps 1.5  lc rgb "#00AAAA" #offener Kreis
set style line 17 lt 1 lw 4 pt 7 ps 1.5  lc rgb "#1100AA"  #blau,solid,Bullet
set style line 18 lt 4 lw 4 pt 8 ps 1.5  lc rgb "#6600AA"  #lila, aufrechtes geschloss. Dreieck
set style line 19 lt 7 lw 4 pt 9 ps 1.5  lc rgb "#999999"  #grau, aufr. gschl. Dreieck


set style line 21 lt 1 lw 20 pt 7 ps 1.9  lc rgb "#000000" #schwarz,solid,bullet
set style line 22 lt 1 lw 20 pt 2 ps 1.5  lc rgb "#CC0022" #rot, dash, Kreuz
set style line 29 lt 7 lw 20 pt 9 ps 1.5  lc rgb "#999999"  #grau, aufr. gschl. Dreieck



############### Beispiele fuer Funktionen ####################

gauss(mu,sigma,x) = exp(-(x-mu)**2/(2* sigma**2)) / (sigma * sqrt(2*pi))
gauss2d(x,y,sx,sy,r) \
 = exp(-(0.5/(1-r*r))*((x/sx)**2+(y/sy)**2-2*r*x*y/(sx*sy)))\
  / (2*pi*sx*sy*sqrt(1-r*r))
gauss2d_invCovMatrix(x1,x2,invs11, invs12,invs22)\
 = sqrt(invs11*invs22-invs12**2) /  (2*pi)\
 * exp(-0.5*(x1*invs11*x1+2*x1*invs12*x2+x2*invs22*x2))


xlimited(x)       =(x<-5) ? -5 : ( (x>5) ? 5 : x)

  # integral of the Standard Normal Distribution
phi(x)            =norm(x)
phiOld(x)            =0.5*(erf(xlimited(x/sqrt(2)))+1)  #equivalent
                
 # Quantil of the Standard Normal Distribution
phiQuantil(q)=invnorm(q)
phiQuantilOld(q)=sqrt(2.)*inverf(-1.+2*q)  #equivalent
                


lorenz(x)         = 1/(pi*(1+x**2))
expo(lambda,x)    = lambda*exp(-lambda*x)
gleich(a,b,x)     = ((x>=a)&&(x<=b)) ? 1/(b-a) : 0

#gamma(n+1)=n! for integer n>=0
# nu=Zahl of the DOF
studnorm(nu)  = gamma(0.5*(nu+1)) / (sqrt(nu*pi)*gamma(0.5*nu))
student(x,nu) =  studnorm(nu) / (1.+x**2/nu)**(0.5*(nu+1))
studentCum(x,nu) = ibeta(0.5*nu,0.5*nu,0.5*(x+sqrt(x**2+nu))/sqrt(x**2+nu))

#Iteration; Startwert: phiQuantil(p)
tQuantil_it(t,p,nu)=t-(studentCum(t,nu)-p)/student(t,nu)
tQuantil(p,nu)=tQuantil_it(tQuantil_it(tQuantil_it(tQuantil_it(tQuantil_it(\
  tQuantil_it(tQuantil_it(tQuantil_it(tQuantil_it(tQuantil_it(\
  phiQuantil(p),p,nu), p,nu), p,nu), p,nu), p,nu), p,nu),\
  p,nu), p,nu), p,nu), p,nu)


betaFun(x,y)=gamma(x)*gamma(y)/(gamma(x+y))
fisher(x,d1,d2)=sqrt( (d1*x)**d1 * d2**d2/( (d1*x+d2)**(d1+d2)))\
 / (x*betaFun(0.5*d1, 0.5*d2))
fisherCum(x,d1,d2)=ibeta(0.5*d1, 0.5*d2, d1*x/(d1*x+d2))


chinorm(nu)   = 1./(2**(0.5*nu)*gamma(0.5*nu))
chi2(x,nu)    =  chinorm(nu) *exp(-0.5*x) * x**(0.5*nu-1)
chi2Cum(x,nu)=igamma(0.5*nu,0.5*x)

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x

t1_0975=12.70
t2_0975=4.303
t3_0975=3.182
t4_0975=2.776
t5_0975=2.571
t6_0975=2.447
t7_0975=2.365
t8_0975=2.306
t9_0975=2.262
t10_0975=2.228
t11_0975=2.201

t1_095=6.314
t2_095=2.920
t3_095=2.353
t4_095=2.132
t5_095=2.015
t6_095=1.943
t7_095=1.895
t8_095=1.860
t9_095=1.833
t10_095=1.812
t11_095=1.796

set term post eps enhanced color solid "Helvetica" 22



#######################################
print "plotting fisherCum_eng.eps"
set out "fisherCum_eng.eps"
#######################################

set param
set key bottom right
set xlabel "Argument f"
set ylabel "Cumulated Fisher-F-Distribution F_F^{(2,n)}(f)"
set yrange [0.6:]
xmax=14.
set xrange [0:xmax]
plot[t=0:1]\
  xmax*t, fisherCum(xmax*t,2,2) t "n=2" w l ls 12,\
  xmax*t, fisherCum(xmax*t,2,3) t "n=3" w l ls 13,\
  xmax*t, fisherCum(xmax*t,2,4) t "n=4" w l ls 14,\
  xmax*t, fisherCum(xmax*t,2,5) t "n=5" w l ls 15,\
  xmax*t, fisherCum(xmax*t,2,10) t "n=10" w l ls 16,\
  xmax*t, fisherCum(xmax*t,2,20) t "n=20" w l ls 17,\
  xmax*t, 0.95 t "F=0.95" w l ls 11

#######################################
print "plotting studentCum_eng.eps"
set out "studentCum_eng.eps"
#######################################

set xlabel "Argument t"
set ylabel "Cumulated Student-t-Distribution F_T^{(n)}(t)"
set yrange [0.8:]
xmax=5.
set xrange [0:xmax]
plot[t=0:1]\
  xmax*t, studentCum(xmax*t,1) t "n=1" w l ls 12,\
  xmax*t, studentCum(xmax*t,2) t "n=2" w l ls 13,\
  xmax*t, studentCum(xmax*t,3) t "n=3" w l ls 14,\
  xmax*t, studentCum(xmax*t,4) t "n=4" w l ls 15,\
  xmax*t, studentCum(xmax*t,5) t "n=5" w l ls 16,\
  xmax*t, studentCum(xmax*t,10) t "n=10" w l ls 17,\
  xmax*t, norm(xmax*t) t "n -> infinity (Gauss)" w l ls 18,\
  xmax*t, 0.95 t "F=0.95" w l ls 11

#######################################
print "plotting gaussCum_eng.eps"
set out "gaussCum_eng.eps"
#######################################

set xlabel "Argument z"
set ylabel "Cumulated Standard Normal Distribution {/Symbol F}(z)"
set yrange [0.8:]
xmax=3.
set xrange [0:xmax]

plot [t=0:1]\
  xmax*t, norm(xmax*t) t "" w l ls 18,\
  xmax*t, 0.95 t "F=0.95" w l ls 11,\
  xmax*t, 0.975 t "F=0.975" w l ls 1


#######################################
print "plotting chi2Cum_eng.eps"
set out "chi2Cum_eng.eps"
#######################################

set xlabel "Argument t"
set ylabel "Cumulated {/Symbol c}^2(n)-Distribution F"
set yrange [0.8:]
xmax=25.
set xrange [0:xmax]
plot[t=0:1]\
  xmax*t, chi2Cum(xmax*t,1) t "n=1" w l ls 12,\
  xmax*t, chi2Cum(xmax*t,2) t "n=2" w l ls 13,\
  xmax*t, chi2Cum(xmax*t,3) t "n=3" w l ls 14,\
  xmax*t, chi2Cum(xmax*t,4) t "n=4" w l ls 15,\
  xmax*t, chi2Cum(xmax*t,5) t "n=5" w l ls 16,\
  xmax*t, chi2Cum(xmax*t,6) t "n=6" w l ls 17,\
  xmax*t, chi2Cum(xmax*t,10) t "n=10" w l ls 18,\
  xmax*t, 0.95 t "F=0.95" w l ls 11


#######################################
print "plotting gaussQuantil_eng.eps"
set out "gaussQuantil_eng.eps"
#######################################

set key top left

xmin=0.8
xmax=1.
ymax=3.

set xlabel "Argument p"
set xrange [xmin:xmax]
set ylabel "Gauss-Quantile Function {/Symbol F}^{-1}(p)"
set yrange [0:ymax]

plot [t=xmin:xmax]\
  t, invnorm(t) t "" w l ls 18,\
  0.95, ymax/(xmax-xmin)*(t-xmin) t "p=0.95" w l ls 11,\
  0.975, ymax/(xmax-xmin)*(t-xmin) t "p=0.975" w l ls 1


#######################################
print "plotting studentQuantil_eng.eps"
set out "studentQuantil_eng.eps"
#######################################

set key top left

xmin=0.8
xmax=1.
ymax=5.

set xlabel "Argument p"
set xrange [xmin:xmax]
set ylabel "Student-Quantile Function {/Symbol F}^{-1}_{T(n)}(p)"
set yrange [0:ymax]

plot [t=xmin:xmax]\
  t, tQuantil(t,1) t "n=1" w l ls 12,\
  t, tQuantil(t,2) t "n=2" w l ls 13,\
  t, tQuantil(t,3) t "n=3" w l ls 14,\
  t, tQuantil(t,4) t "n=4" w l ls 15,\
  t, tQuantil(t,5) t "n=5" w l ls 16,\
  t, tQuantil(t,10) t "n=10" w l ls 17,\
  t, invnorm(t) t "Gauss" w l ls 18,\
  0.95, ymax/(xmax-xmin)*(t-xmin) t "p=0.95" w l ls 11,\
  0.975, ymax/(xmax-xmin)*(t-xmin) t "p=0.975" w l ls 1


#######################################
print "plotting gaussDensity_eng.eps"
set out "gaussDensity_eng.eps"
#######################################

set key top right
set noparam
xmin=-3
xmax=3
set xlabel "z"
set xrange [xmin:xmax]
set ylabel " Density f of the Standard Normal Distribution"
set auto y
plot[x=xmin:xmax] gauss(0,1,x) t "" w l ls 18 


#######################################
print "plotting studentDensity_eng.eps"
set out "studentDensity_eng.eps"
#######################################
xmin=-5
xmax=5
set xlabel "t"
set xrange [xmin:xmax]
set ylabel "Density f_{T(n)}(t)  of the Student-t Distributions"
set auto y
plot[x=xmin:xmax]\
 student(x,1) t "n=1" w l ls 12,\
 student(x,2) t "n=2" w l ls 13,\
 student(x,3) t "n=3" w l ls 14,\
 student(x,4) t "n=4" w l ls 15,\
 student(x,5) t "n=5" w l ls 16,\
 student(x,10) t "n=10" w l ls 17,\
 gauss(0,1,x) t "Gauss" w l ls 18

#######################################
print "plotting chi2Density_eng.eps"
set out "chi2Density_eng.eps"
#######################################
xmin=0
xmax=15
set xlabel "x"
set xrange [xmin:xmax]
set ylabel "Density f_{{/Symbol c}^2}^{(n)}(x) of the {/Symbol c}^2-Distributions"
set yrange [:0.7]
plot[x=xmin:xmax]\
 chi2(x,1) t "n=1" w l ls 12,\
 chi2(x,2) t "n=2" w l ls 13,\
 chi2(x,3) t "n=3" w l ls 14,\
 chi2(x,4) t "n=4" w l ls 15,\
 chi2(x,5) t "n=5" w l ls 16,\
 chi2(x,6) t "n=6" w l ls 17,\
 chi2(x,10) t "n=10" w l ls 18

