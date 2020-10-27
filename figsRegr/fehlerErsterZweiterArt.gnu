
# gnuplot4-Syntax!

# von  ~/info/gnuTemplate.gnu
# siehe auch
# ~/info/gnuplot,   
# ~/info/gnuTemplate42.gnu,  
# ~/info/gnuColoredContour/*.gnu



set encoding iso_8859_1   # dann äöüßÄÖÜ durch woertl. Eingabe korrekt
#\304     "A
#\326     "O
#\334     "U
#\344     "a
#\366     "o
#\374     "u
#\337     &szlig;

##########################################################
# (Line- und Punkttypen Oct. 2009)
#geordnet nach hue (ps)
# BUG: Manchmal
#     obskure Abhaengigkeit der Symbole von Reihenfolge bzw linestyles
##########################################################

# set style line: linetype lt, point type pt

# if independently color, dash-style, point-style gnuplot 4.2 and later
set style line 99 lt 1 lw 3 pt 4 ps 1.5 lc rgb "#1100EE" #blau, solid, open box

set style line 1 lt 1 lw 2 pt 1 ps 1.5  lc rgb "#000000" #schwarz,solid,plus sign
set style line 2 lt 7 lw 8 pt 2 ps 1.5  lc rgb "#CC0022" #rot, dash, Kreuz
set style line 3 lt 8 lw 5 pt 3 ps 1.5 #blassrot, offener Stern
set style line 4 lt 6 lw 5 pt 4 ps 1.5 #gelb, offenes Quadrat
set style line 5 lt 1 lw 12 pt 5 ps 1.5  lc rgb "#00DD22"  #gruen,solid,ClosedBox
set style line 6 lt 5 lw 5 pt 6 ps 1.5 #blasstuerkisblau, offener Kreis
set style line 7 lt 1 lw 2 pt 7 ps 1.5  lc rgb "#1100AA"  #blau,solid,Bullet
set style line 8 lt 4 lw 5 pt 8 ps 1.5 #lila, aufrechtes geschloss. Dreieck
set style line 9 lt 7 lw 5 pt 9 ps 1.5 #schwarz, aufrechtes geschl. Dreieck
set style line 10 lt 1 lw 5 pt 10 ps 1.5 #rot, upside-down offenes Dreieck
set style line 11 lt 7 lw 10 pt 11 ps 1.5 #schwarz, upside-down geschl. Dreieck
set style line 12 lt 1 lw 5 pt 12 ps 1.5 #rot, offene Raute
set style line 13 lt 8 lw 5 pt 13 ps 1.5 #blassrot, geschl. Raute
set style line 14 lt 6 lw 5 pt 14 ps 1.5 #gelb, "Sonderzeichen" ...
set style line 15 lt 2 lw 5 pt 15 ps 1.5 #gruen, geschl. Fuenfeck
set style line 16 lt 5 lw 5 pt 16 ps 1.5 #blasstuerkisblau, 
set style line 17 lt 3 lw 5 pt 17 ps 1.5 #blau, 
set style line 18 lt 4 lw 5 pt 18 ps 1.5 #lila,
set style line 19 lt 7 lw 5 pt 19 ps 1.5
set style line 20 lt 1 lw 5 pt 20 ps 1.5



############### Beispiele fuer Funktionen ####################

gauss(mu,sigma,x) = exp(-(x-mu)**2/(2* sigma**2)) / (sigma * sqrt(2*pi))

xlimited(x)       =(x<-5) ? -5 : ( (x>5) ? 5 : x)

  # integral der Standardnormalverteilung
phi(x)            =0.5*(erf(xlimited(x/sqrt(2)))+1) 
phiTest(x)=norm(x)
                
 # Quantil der Standardnormalverteilung
phiQuantil(q)=sqrt(2.)*inverf(-1.+2*q)
phiQuantilTest(q)=invnorm(q)

#print phi(0.95)
#print phiTest(0.95)


lorenz(x)         = 1/(pi*(1+x**2))
expo(lambda,x)    = lambda*exp(-lambda*x)
gleich(a,b,x)     = ((x>=a)&&(x<=b)) ? 1/(b-a) : 0

#gamma(n+1)=n! for integer n>=0
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

############################
 # Chi-Quadrat-Verteilung und Derivate, nu=Zahl der FG
 ############################

chinorm(nu)   = 1./(2**(0.5*nu)*gamma(0.5*nu))
chi2(x,nu)    =  chinorm(nu) *exp(-0.5*x) * x**(0.5*nu-1)
chi2Cum(x,nu)=igamma(0.5*nu,0.5*x)

#Newton-Iteration fuer Quantile; Startwert: nu
cQuantil_it(q,p,nu)=q-(chi2Cum(q,nu)-p)/chi2(q,nu)
chi2Quantil(p,nu)=cQuantil_it(cQuantil_it(cQuantil_it(cQuantil_it(cQuantil_it(\
  cQuantil_it(cQuantil_it(cQuantil_it(cQuantil_it(cQuantil_it(\
  nu,p,nu), p,nu), p,nu), p,nu), p,nu), p,nu),\
  p,nu), p,nu), p,nu), p,nu)



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

t1_09=3.078
t2_09=1.886
t3_09=1.638
t4_09=1.533

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x

####################################################

#set term post eps enhanced color solid "Helvetica" 20
set term post eps enhanced color dashed "Helvetica" 20
set size 1,0.7
set param 
set key

# Guetefunktionen fuer parametr. Tests auf Mittelwert mit sigma2=bekannt

alpha=0.05  # Fehlerwahrcheinlichkeit
m=2       # Zahl der freiheitsgrade der Student-t-Verteilung
m_chi2=5   # Zahl der freiheitsgrade der chi^2-Verteilung

#############################################
# Guetefunktion fuer Punkttest H0: mu=mu0; z=(beta-beta0)/sigma_beta
# Guetefunktion fuer Punkttest H0: mu=mu0; t=(beta-beta0)/hat(sigma_beta)
#############################################

zHalb(alpha)=phiQuantil(0.5*alpha)
tHalb(alpha,m)=tQuantil(0.5*alpha,m)

#############################################

GfunEqGauss(z,alpha)=phi(phiQuantil(0.5*alpha)-z)\
    + (1.-phi(-phiQuantil(0.5*alpha)-z))

GfunEqStudent(t,alpha,m)=studentCum(tQuantil(0.5*alpha,m)-t,m)\
    + (1.-studentCum(-tQuantil(0.5*alpha,m)-t,m))

GfunLeGauss(z,alpha)=1-phi(phiQuantil(1-alpha)-z)

GfunLeStudent(t,alpha,m)=1-studentCum(tQuantil(1-alpha,m)-t,m)

GfunGeGauss(z,alpha)=phi(phiQuantil(alpha)-z)

GfunGeStudent(t,alpha,m)=studentCum(tQuantil(alpha,m)-t,m)

#############################################



set out "fehler12art_eqGauss.eps"
print "plotting fehler12art_eqGauss.eps"
set key center
set xlabel "{/Symbol D}Z={/Symbol (b-b_0)/s_b}"
zmin=-5.
zmax=5.
set xrange [zmin:zmax]

set ylabel "Wahrscheinlichkeit"
set yrange [0:1.01]
set ytics 0.2

plot[t=0:1]\
 0, GfunEqGauss(0,alpha) t "{/Symbol a}-Fehler" w p ls 5,\
 0.1+t*(zmax-0.1), 1-GfunEqGauss(0.1+t*(zmax-0.1),alpha)\
      t "{/Symbol b}-Fehler" w l ls 2,\
 -0.1+t*(zmin+0.1), 1-GfunEqGauss(-0.1+t*(zmin+0.1),alpha)\
      t "" w l ls 2,\
 zmin+t*(zmax-zmin), GfunEqGauss(zmin+t*(zmax-zmin),alpha)\
        t "G\374tefunktion" w l ls 7

#############################################
set out "fehler12art_eqStudent_2FG.eps"
print "plotting fehler12art_eqStudent_2FG.eps"
#############################################

tmin=-5.
tmax=5.
set xrange [tmin:tmax]
set xlabel "{/Symbol (b-b_0)/s_b}"  # wie bei bekannter Varianz; sind ja WAHRE Werte

plot[t=0:1]\
  0, GfunEqStudent(0,alpha,m) t "{/Symbol a}-Fehler" w p ls 5,\
  0.1+t*(tmax-0.1), 1-GfunEqStudent(0.1+t*(tmax-0.1),alpha,m)\
      t "{/Symbol b}-Fehler" w l ls 2,\
 -0.1+t*(tmin+0.1), 1-GfunEqStudent(-0.1+t*(tmin+0.1),alpha,m)\
      t "" w l ls 2,\
 tmin+t*(tmax-tmin), GfunEqStudent(tmin+t*(tmax-tmin),alpha,m)\
        t "G\374tefunktion" w l ls 7


#############################################
# Guetefunktion fuer Intervalltest H0: mu<=mu0; z=(beta-beta0)/sigma_beta
# Guetefunktion fuer Intervalltest H0: mu<=mu0; t=(beta-beta0)/hat(sigma)_beta
#############################################

set out "fehler12art_leGauss.eps"
print "plotting fehler12art_leGauss.eps"
set key left top


set xlabel "{/Symbol (b-b_0)/s_b}"
set xrange [zmin:zmax]

plot[t=0:1]\
 zmin+t*(0.-zmin), GfunLeGauss(zmin+t*(0.-zmin),alpha)\
      t "{/Symbol a}-Fehler" w l ls 5,\
 0.02+t*(zmax-0.02), 1.-GfunLeGauss(0.02+t*(zmax-0.02),alpha)\
      t "{/Symbol b}-Fehler" w l ls 2,\
 zmin+t*(zmax-zmin), GfunLeGauss(zmin+t*(zmax-zmin),alpha)\
        t "G\374tefunktion" w l ls 7

#############################################
set out "fehler12art_leStudent_2FG.eps"
print "plotting fehler12art_leStudent_2FG.eps"
#############################################

set xrange [tmin:tmax]
set xlabel "{/Symbol (b-b_0)/s_b}"  # wie bei bekannter Varianz; sind ja WAHRE Werte

plot[t=0:1]\
 tmin+t*(0.-tmin), GfunLeStudent(tmin+t*(0.-tmin),alpha,m)\
      t "{/Symbol a}-Fehler" w l ls 5,\
 0.02+t*(tmax-0.02), 1.-GfunLeStudent(0.02+t*(tmax-0.02),alpha,m)\
      t "{/Symbol b}-Fehler" w l ls 2,\
 tmin+t*(tmax-tmin), GfunLeStudent(tmin+t*(tmax-tmin),alpha,m)\
        t "G\374tefunktion" w l ls 7

#############################################
# Guetefunktion fuer Intervalltest H0: mu>=mu0; y=(mu-mu0)/sigma
#############################################

set key right
set out "fehler12art_geGauss.eps"
print "plotting fehler12art_geGauss.eps"


set xrange [zmin:zmax]
set xlabel "{/Symbol D}Z={/Symbol (b-b_0)/s_b}"

plot[t=0:1]\
 0.02+t*(zmax-0.02), GfunGeGauss(0.02+t*(zmax-0.02),alpha)\
      t "{/Symbol a}-Fehler" w l ls 5,\
 zmin+t*(0.-zmin), 1.-GfunGeGauss(zmin+t*(0.-zmin),alpha)\
      t "{/Symbol b}-Fehler" w l ls 2,\
 zmin+t*(zmax-zmin), GfunGeGauss(zmin+t*(zmax-zmin),alpha)\
        t "G\374tefunktion" w l ls 7

#############################################
set out "fehler12art_geStudent_2FG.eps"
print "plotting fehler12art_geStudent_2FG.eps"
#############################################

set xrange [tmin:tmax]
set xlabel "{/Symbol (b-b_0)/s_b}"  # wie bei bekannter Varianz; sind ja WAHRE Werte


plot[t=0:1]\
 0.02+t*(tmax-0.02), GfunGeStudent(0.02+t*(tmax-0.02),alpha,m)\
      t "{/Symbol a}-Fehler" w l ls 5,\
 tmin+t*(0.-tmin), 1.-GfunGeStudent(tmin+t*(0.-tmin),alpha,m)\
      t "{/Symbol b}-Fehler" w l ls 2,\
 tmin+t*(tmax-tmin), GfunGeStudent(tmin+t*(tmax-tmin),alpha,m)\
        t "G\374tefunktion" w l ls 7




#############################################
# Guetefunktion fuer Punkttest H0: sigma^2=sigma_0^2; y=sigma/sigma_0
#############################################

set out "fehler12art_eqChi2_5FG.eps"
print "plotting fehler12art_eqChi2_5FG.eps"
set key center center
set label 1 "5 Freiheitsgrade" at screen 0.7,0.21

#qHalb(alpha,m)=chi2Quantil(0.5*alpha,m) # 1.145        #q^(5)_{alpha/2} Qantil
#qAlpha(alpha,m)=chi2Quantil(alpha,m) #1.610       #q^(5)_{alpha} Qantil
#q1mAlpha(alpha,m)=chi2Quantil(1-alpha,m)  #9.236  #q^(5)_{1-alpha} Qantil
#q1mHalb(alpha,m)= chi2Quantil(1-0.5*alpha,m)# 11.07   #q^(5)_{1-alpha/2} Qantil

#############################################

GfunVarEqual(y,alpha,m)=chi2Cum(chi2Quantil(0.5*alpha,m)/y**2,m)\
  +(1.-chi2Cum(chi2Quantil(1-0.5*alpha,m)/y**2,m))

GfunVarLe(y,alpha,m)=1-chi2Cum(chi2Quantil(1-alpha,m)/y**2,m)

GfunVarGe(y,alpha,m)=chi2Cum(chi2Quantil(alpha,m)/y**2,m)

#############################################


set xlabel "{/Symbol s/s_0}"
ymin=0.2
ymax=2.
set xrange [ymin:ymax]

set ylabel "Wahrscheinlichkeit"

plot[t=0:1]\
 1, GfunVarEqual(1,alpha,m_chi2) t "{/Symbol a}-Fehler" w p ls 5,\
 1.04+t*(ymax-1.04), 1-GfunVarEqual(1.04+t*(ymax-1.04),alpha,m_chi2)\
      t "{/Symbol b}-Fehler" w l ls 2,\
 0.96+t*(ymin-0.96), 1-GfunVarEqual(0.96+t*(ymin-0.96),alpha,m_chi2)\
      t "" w l ls 2,\
 ymin+t*(ymax-ymin), GfunVarEqual(ymin+t*(ymax-ymin),alpha,m_chi2)\
        t "G\374tefunktion" w l ls 7


#############################################
# Guetefunktion fuer Intervalltest H0: sigma^2<=sigma0^2; y=sigma/sigma0
#############################################

set out "fehler12art_leChi2_5FG.eps"
print "plotting fehler12art_leChi2_5FG.eps"
set key left top



plot[t=0:1]\
 ymin+t*(1.-ymin), GfunVarLe(ymin+t*(1.-ymin),alpha,m_chi2)\
      t "{/Symbol a}-Fehler" w l ls 5,\
 1+t*(ymax-1), 1.-GfunVarLe(1+t*(ymax-1),alpha,m_chi2)\
      t "{/Symbol b}-Fehler" w l ls 2,\
 ymin+t*(ymax-ymin), GfunVarLe(ymin+t*(ymax-ymin),alpha,m_chi2)\
        t "G\374tefunktion" w l ls 7

#############################################
# Guetefunktion fuer Intervalltest H0: sigma^2>=sigma0^2; y=sigma/sigma0
#############################################

set out "fehler12art_geChi2_5FG.eps"
print "plotting fehler12art_geChi2_5FG.eps"

set key right top


plot[t=0:1]\
 ymin+t*(1.-ymin), 1-GfunVarGe(ymin+t*(1.-ymin),alpha,m_chi2)\
      t "{/Symbol b}-Fehler" w l ls 2,\
 1+t*(ymax-1), GfunVarGe(1+t*(ymax-1),alpha,m_chi2)\
      t "{/Symbol a}-Fehler" w l ls 5,\
 ymin+t*(ymax-ymin), GfunVarGe(ymin+t*(ymax-ymin),alpha,m_chi2)\
        t "G\374tefunktion" w l ls 7

