
# gnuplot4-Syntax!

# von  ~/info/gnuTemplate.gnu
# siehe auch
# ~/info/gnuplot,   
# ~/info/gnuTemplate42.gnu,  
# ~/info/gnuColoredContour/*.gnu
#http://www.chemie.fu-berlin.de/chemnet/use/info/gnuplot/gnuplot_27.html


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


set style line 1 lt 1 lw 10 pt 7 ps 1.9  lc rgb "#000000" #schwarz,solid,bullet
set style line 11 lt 1 lw 2 pt 1 ps 1.5  lc rgb "#000000" #schwarz,solid,plus sign
set style line 21 lt 1 lw 40 pt 7 ps 1.5  lc rgb "#999999" #schwarz,solid,plus sign
set style line 2 lt 1 lw 5 pt 2 ps 1.5  lc rgb "#CC0022" #rot, dash, Kreuz
set style line 3 lt 8 lw 2 pt 3 ps 1.2 #blassrot, offener star
set style line 4 lt 6 lw 5 pt 4 ps 1.5 #gelb, offenes Quadrat
set style line 5 lt 1 lw 5 pt 5 ps 1.5  lc rgb "#00DD22"  #gruen,solid,ClosedBox
set style line 6 lt 5 lw 7 pt 7 ps 1.5  lc rgb "#00DDDD" #blasstuerkisblau, offener Kreis
set style line 7 lt 1 lw 10 pt 7 ps 1.5  lc rgb "#1100AA"  #blau,solid,Bullet
set style line 8 lt 4 lw 2 pt 8 ps 1.5 #lila, aufrechtes geschloss. Dreieck
set style line 9 lt 7 lw 2 pt 9 ps 1.5 #schwarz, aufrechtes geschl. Dreieck

set style line 12 lt 1 lw 3 pt 12 ps 1.5  lc rgb "#CC0022" #rot, offene Raute
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

chinorm(nu)   = 1./(2**(0.5*nu)*gamma(0.5*nu))
chi2(x,nu)    =  chinorm(nu) *exp(-0.5*x) * x**(0.5*nu-1)
chi2Cum(x,nu)=igamma(0.5*nu,0.5*x)

#Newton-Iteration fuer Quantile; Startwert: nu
cQuantil_it(q,p,nu)=q-(chi2Cum(q,nu)-p)/chi2(q,nu)
chi2Quantil(p,nu)=cQuantil_it(cQuantil_it(cQuantil_it(cQuantil_it(cQuantil_it(\
  cQuantil_it(cQuantil_it(cQuantil_it(cQuantil_it(cQuantil_it(\
  nu,p,nu), p,nu), p,nu), p,nu), p,nu), p,nu),\
  p,nu), p,nu), p,nu), p,nu)

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x

set term post eps enhanced color solid "Helvetica" 18
#set term post eps enhanced color dashed "Helvetica" 20

####################################################
####################################################
# Data PT-Beispiel
# x1=Preis
# x2=Geschw.
# y=Zahl PT-Trips/Year/Person

n=10.
J=2.
alpha=0.05

t_nm1_0975=tQuantil(1-0.5*alpha,n-1)
t_nm2_0975=tQuantil(1-0.5*alpha,n-2)
t_nm3_0975=tQuantil(1-0.5*alpha,n-3)

x11=1;		x21=30;		y1=280;
x12=1.1;	x22=22;		y2=230;
x13=1.5;	x23=25;		y3=200;
x14=1.6;	x24=21;		y4=140;
x15=2.0;	x25=29;		y5=190;
x16=2.2;	x26=37;		y6=200;
x17=2.4;	x27=35;		y7=200;
x18=2.5;	x28=32;		y8=180;
x19=3.0;	x29=41;		y9=210;
x10=3.2;	x20=28;		y0=110;

x1datamin=0.8
x1datamax=3.3
x2datamin=10
x2datamax=45
ydatamin=100
ydatamax=285

####################################################
####################################################

x1bar=(x11+x12+x13+x14+x15+x16+x17+x18+x19+x10)/n
x2bar=(x21+x22+x23+x24+x25+x26+x27+x28+x29+x20)/n
ybar=(y1+y2+y3+y4+y5+y6+y7+y8+y9+y0)/n

sum11=x11**2+x12**2+x13**2+x14**2+x15**2+x16**2+x17**2+x18**2+x19**2+x10**2
sum22=x21**2+x22**2+x23**2+x24**2+x25**2+x26**2+x27**2+x28**2+x29**2+x20**2
sumyy=y1**2+y2**2+y3**2+y4**2+y5**2+y6**2+y7**2+y8**2+y9**2+y0**2

sum12=x11*x21+x12*x22+x13*x23+x14*x24+x15*x25+x16*x26+x17*x27+x18*x28+x19*x29+x10*x20
sum1y=x11*y1+x12*y2+x13*y3+x14*y4+x15*y5+x16*y6+x17*y7+x18*y8+x19*y9+x10*y0
sum2y=x21*y1+x22*y2+x23*y3+x24*y4+x25*y5+x26*y6+x27*y7+x28*y8+x29*y9+x20*y0

s11=sum11/n-x1bar**2
s22=sum22/n-x2bar**2
syy=sumyy/n-ybar**2

s12=sum12/n-x1bar*x2bar
s1y=sum1y/n-x1bar*ybar
s2y=sum2y/n-x2bar*ybar

#########################################################
print ""
print "Mittelwerte, Variance- and Kovarianzmatrizen"
print "====================================="
#########################################################

print "x1bar=",x1bar
print "x2bar=",x2bar
print "ybar=",ybar

print "s11=",s11,"\t s12=",s12
print "s21=",s12,"\t s22=",s22

print "s1y=",s1y
print "s2y=",s2y
print "syy=",syy

#########################################################
print ""
print "Regressionskoeffizienten beta0, beta1, beta2, Elastizitaeten elast1,elast2"
print "Deskr. (=normales) Bestimmtheitsma3 B"
print "=================================================="
#########################################################

dets12=s11*s22-s12**2

beta1=(s22*s1y-s12*s2y)/dets12
beta2=(s11*s2y-s12*s1y)/dets12
beta0=ybar-beta1*x1bar-beta2*x2bar

elast1=x1bar/ybar*beta1
elast2=x2bar/ybar*beta2

r12=s12/sqrt(s11*s22)
s_hatyhaty=beta1*s1y+beta2*s2y 
s_epseps=syy-s_hatyhaty  
B=s_hatyhaty/syy  #deskriptives (normales) Best.-ma3
Bquer(B,n,J)=1-(n-1)/(n-J-1) * (1-B)

print "detS=",dets12
print "beta1=",beta1," elast1=",elast1
print "beta2=",beta2," elast2=",elast2
print "beta0=",beta0

print "Korrelationen and deskr. Bestimmtheitsma3: r12=",r12," B=r^2(y, haty)=",B
print "deskriptive erklaerte Variance s_hatyhaty=",s_hatyhaty
print "deskr. nichterkl. Variance (Var.-Additionsregel) s_epseps=",s_epseps
print "induktives Bestimmtheitsma3 Bquer(B,n,2)=",Bquer(B,n,2)

#########################################################
print ""
print "Vergleich: Modelle nur einer exog. Var x1 bzw. x2"
print "==================================================="
#########################################################

B_m1=s1y**2/(s11*syy)
B_m2=s2y**2/(s22*syy)

print "Modell 1: Preis x1 allein:"
print "beta1_m1=s1y/s11=",s1y/s11," B_m1=",B_m1," Bquer_m1=",Bquer(B_m1,n,1)

print "Geschw. x2 allein:"
print "beta2_m2=s2y/s22=",s2y/s22," B_m2=",B_m2," Bquer_m2=",Bquer(B_m2,n,1)



#########################################################
print ""
print "Geschaetzte Kovarianzmatrix of the beta_j  (2 exogene Var x1,x2)"
print "======================================"
#########################################################

sig2_est=n/(n-3) * s_epseps  #Residualvarianz-Schaetzer 
Vbeta11=1./n * sig2_est/(s11*(1-r12**2))
Vbeta22=1./n * sig2_est/(s22*(1-r12**2))
Vbeta12=- 1./n * sig2_est*r12**2/(s12*(1-r12**2)) # MINUS sign!!!
rbeta12=Vbeta12/sqrt(Vbeta11*Vbeta22)
Vbeta_00=sig2_est/n+ x1bar**2*Vbeta11+x2bar**2*Vbeta22\
 +2*x1bar*x2bar*Vbeta12

print "Residualvarianz-Schaetzer sig2_est=",sig2_est
print "Vbeta_00=",Vbeta_00," Vbeta_01=n.a.  Vbeta_02=n.a"
print  "Vbeta_10=n.a., Vbeta11=",Vbeta11," Vbeta12=",Vbeta12
print  "Vbeta_20=n.a., Vbeta_21=",Vbeta12," Vbeta22=",Vbeta22
print "rbeta12=",rbeta12

#########################################################
print ""
print "Confidence Intervals for alpha=0.05"
print "======================================"
#########################################################


dbeta1=t_nm3_0975*sqrt(Vbeta11)
dbeta2=t_nm3_0975*sqrt(Vbeta22)

print "Quantil t^(",n-3,")_{0.975}=",t_nm3_0975
print "Halbe Breiten of the KI and KI selbst:"
print "dbeta1=",dbeta1,"   KI(beta_1)=[",beta1-dbeta1,",",beta1+dbeta1,"]"
print "dbeta2=",dbeta2,"   KI(beta_2)=[",beta2-dbeta2,",",beta2+dbeta2,"]"


###### m1=Modell nur x1, m2=Modell nur x2 #########

sig2_est_m1=n/(n-2) * s_epseps  #Residualvarianz-Schaetzer 
sig2_est_m2=sig2_est_m1
Vbeta11_m1=1./n * sig2_est_m1/s11
Vbeta22_m2=1./n * sig2_est_m2/s22
dbeta1_m1=t_nm2_0975*sqrt(Vbeta11_m1)
dbeta2_m2=t_nm2_0975*sqrt(Vbeta22_m2)


print "\nHalbe Breiten of the KI and KI selbst bei Modellen nur x1 bzw. x2:"
print "Quantil t^(",n-2,")_{0.975}=",t_nm2_0975
print "sqrt(Vbeta11_m1)=",sqrt(Vbeta11_m1)
print "dbeta1_m1=",dbeta1_m1,"   KI(beta1_m1)=[",s1y/s11-dbeta1_m1,",",s1y/s11+dbeta1_m1,"]"
print "dbeta2_m2=",dbeta2_m2,"   KI(beta2_m2)=[",s2y/s22-dbeta2_m2,",",s2y/s22+dbeta2_m2,"]"

#######################################
print "plotting OEV_f_gaussStudent_eng.eps"
set out "OEV_f_gaussStudent_eng.eps"
#######################################
set noparam
xnormmin=-2.
xnormmax=6.5
set xlabel "z bzw. t"
set xrange [xnormmin:xnormmax]
set ylabel "Density f_z(z) bzw. f_t(t)"
set auto y
set key

plot\
 gauss(0,1,x) t "Standard Normal Distribution" w l ls 1,\
 student(x,1) t "Student-Distribution {/Symbol n}=1 DOF" w l ls 2,\
 student(x,2) t "Student-Distribution {/Symbol n}=2 DOF" w l ls 3


#######################################
print "plotting OEV_F_gaussStudent_eng.eps"
set out "OEV_F_gaussStudent_eng.eps"
#######################################

set ylabel "Distribution Function F_z(z) bzw. F_t(t)"
set auto y
set key center right
set param

plot[x=xnormmin:xnormmax]\
 x,phi(x) t "Standard Normal Distribution" w l ls 1,\
 x,studentCum(x,1) t "Student-Distribution {/Symbol n}=1 DOF" w l ls 2,\
 x,studentCum(x,2) t "Student-Distribution {/Symbol n}=2 DOF" w l ls 3,\
 x,0.95 t "" w l ls 99,\
 phiQuantil(0.95), (x-xnormmin)/(xnormmax-xnormmin) t "z_{0.95}" w l ls 99,\
 6.314, (x-xnormmin)/(xnormmax-xnormmin) t "t^{(1)}_{0.95}" w l ls 2

# 0, (x-xnormmin)/(xnormmax-xnormmin) t "" w l ls 99


#######################################
print "plotting OEV_f_hatbeta1_eng.eps"
set out "OEV_f_hatbeta1_eng.eps"
#######################################

set title "Confidence Interval for H_0: {/Symbol b}_1 and {/Symbol s}  as estimated and {/Symbol a}=0.05"
tmin=-1.5*t_nm3_0975
tmax=2.1*t_nm3_0975
sigbeta1=sqrt(Vbeta11)
sigbeta2=sqrt(Vbeta22)
beta1min=beta1+tmin*sigbeta1
beta1max=beta1+tmax*sigbeta1
beta2min=beta2+tmin*sigbeta2
beta2max=beta2+tmax*sigbeta2

studentNoNorm(x,mux,hatsigx,nu)=1./hatsigx*student( (x-mux)/hatsigx,nu)

set xlabel "hat {/Symbol b}_1"
set xrange [beta1min:beta1max]
set ylabel "Density f(hat {/Symbol b}_1)"
set param
set key at screen 0.94,0.78



plot [t=tmin:tmax]\
  beta1+t*sigbeta1,studentNoNorm(beta1+t*sigbeta1,beta1,sigbeta1,n-3)\
    t "Density" w l ls 2,\
  beta1-dbeta1 + (t-tmin)/(tmax-tmin)*2*dbeta1,0 t "Confidence Interval" w l ls 21


#######################################
print "plotting OEV_F_hatbeta1_eng.eps"
set out "OEV_F_hatbeta1_eng.eps"
#######################################

set key right center
set ylabel "Distribution Function F(hat {/Symbol b}_1)"

plot [t=tmin:tmax]\
  beta1+t*sigbeta1,studentCum(t,n-3) t "F" w l ls 2,\
  beta1-dbeta1 + (t-tmin)/(tmax-tmin)*2*dbeta1,0 t "KI" w l ls 21,\
  beta1-dbeta1, (t-tmin)/(tmax-tmin) t "" w l ls 11,\
  beta1+dbeta1, (t-tmin)/(tmax-tmin) t "" w l ls 11,\
  beta1+t*sigbeta1,0.025 t "F=0.025" w l ls 16,\
  beta1+t*sigbeta1,0.975 t "F=0.975" w l ls 17

#######################################
print "plotting OEV_f_hatbeta2_eng.eps"
set out "OEV_f_hatbeta2_eng.eps"
#######################################

set title "Confidence Interval for H_0: {/Symbol b}_2 and {/Symbol s} as estimated and {/Symbol a}=0.05"

set key right top
set xlabel "hat {/Symbol b}_2"
set xrange [beta2min:beta2max]
set ylabel "Density f(hat {/Symbol b}_2)"
plot [t=tmin:tmax]\
  beta2+t*sigbeta2,studentNoNorm(beta2+t*sigbeta2,beta2,sigbeta2,n-3)\
    t "Density" w l ls 2,\
  beta2-dbeta2 + (t-tmin)/(tmax-tmin)*2*dbeta2,0 t "Confidence Interval" w l ls 21

#######################################
print "plotting OEV_F_hatbeta2_eng.eps"
set out "OEV_F_hatbeta2_eng.eps"
#######################################

set key right center
set ylabel "Distribution Function F(hat {/Symbol b}_1)"

plot [t=tmin:tmax]\
  beta2+t*sigbeta2,studentCum(t,n-3) t "F" w l ls 2,\
  beta2-dbeta2 + (t-tmin)/(tmax-tmin)*2*dbeta2,0 t "KI" w l ls 21,\
  beta2-dbeta2, (t-tmin)/(tmax-tmin) t "" w l ls 11,\
  beta2+dbeta2, (t-tmin)/(tmax-tmin) t "" w l ls 11,\
  beta2+t*sigbeta2,0.025 t "F=0.025" w l ls 16,\
  beta2+t*sigbeta2,0.975 t "F=0.975" w l ls 17


#######################################
print "plotting OEV_f_studentTestvar_eng.eps"
set out "OEV_f_studentTestvar_eng.eps"
#######################################

set title "Confidence Interval for H_0: {/Symbol b}_{1/2} and {/Symbol s} as estimated unter  {/Symbol a}=0.05"

set key right top
set xlabel "Test Variable t"
set xrange [tmin:tmax]
set ylabel "Density f(t)"
plot [t=tmin:tmax]\
  t,student(t,n-3) t "Density" w l ls 2,\
  -t_nm3_0975 + (t-tmin)/(tmax-tmin)*2*t_nm3_0975,0 t "Confidence Interval" w l ls 21

#######################################
print "plotting OEV_f_student_KI_eng.eps"
set out "OEV_f_student_KI_eng.eps"
#######################################
tmin=-4
tmax=4
set xrange [tmin:tmax]

plot [t=tmin:tmax]\
  t,student(t,2) t "Density" w l ls 2,\
  -tQuantil(0.95,2) + (t-tmin)/(tmax-tmin)*2*tQuantil(0.95,2),0 t "Confidence Interval" w l ls 21

tmin=-1.5*t_nm3_0975
tmax=2.1*t_nm3_0975
set xrange [tmin:tmax]

#######################################
print "plotting OEV_F_studentTestvar_eng.eps"
set out "OEV_F_studentTestvar_eng.eps"
#######################################

set ylabel "Distribution Function F(t)"

plot [t=tmin:tmax]\
  t,studentCum(t,n-3) t "F" w l ls 2,\
  -t_nm3_0975 + (t-tmin)/(tmax-tmin)*2*t_nm3_0975,0 t "KI" w l ls 21,\
  -t_nm3_0975, (t-tmin)/(tmax-tmin) t "" w l ls 11,\
  t_nm3_0975, (t-tmin)/(tmax-tmin) t "" w l ls 11,\
  -t_nm3_0975,0.025 t "F=0.025" w l ls 16,\
  t_nm3_0975,0.975 t "F=0.975" w l ls 17

#######################################
print "plotting OEV_F_student_KI_eng.eps"
set out "OEV_F_student_KI_eng.eps"
#######################################
tmin=-4
tmax=4
set xrange [tmin:tmax]
set key right center
t2_095=tQuantil(0.95,2)
plot [t=tmin:tmax]\
  t,studentCum(t,2) t "F" w l ls 2,\
  -t2_095 + (t-tmin)/(tmax-tmin)*2*t2_095,0 t "KI" w l ls 21,\
  -t2_095, (t-tmin)/(tmax-tmin) t "" w l ls 11,\
  t2_095, (t-tmin)/(tmax-tmin) t "" w l ls 11,\
  -t2_095,0.05 t "F=0.05" w l ls 16,\
  t2_095,0.95 t "F=0.95" w l ls 17

tmin=-1.5*t_nm3_0975
tmax=2.1*t_nm3_0975
set xrange [tmin:tmax]

#######################################
print "plotting OEV_scatterplot_x1x2_eng.eps"
set out "OEV_scatterplot_x1x2_eng.eps"
#######################################

set xlabel "Price x_1 (Euro)"
set ylabel "Speed x_2 (km/h)"
set nokey
set xrange [x1datamin:x1datamax]
set notitle

plot\
  x11,x21 w p ls 1,\
  x12,x22 w p ls 1,\
  x13,x23 w p ls 1,\
  x14,x24 w p ls 1,\
  x15,x25 w p ls 1,\
  x16,x26 w p ls 1,\
  x17,x27 w p ls 1,\
  x18,x28 w p ls 1,\
  x19,x29 w p ls 1,\
  x10,x20 w p ls 1

#######################################
print "plotting OEV_scatterplot_x1y_eng.eps"
set out "OEV_scatterplot_x1y_eng.eps"
#######################################

set xlabel "Price x_1 (Euro)"
set ylabel "Demand y (Trips/Person/Year)"
set key
set yrange [ydatamin:ydatamax]
set notitle
haty(x1,x2)=beta0+beta1*x1+beta2*x2
hatymod(x1,x2)=min(ydatamax,max(ydatamin,beta0+beta1*x1+beta2*x2))
hatymod2(x1,x2)=max(ydatamin,beta0+beta1*x1+beta2*x2)

plot[x1=x1datamin:x1datamax]\
  x11,y1 t "Data" w p ls 1,\
  x12,y2 t "" w p ls 1,\
  x13,y3 t "" w p ls 1,\
  x14,y4 t "" w p ls 1,\
  x15,y5 t "" w p ls 1,\
  x16,y6 t "" w p ls 1,\
  x17,y7 t "" w p ls 1,\
  x18,y8 t "" w p ls 1,\
  x19,y9 t "" w p ls 1,\
  x10,y0 t "" w p ls 1,\
  x1, ybar+s1y/s11*(x1-x1bar) t "Simple Regression(x_1)" w l ls 2,\
  x11,haty(x11,x21) t "Multiple Regression(x_1,x_2)" w p ls 3,\
  x12,haty(x12,x22) t "" w p ls 3,\
  x13,haty(x13,x23) t "" w p ls 3,\
  x14,haty(x14,x24) t "" w p ls 3,\
  x15,haty(x15,x25) t "" w p ls 3,\
  x16,haty(x16,x26) t "" w p ls 3,\
  x17,haty(x17,x27) t "" w p ls 3,\
  x18,haty(x18,x28) t "" w p ls 3,\
  x19,haty(x19,x29) t "" w p ls 3,\
  x10,haty(x10,x20) t "" w p ls 3

#######################################
print "plotting OEV_scatterplot_x2y_eng.eps"
set out "OEV_scatterplot_x2y_eng.eps"
#######################################

set xlabel "Speed x_2 (km/h)"
set key
set xrange [x2datamin:x2datamax]
set yrange [ydatamin:ydatamax]
set notitle


plot[x2=x2datamin:x2datamax]\
  x21,y1 t "Data" w p ls 1,\
  x22,y2 t "" w p ls 1,\
  x23,y3 t "" w p ls 1,\
  x24,y4 t "" w p ls 1,\
  x25,y5 t "" w p ls 1,\
  x26,y6 t "" w p ls 1,\
  x27,y7 t "" w p ls 1,\
  x28,y8 t "" w p ls 1,\
  x29,y9 t "" w p ls 1,\
  x20,y0 t "" w p ls 1,\
  x2, ybar+s2y/s22*(x2-x2bar) t "Simple Regression(x_1)" w l ls 2,\
  x21,haty(x11,x21) t "Multiple Regression(x_1,x_2)" w p ls 3,\
  x22,haty(x12,x22) t "" w p ls 3,\
  x23,haty(x13,x23) t "" w p ls 3,\
  x24,haty(x14,x24) t "" w p ls 3,\
  x25,haty(x15,x25) t "" w p ls 3,\
  x26,haty(x16,x26) t "" w p ls 3,\
  x27,haty(x17,x27) t "" w p ls 3,\
  x28,haty(x18,x28) t "" w p ls 3,\
  x29,haty(x19,x29) t "" w p ls 3,\
  x20,haty(x10,x20) t "" w p ls 3


#######################################
print "plotting OEV_f2_hatbeta1_hatbeta2_eng.eps"
set out "OEV_f2_hatbeta1_hatbeta2_eng.eps"
#######################################

set param
set view map
set multiplot
set contour surface
set cntrparam bspline
set cntrparam levels 10 
set isosamples 30,30
set palette defined ( 0 "white", 5 "yellow", 30 "orange",\
  70 "red",  100 "#99003") 

xmin=beta1-4*sigbeta1
#xmin=0
#xmax=beta1+4*sigbeta1
xmax=0

#ymin=beta2-4*sigbeta2
ymin=0
ymax=beta2+4*sigbeta2
#ymax=0

set xlabel "hat {/Symbol b}_1"
set xrange [xmin:xmax]

set ylabel "hat {/Symbol b}_2"
set yrange [ymin:ymax]

set pm3d
unset clabel
set title "2d-Density hat {/Symbol b}_j unter H_0: {/Symbol s} and {/Symbol b}_j as gemessen"
print "rbeta12=",rbeta12

set nokey
splot[x=xmin:xmax][y=ymin:ymax]\
  x,y,gauss2d(x-beta1,y-beta2,sigbeta1,sigbeta2,rbeta12)\
   t "f_2(hat {/Symbol b}_1, hat {/Symbol b}_2)" w l ls 98

set cntrparam levels discrete 0.0005  # eine Contour !
set key at screen 0.7,0.31
unset pm3d; unset surface
splot[x=xmin:xmax][y=ymin:ymax]\
  x,y,gauss2d(x-beta1,y-beta2,sigbeta1,sigbeta2,rbeta12)\
   t "H_0 - Confidence Region F-Test" w l ls 7

set surface
set key at screen 0.7,0.26
splot[t=0:1]\
 beta1-dbeta1,beta2-dbeta2+2*t*dbeta2,100+0.01*t\
      t "H_0 - Confidence Region t-Test" w l ls 5,\
 beta1+dbeta1,beta2-dbeta2+2*t*dbeta2,100+0.01*t t "" w l ls 5,\
 beta1-dbeta1+2*t*dbeta1,beta2-dbeta2,100+0.01*t t "" w l ls 5,\
 beta1-dbeta1+2*t*dbeta1,beta2+dbeta2,100+0.01*t t "" w l ls 5
set nomultiplot
set notitle

#######################################
print "plotting OEV_hyperbel_M1_eng.eps"
set out "OEV_hyperbel_M1_eng.eps"
#######################################

haty_m1(x1)=ybar+s1y/s11*(x1-x1bar)

haty_m1_sig2(x1)=sig2_est_m1/n*(1.+(x1-x1bar)**2/s11)
haty_m1_sig(x1)=sqrt(haty_m1_sig2(x1))
#fdens_haty_m1(y,x1)=gauss(haty_m1(x1),haty_m1_sig(x1),y)
fdens_haty_m1(y,x1)=1./haty_m1_sig(x1)*student( (y-haty_m1(x1))/haty_m1_sig(x1), n-2)
#Fdens_haty_m1(y,x1)=phi( (y-haty_m1(x1))/haty_m1_sig(x1))
Fdens_haty_m1(y,x1)=studentCum( (y-haty_m1(x1))/haty_m1_sig(x1), n-2)

print "t_nm2_0975=",t_nm2_0975," studentCum(t_nm2_0975,n-2)=", studentCum(t_nm2_0975,n-2)
print "studentCum(4,n-2)=", studentCum(4,n-2)
print "haty_m1_sig(100)/100.=",haty_m1_sig(100)/100.
print "sqrt(1./n * sig2_est_m1/s11)=",sqrt(1./n * sig2_est_m1/s11)


set multiplot

x1min=0.8
x1max=x1datamax+0.8*(x1datamax-x1datamin)


set xlabel "x_1 (Euro/Fahrt)"
set xrange [x1min:x1max]
set ylabel "y (Trips/Year/Person)"
set yrange [ydatamin:ydatamax]
unset cntrparam; set cntrparam levels 12 
set isosamples 30,40



unset surface
unset pm3d; set pm3d; set pm3d map; 
set surface; set contour surface
set key at screen 0.4,0.97
splot [x1=x1min:x1max][y=ydatamin:ydatamax]\
  x1,y,fdens_haty_m1(y,x1) t "Density f(hat y|x_1)" w l ls 98


unset pm3d; unset surface
unset cntrparam; set cntrparam levels discrete 0.025,0.975
set key at screen 0.4,0.92
splot [x1=x1min:x1max][y=ydatamin:ydatamax]\
  x1,y,Fdens_haty_m1(y,x1) t "2.5% and 97.5%-Quantile" w l ls 7

set surface
set key at screen 0.8,0.97
splot\
  x11,y1,101 t "Data" w p ls 1,\
  x12,y2,100 t "" w p ls 1,\
  x13,y3,100 t "" w p ls 1,\
  x14,y4,100 t "" w p ls 1,\
  x15,y5,100 t "" w p ls 1,\
  x16,y6,100 t "" w p ls 1,\
  x17,y7,100 t "" w p ls 1,\
  x18,y8,100 t "" w p ls 1,\
  x19,y9,100 t "" w p ls 1,\
  x10,y0,100 t "" w p ls 1

set key at screen 0.8,0.92
splot[t=0:1]\
 x1min+t*(x1max-x1min), haty_m1(x1min+t*(x1max-x1min)),100+t\
  t "hat y(x_1|Modell 1)" w l ls 2,\
 x1min+t*(x1max-x1min), ybar,100 t "x1quer,yquer" w l ls 99,\
 x1bar, ydatamin+t*(ydatamax-ydatamin), 100 t "" w l ls 99

set nomultiplot


#######################################
print "plotting OEV_scatter3d_eng.eps"
set out "OEV_scatter3d_eng.eps"
#######################################
set nogrid
unset colorbox

set param
set view 30,30
set multiplot
set contour surface
set cntrparam bspline
unset cntrparam; set cntrparam levels 20 
unset clabel
set isosamples 15,15
set palette defined ( 0 "white", 5 "yellow", 50 "orange",\
  80 "#FF6666",  99 "#FF44BB", 100 "white") 


set xlabel "x_1"
set xrange [x1datamin:x1datamax]
set ylabel "x_2"
set yrange [x2datamin:x2datamax]
set zlabel "y" offset 6,4
set zrange [ydatamin:ydatamax]

unset surface
set pm3d  
set contour surface

set pm3d hidden3d 98


set key at screen 0.7,0.91
splot[x1=x1datamin:x1datamax][x2=x2datamin:x2datamax]\
  x1,x2,hatymod(x1,x2)\
   t "{/Symbol b}_0+{/Symbol b}_1 x_1+{/Symbol b}_2 x_2" w l ls 98

set key at screen 0.7,0.86
unset pm3d; 
set surface
splot\
 x11,x21,y1 t "Data" w p ls 1,\
 x12,x22,y2 t "" w p ls 1,\
 x13,x23,y3 t "" w p ls 1,\
 x14,x24,y4 t "" w p ls 1,\
 x15,x25,y5 t "" w p ls 1,\
 x16,x26,y6 t "" w p ls 1,\
 x17,x27,y7 t "" w p ls 1,\
 x18,x28,y8 t "" w p ls 1,\
 x19,x29,y9 t "" w p ls 1,\
 x10,x20,y0 t "" w p ls 1



unset pm3d; 
set surface
set key at screen 0.7,0.81
splot[t=0:1]\
 x11,x21,haty(x11,x21)+t*(y1-haty(x11,x21)) t "{/Symbol e}_i" w l ls 6,\
 x12,x22,haty(x12,x22)+t*(y2-haty(x12,x22)) t "" w l ls 6,\
 x13,x23,haty(x13,x23)+t*(y3-haty(x13,x23)) t "" w l ls 6,\
 x14,x24,haty(x14,x24)+t*(y4-haty(x14,x24)) t "" w l ls 6,\
 x15,x25,haty(x15,x25)+t*(y5-haty(x15,x25)) t "" w l ls 6,\
 x16,x26,haty(x16,x26)+t*(y6-haty(x16,x26)) t "" w l ls 6,\
 x17,x27,haty(x17,x27)+t*(y7-haty(x17,x27)) t "" w l ls 6,\
 x18,x28,haty(x18,x28)+t*(y8-haty(x18,x28)) t "" w l ls 6,\
 x19,x29,haty(x19,x29)+t*(y9-haty(x19,x29)) t "" w l ls 6,\
 x10,x20,haty(x10,x20)+t*(y0-haty(x10,x20)) t "" w l ls 6

set nomultiplot

