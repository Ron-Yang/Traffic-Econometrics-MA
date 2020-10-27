# (feb11) aus klausur_WS1011_Ma_regr.gnu
# Streamlined im Vgl zu Matrix1

###########################
# Template:
# von ~/vorlesungen/Verkehrsoekonometrie_Ma/skript/figsRegr/regr_2exogVarMatrix.gnu
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
# (Line- und Punkttypen Oct. 2009)
#geordnet nach hue (ps)
# BUG: Manchmal
#     obskure Abhaengigkeit der Symbole von Reihenfolge bzw linestyles
##########################################################

# set style line: linetype lt, point type pt

# if independently color, dash-style, point-style gnuplot 4.2 and later
set style line 98 lt 7 lw 0 linecolor rgb "#000000" # beliebige Farben:Schwarz
set style line 99 lt 1 lw 5 linecolor rgb "#000000" # beliebige Farben:Schwarz


set style line 1 lt 1 lw 2 pt 7 ps 1.9  lc rgb "#000000" #schwarz,solid,bullet
set style line 2 lt 1 lw 2 pt 2 ps 1.5  lc rgb "#CC0022" #rot, dash, Kreuz
set style line 3 lt 8 lw 2 pt 4 ps 1.2 #blassrot, offenes Quadrat
set style line 4 lt 6 lw 2 pt 4 ps 1.5  lc rgb "#FFAA00"  #gelb, offenes Quadrat
set style line 5 lt 1 lw 2 pt 5 ps 1.5  lc rgb "#00DD22"  #gruen,solid,ClosedBox
set style line 6 lt 5 lw 2 pt 4 ps 1.5  lc rgb "#00AAAA" #offenes Quadrat
set style line 7 lt 1 lw 2 pt 4 ps 2.0  lc rgb "#1100AA"  #blau,offenes Quadrat
set style line 8 lt 4 lw 2 pt 8 ps 1.5 #lila, aufrechtes geschloss. Dreieck
set style line 9 lt 7 lw 2 pt 9 ps 1.5  lc rgb "#999999"  #grau, aufr. gschl. Dreieck

set style line 11 lt 1 lw 6 pt 8 ps 1.9  lc rgb "#000000" #schwarz,dreieck
set style line 12 lt 1 lw 6 pt 2 ps 1.5  lc rgb "#CC0022" #rot, dash, Kreuz
set style line 13 lt 8 lw 6 pt 3 ps 1.2 #blassrot, offener Stern
set style line 14 lt 6 lw 6 pt 4 ps 1.5  lc rgb "#FFAA00"  #gelb, offenes Quadrat
set style line 15 lt 1 lw 6 pt 5 ps 1.5  lc rgb "#00DD22"  #gruen,solid,ClosedBox
set style line 16 lt 5 lw 6 pt 7 ps 1.5  lc rgb "#00AAAA" #offener Kreis
set style line 17 lt 1 lw 6 pt 7 ps 1.5  lc rgb "#1100AA"  #blau,solid,Bullet
set style line 18 lt 4 lw 6 pt 8 ps 1.5  lc rgb "#6600AA"  #lila, aufrechtes geschloss. Dreieck
set style line 19 lt 7 lw 6 pt 9 ps 1.5  lc rgb "#999999"  #grau, aufr. gschl. Dreieck


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

  # integral der Standardnormalverteilung
phi(x)            =norm(x)
phiOld(x)            =0.5*(erf(xlimited(x/sqrt(2)))+1)  #equivalent
                
 # Quantil der Standardnormalverteilung
phiQuantil(q)=invnorm(q)
phiQuantilOld(q)=sqrt(2.)*inverf(-1.+2*q)  #equivalent
                


lorenz(x)         = 1/(pi*(1+x**2))
expo(lambda,x)    = lambda*exp(-lambda*x)
gleich(a,b,x)     = ((x>=a)&&(x<=b)) ? 1/(b-a) : 0

#gamma(n+1)=n! for integer n>=0
# nu=Zahl der FG
studnorm(nu)  = gamma(0.5*(nu+1)) / (sqrt(nu*pi)*gamma(0.5*nu))
student(x,nu) =  studnorm(nu) / (1.+x**2/nu)**(0.5*(nu+1))
studentCum(x,nu) = ibeta(0.5*nu,0.5*nu,0.5*(x+sqrt(x**2+nu))/sqrt(x**2+nu))

betaFun(x,y)=gamma(x)*gamma(y)/(gamma(x+y))
fisher(x,d1,d2)=sqrt( (d1*x)**d1 * d2**d2/( (d1*x+d2)**(d1+d2)))\
 / (x*betaFun(0.5*d1, 0.5*d2))
fisherCum(x,d1,d2)=ibeta(0.5*d1, 0.5*d2, d1*x/(d1*x+d2))


chinorm(nu)   = 1./(2**(0.5*nu)*gamma(0.5*nu))
chi2(x,nu)    =  chinorm(nu) *exp(-0.5*x) * x**(0.5*nu-1)

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

####################################################
####################################################
# Daten 
# y=MIV-Weglaenge pro Tag
# x1=Einwohner/km^2
# x2=Dummy D=0, USA=1

n=10.
J=2.

t_nm1_0975=t9_0975
t_nm2_0975=t8_0975
t_nm3_0975=t7_0975

x11=140.;	x21=0;		y1=35;
x15=320.;	x25=0;		y5=24;
x12=180.;	x22=0;		y2=32;
x14=300.;	x24=0;		y4=28;
x13=240.;	x23=0;		y3=33;
x16=30;		x26=1;		y6=55;
x17=180.;	x27=1;		y7=44;
x18=70.;		x28=1;		y8=51;
x19=130.;	x29=1;		y9=41;
x110=95.;	x210=1;		y10=50;

# nur Pseudo, um danach Summen nicht aendern zu muessen
x111=0.;		x211=0;		y11=0;
x112=0.;		x212=0;		y12=0;

x1datamin=30
x1datamax=350
x2datamin=-0.1
x2datamax=1.1
ydatamin=20
ydatamax=60

####################################################
####################################################

x0bar=1.
x1bar=(x11+x12+x13+x14+x15+x16+x17+x18+x19+x110+x111+x112)/n
x2bar=(x21+x22+x23+x24+x25+x26+x27+x28+x29+x210+x211+x212)/n
ybar=(y1+y2+y3+y4+y5+y6+y7+y8+y9+y10+y11+y12)/n

# 3 x 3 - Matrix (\m{X}\tr\cdot\m{X})_{jk}=sumjk

sum00=n
sum01=n*x1bar
sum02=n*x2bar

sum11=x11**2+x12**2+x13**2+x14**2+x15**2+x16**2+x17**2+x18**2+x19**2+x110**2+x111**2+x112**2

sum12=x11*x21+x12*x22+x13*x23+x14*x24+x15*x25+x16*x26+x17*x27+x18*x28+x19*x29+x110*x210+x111*x211+x112*x212

sum22=x21**2+x22**2+x23**2+x24**2+x25**2+x26**2+x27**2+x28**2+x29**2+x210**2+x211**2+x212**2

sum10=sum01
sum20=sum02
sum21=sum12


# 3-Vektor (\m{X}\tr\cdot \vec{y})_j=sumjy

sum0y=n*ybar
sum1y=x11*y1+x12*y2+x13*y3+x14*y4+x15*y5+x16*y6+x17*y7+x18*y8+x19*y9+x110*y10+x111*y11+x112*y12
sum2y=x21*y1+x22*y2+x23*y3+x24*y4+x25*y5+x26*y6+x27*y7+x28*y8+x29*y9+x210*y10+x211*y11+x212*y12

# Quadratsumme y

sumyy=y1**2+y2**2+y3**2+y4**2+y5**2+y6**2+y7**2+y8**2+y9**2+y10**2+y11**2+y12**2

# "klassische" Kovarianzmatrix der "richtigen" Variablen x1,x2 und Kovarianzvektor y

s11=sum11/n - x1bar**2
s22=sum22/n - x2bar**2
s12=sum12/n - x1bar*x2bar
s21=s12
s1y=sum1y/n - x1bar*ybar
s2y=sum2y/n - x2bar*ybar



# Matrixinverse (X'X)^{-1}

detXX=sum00*(sum11*sum22-sum21*sum12)\
        - sum01*(sum10*sum22-sum20*sum12)\
        + sum02*(sum10*sum21-sum20*sum11)

inv00=(sum11*sum22-sum12*sum21)/detXX
inv01=(sum02*sum21-sum01*sum22)/detXX
inv02=(sum01*sum12-sum02*sum11)/detXX

inv11=(sum00*sum22-sum02*sum20)/detXX
inv12=(sum02*sum10-sum00*sum12)/detXX

inv22=(sum00*sum11-sum01*sum10)/detXX

inv10=inv01
inv20=inv02
inv21=inv12

#########################################################
print ""
print "Daten: Mittelwerte und Matrizen"
print "========================"
#########################################################

print "\nx1bar=",x1bar
print "x2bar=",x2bar
print "ybar=",ybar

print "\nMatrix X'X"
print "(",sum00,"\t",sum01,"\t",sum02,")"
print "(",sum10,"\t",sum11,"\t",sum12,")"
print "(",sum20,"\t",sum21,"\t",sum22,")"

print "\nMatrix (X'X)^{-1}"
print "(",inv00,"\t",inv01,"\t",inv02,")"
print "(",inv10,"\t",inv11,"\t",inv12,")"
print "(",inv20,"\t",inv21,"\t",inv22,")"

print "\nVektor X'y"
print "(",sum0y,"\t",sum1y,"\t",sum2y,")^T"

print "\nKovarianzmatrix exogene Var:"
print "s11=",s11," s12=",s12," s22=",s22

print "\nKovarianzvektor exogen-endogene Var:"
print "s1y=",s1y," s2y=",s2y


#########################################################
print "\nParametersch\"atzer und stat. Eigenschaften"
print "=================================="
#########################################################

beta0=inv00*sum0y+inv01*sum1y+inv02*sum2y
beta1=inv10*sum0y+inv11*sum1y+inv12*sum2y
beta2=inv20*sum0y+inv21*sum1y+inv22*sum2y

beta1Test=(s1y*s22-s2y*s12)/(s11*s22-s12**2) # mit klass Methode
beta2Test=(s2y*s11-s1y*s21)/(s11*s22-s12**2)

s_1y=sum1y/n - x1bar*ybar
s_2y=sum2y/n - x2bar*ybar

s_yy=sumyy/n - ybar**2   # Deskriptive Gesamtvarianz
hatsig_yy=n/(n-1) *s_yy    # Induktive Gesamtvarianz

s_hatyhaty=beta1*s_1y+beta2*s_2y       # Deskr. erkl. Var (OHNE beta0, da s0y=0!)
s_epseps=s_yy-s_hatyhaty                       # Deskr. Residualvarianz
hatsig_epseps=n/(n-1-J) * s_epseps       # Indukt. Residualvarianz

B=s_hatyhaty/s_yy                               # Deskriptives Best.-ma3
Bquer(B,n,J)=1-(n-1)/(n-1-J) * (1-B)     # Induktives Best.-ma3

covbeta00=hatsig_epseps*inv00
covbeta01=hatsig_epseps*inv01
covbeta02=hatsig_epseps*inv02
covbeta11=hatsig_epseps*inv11
covbeta12=hatsig_epseps*inv12
covbeta22=hatsig_epseps*inv22
covbeta10=covbeta01
covbeta20=covbeta02
covbeta21=covbeta12

r_beta0beta1=covbeta01/sqrt(covbeta00*covbeta11)
r_beta0beta2=covbeta02/sqrt(covbeta00*covbeta22)
r_beta1beta2=covbeta12/sqrt(covbeta11*covbeta22)

###########

print "\nVektor beta"
print "(",beta0,"\t",beta1,"\t",beta2,")^T"
print "Test klassisch: beta1Test=",beta1Test," beta2Test=",beta2Test

print "\nInduktiver Gesamtvarianzschaetzer hatsig_yy=",hatsig_yy
print "Induktiver Residualvarianzschaetzer hatsig_epseps=",hatsig_epseps
print "\nDeskriptives Best.-ma3 B=",B
print "Induktives Best.-ma3 Bquer=",Bquer(B,n,J)
print "\nCovarianzmatrix  hatsig_epseps*(X'X)^{-1} der beta-Schaetzer"
print "(",covbeta00,"\t",covbeta01,"\t",covbeta02,")"
print "(",covbeta10,"\t",covbeta11,"\t",covbeta12,")"
print "(",covbeta20,"\t",covbeta21,"\t",covbeta22,")"

print "Korrelationen"

print "r_beta0beta1=",r_beta0beta1," r_beta0beta2=",r_beta0beta2
print "r_beta1beta2=",r_beta1beta2

#########################################################
print "\nKonfidenzintervalle zu alpha=0.05"
print "======================================"
#########################################################

dbeta0=t_nm3_0975*sqrt(covbeta00)
dbeta1=t_nm3_0975*sqrt(covbeta11)
dbeta2=t_nm3_0975*sqrt(covbeta22)

print "Quantil t^(",n-3,")_{0.975}=",t_nm3_0975
print "Halbe Breiten der KI und KI selbst:"
print "dbeta0=",dbeta0,"   KI(beta_0)=[",beta0-dbeta0,",",beta0+dbeta0,"]"
print "dbeta1=",dbeta1,"   KI(beta_1)=[",beta1-dbeta1,",",beta1+dbeta1,"]"
print "dbeta2=",dbeta2,"   KI(beta_2)=[",beta2-dbeta2,",",beta2+dbeta2,"]"



#########################################################
print "\nSymmetrische Tests auf Parameterwerte=0: t- und p-Werte"
print "======================================"
t0=abs(beta0/sqrt(covbeta00))
t1=abs(beta1/sqrt(covbeta11))
t2=abs(beta2/sqrt(covbeta22))
p0=2*(1-studentCum(t0,n-3))
p1=2*(1-studentCum(t1,n-3))
p2=2*(1-studentCum(t2,n-3))
print "beta0: t-Wert=",t0," p-Wert=",p0
print "beta1: t-Wert=",t1," p-Wert=",p1
print "beta2: t-Wert=",t2," p-Wert=",p2

#########################################################
print "\nKonfidenzintervalle der Linearkombination U=beta2+c beta1"
print "======================================"
#########################################################
 
c1=dbeta2/dbeta1
c2=-dbeta2/dbeta1
erwU1=beta2+c1*beta1
erwU2=beta2+c2*beta1
varU1=covbeta22 + c1**2*covbeta11+2*c1*covbeta12
varU2=covbeta22 + c2**2*covbeta11+2*c2*covbeta12
dU1=t_nm3_0975*sqrt(varU1)
dU2=t_nm3_0975*sqrt(varU2)
print "U1=beta2+c1*beta1 = beta2 +dbeta2/dbeta1*beta1"
print "U2=beta2+c2*beta1 = beta2 -dbeta2/dbeta1*beta1"
print "Halbe Breiten von KI(U1/2) und KI(U1/2) selbst"
print "dU1=",dU1,"   KI(U1)=[",erwU1-dU1,",",erwU1+dU1,"]"
print "dU2=",dU2,"   KI(U2)=[",erwU2-dU2,",",erwU2+dU2,"]"


#########################################################
print "\nTest der verbundenen Nullhypothese beta1=beta10, beta2=beta20"
print "==================================================="
#########################################################

beta10=0.
beta20=0.

# Nullmodell yNull=betaNull definiert durch
# y (x1,x2;betaNull)=betaNull+beta10*x1+beta20*x2

#y-"Messwerte" des effektiven Trivialmodells yNull=betaNull
yNull1=y1-beta10*x11-beta20*x21
yNull2=y2-beta10*x12-beta20*x22
yNull3=y3-beta10*x13-beta20*x23
yNull4=y4-beta10*x14-beta20*x24
yNull5=y5-beta10*x15-beta20*x25
yNull6=y6-beta10*x16-beta20*x26
yNull7=y7-beta10*x17-beta20*x27
yNull8=y8-beta10*x18-beta20*x28
yNull9=y9-beta10*x19-beta20*x29
yNull10=y10-beta10*x110-beta20*x210
yNull11=y11-beta10*x111-beta20*x211
yNull12=y12-beta10*x112-beta20*x212

sumyNullyNull=yNull1**2+yNull2**2+yNull3**2+yNull4**2+yNull5**2+yNull6**2+yNull7**2+yNull8**2+yNull9**2+yNull10**2+yNull11**2+yNull12**2

betaNull=ybar-beta10*x1bar-beta20*x2bar

FminNull=sumyNullyNull-n*betaNull**2
Fmin=n*s_epseps
f_realis1=(FminNull-Fmin)/Fmin * (n-3.)/(3.-1.)

print  "\nH01: beta10=",beta10," beta20=",beta20," f_realis1=",f_realis1
print "Fehlerquadratsumme volles Modell: Fmin=",Fmin
print "Kalibrierter Parameter des Nullmodells: betaNull=",betaNull
print "Fehlerquadratsumme Nullmodell: FminNull=",FminNull

print  "f_realis1=",f_realis1," sim F(3-1,n-3) distributed"
print "fisherCum(f_realis1,2,n-3)=",fisherCum(f_realis1,2,n-3)
print "p-Wert=1-fisherCum(f_realis1,2,n-3)=",1-fisherCum(f_realis1,2,n-3)
print "Veransch: (FminNull-Fmin)=",FminNull-Fmin," mit zwei FG entsorgt"
print " => ",0.5*(FminNull-Fmin)," pro FG"
print " Es verbleibt Fmin/(n-3)=",Fmin/(n-3)," Fehlerquadratsumme pro FG"

##############################



#########################################################
print ""
print "Vergleich: Modelle mit nur einer exog. Var x1 bzw. x2"
print "==================================================="
#########################################################

s_11=sum11/n - x1bar**2 
s_22=sum22/n - x2bar**2 

beta1_M1=s_1y/s_11
beta2_M2=s_2y/s_22

s_hatyhaty_M1=beta1_M1*s_1y
s_hatyhaty_M2=beta2_M2*s_2y

sig2_est_M1=(s_yy-s_hatyhaty_M1)*n/(n-2)
sig2_est_M2=(s_yy-s_hatyhaty_M2)*n/(n-2)

covbeta11_M1 = 1./n * sig2_est_M1/s_11
covbeta22_M2 = 1./n * sig2_est_M2/s_22

dbeta1_M1=t_nm2_0975*sqrt(covbeta11_M1)
dbeta2_M2=t_nm2_0975*sqrt(covbeta22_M2)

F_M1=(s_yy-s_hatyhaty_M1)*n
F_M2=(s_yy-s_hatyhaty_M2)*n

fRealis_M1vsM=(F_M1-Fmin)/Fmin * (n-3)/1.
fRealis_M2vsM=(F_M2-Fmin)/Fmin * (n-3)/1.
p_M1vsM=1-fisherCum(fRealis_M1vsM,1,n-3)
p_M2vsM=1-fisherCum(fRealis_M2vsM,1,n-3)
print "Modell M1: Sterne allein:"
print "beta1_M1=",beta1_M1, " KI=[",beta1_M1-dbeta1_M1,",",beta1_M1+dbeta1_M1,"]"
print "Modell M2: Preis allein:"
print "beta2_M2=",beta2_M2, " KI=[",beta2_M2-dbeta2_M2,",",beta2_M2+dbeta2_M2,"]"

print "Fisher-Test M1 vs volles Modell M: F_M1=",F_M1
print "  fRealis_M1=",fRealis_M1vsM," p_M1vsM=",p_M1vsM
print "Fisher-Test M2 vs volles Modell M: F_M2=",F_M2
print "  fRealis_M2=",fRealis_M2vsM," p_M2vsM=",p_M2vsM

#########################################################
#########################################################



#######################################
print "plotting weglaenge_Ftest.eps"
set out "weglaenge_Ftest.eps"
#######################################
set param
set key bottom left
set xlabel "f"
set ylabel "Kumulierte Fisher-F-Verteilung F^{2,n-3}(f)"
xmax=14.
set xrange [0:xmax]
plot[t=0:1]\
  xmax*t, fisherCum(xmax*t,2,n-3) t "Kumulierte Fisher-Verteilung" w l ls 12,\
  f_realis1, t\
   t "Realisierter f-Wert bei {/Symbol b}_{10}=0,{/Symbol b}_{20}=0" w l ls 17,\
  t*xmax, 0.95 t "F=0.95" w l ls 1
#  f_realis2, t\
#   t "Realisierter f-Wert bei {/Symbol b}_{10}=34,{/Symbol b}_{20}=-1" w l ls 15,\

#######################################
print "plotting weglaenge_f_hatbeta1.eps"
set out "weglaenge_f_hatbeta1.eps"
#######################################

set title "Konfidenzintervall ({/Symbol b}_1) f\374r {/Symbol a}=0.05"
tmin=-1.5*t_nm3_0975
tmax=2.1*t_nm3_0975
sigbeta1=sqrt(covbeta11)
sigbeta2=sqrt(covbeta22)
beta1min=beta1+tmin*sigbeta1
beta1max=beta1+tmax*sigbeta1
beta2min=beta2+tmin*sigbeta2
beta2max=beta2+tmax*sigbeta2

studentNoNorm(x,mux,hatsigx,nu)=1./hatsigx*student( (x-mux)/hatsigx,nu)

set xlabel "{/Symbol b}_1 - Sch\344tzer"
set xrange [beta1min:beta1max]
set ylabel "Dichtefunktion f"
set param
set key at screen 0.94,0.78



plot [t=tmin:tmax]\
  beta1+t*sigbeta1,studentNoNorm(beta1+t*sigbeta1,beta1,sigbeta1,n-3)\
    t "Dichte" w l ls 12,\
  beta1-dbeta1 + (t-tmin)/(tmax-tmin)*2*dbeta1,0 t "KI" w l ls 29


#######################################
print "plotting weglaenge_F_hatbeta1.eps"
set out "weglaenge_F_hatbeta1.eps"
#######################################

set ylabel "Verteilungsfunktion F"

plot [t=tmin:tmax]\
  beta1+t*sigbeta1,studentCum(t,n-3) t "F" w l ls 12,\
  beta1-dbeta1 + (t-tmin)/(tmax-tmin)*2*dbeta1,0 t "KI" w l ls 29,\
  beta1-dbeta1, (t-tmin)/(tmax-tmin) t "" w l ls 1,\
  beta1+dbeta1, (t-tmin)/(tmax-tmin) t "" w l ls 1,\
  beta1+t*sigbeta1,0.025 t "F=0.025" w l ls 16,\
  beta1+t*sigbeta1,0.975 t "F=0.975" w l ls 17

#######################################
print "plotting weglaenge_f_hatbeta2.eps"
set out "weglaenge_f_hatbeta2.eps"
#######################################

set title "Konfidenzintervall zu H_0: {/Symbol b}_2 und {/Symbol s} wie gesch\344tzt und {/Symbol a}=0.05"

set xlabel "{/Symbol b}_2 - Sch\344tzer"
set xrange [beta2min:beta2max]
set ylabel "Dichtefunktion f "
plot [t=tmin:tmax]\
  beta2+t*sigbeta2,studentNoNorm(beta2+t*sigbeta2,beta2,sigbeta2,n-3)\
    t "Dichte" w l ls 12,\
  beta2-dbeta2 + (t-tmin)/(tmax-tmin)*2*dbeta2,0 t "KI" w l ls 29

#######################################
print "plotting weglaenge_F_hatbeta2.eps"
set out "weglaenge_F_hatbeta2.eps"
#######################################

set ylabel "Verteilungsfunktion F"

plot [t=tmin:tmax]\
  beta2+t*sigbeta2,studentCum(t,n-3) t "F" w l ls 12,\
  beta2-dbeta2 + (t-tmin)/(tmax-tmin)*2*dbeta2,0 t "KI" w l ls 29,\
  beta2-dbeta2, (t-tmin)/(tmax-tmin) t "" w l ls 1,\
  beta2+dbeta2, (t-tmin)/(tmax-tmin) t "" w l ls 1,\
  beta2+t*sigbeta2,0.025 t "F=0.025" w l ls 16,\
  beta2+t*sigbeta2,0.975 t "F=0.975" w l ls 17


#######################################
print "plotting weglaenge_scatterplot_x1x2.eps"
set out "weglaenge_scatterplot_x1x2.eps"
#######################################

set xlabel "Exogene Variable x_1"
set ylabel "Exogene Variable x_2"
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
  x110,x210 w p ls 1,\
  x111,x211 w p ls 1,\
  x112,x212 w p ls 1

#######################################
print "plotting weglaenge_scatterplot_x1y.eps"
set out "weglaenge_scatterplot_x1y.eps"
#######################################

set xlabel "Exogene Variable x_1"
set ylabel "Endogene Variable y"
set key at screen 0.95,0.95

set yrange [ydatamin:ydatamax]
set notitle
haty(x1,x2)=beta0+beta1*x1+beta2*x2
hatymod(x1,x2)=min(ydatamax,max(ydatamin,beta0+beta1*x1+beta2*x2))
hatymod2(x1,x2)=max(ydatamin,beta0+beta1*x1+beta2*x2)

plot[x1=x1datamin:x1datamax]\
  x11,y1 t "Daten" w p ls 1,\
  x12,y2 t "" w p ls 1,\
  x13,y3 t "" w p ls 1,\
  x14,y4 t "" w p ls 1,\
  x15,y5 t "" w p ls 1,\
  x16,y6 t "" w p ls 1,\
  x17,y7 t "" w p ls 1,\
  x18,y8 t "" w p ls 1,\
  x19,y9 t "" w p ls 1,\
  x110,y10 t "" w p ls 1,\
  x111,y11 t "" w p ls 1,\
  x112,y12 t "" w p ls 1,\
  x1, ybar+beta1_M1*(x1-x1bar) t "Einfachregression (x_1)" w l ls 12,\
  x11,haty(x11,x21) t "Mehrfachregression (x_1,x_2)" w p ls 7,\
  x12,haty(x12,x22) t "" w p ls 7,\
  x13,haty(x13,x23) t "" w p ls 7,\
  x14,haty(x14,x24) t "" w p ls 7,\
  x15,haty(x15,x25) t "" w p ls 7,\
  x16,haty(x16,x26) t "" w p ls 7,\
  x17,haty(x17,x27) t "" w p ls 7,\
  x18,haty(x18,x28) t "" w p ls 7,\
  x19,haty(x19,x29) t "" w p ls 7,\
  x110,haty(x110,x210) t "" w p ls 7,\
  x111,haty(x111,x211) t "" w p ls 7,\
  x112,haty(x112,x212) t "" w p ls 7

#######################################
print "plotting weglaenge_scatterplot_x2y.eps"
set out "weglaenge_scatterplot_x2y.eps"
#######################################

set key at screen 0.60,0.95

set xlabel "Exogene Variable x_2"

set xrange [x2datamin:x2datamax]
set yrange [ydatamin:ydatamax]
set notitle


plot[x2=x2datamin:x2datamax]\
  x21,y1 t "Daten" w p ls 1,\
  x22,y2 t "" w p ls 1,\
  x23,y3 t "" w p ls 1,\
  x24,y4 t "" w p ls 1,\
  x25,y5 t "" w p ls 1,\
  x26,y6 t "" w p ls 1,\
  x27,y7 t "" w p ls 1,\
  x28,y8 t "" w p ls 1,\
  x29,y9 t "" w p ls 1,\
  x210,y10 t "" w p ls 1,\
  x211,y11 t "" w p ls 1,\
  x212,y12 t "" w p ls 1,\
  x2, ybar+beta2_M2*(x2-x2bar) t "Einfachregression (x_2)" w l ls 12,\
  x21,haty(x11,x21) t "Mehrfachregression (x_1,x_2)" w p ls 7,\
  x22,haty(x12,x22) t "" w p ls 7,\
  x23,haty(x13,x23) t "" w p ls 7,\
  x24,haty(x14,x24) t "" w p ls 7,\
  x25,haty(x15,x25) t "" w p ls 7,\
  x26,haty(x16,x26) t "" w p ls 7,\
  x27,haty(x17,x27) t "" w p ls 7,\
  x28,haty(x18,x28) t "" w p ls 7,\
  x29,haty(x19,x29) t "" w p ls 7,\
  x210,haty(x110,x210) t "" w p ls 7,\
  x211,haty(x111,x211) t "" w p ls 7,\
  x212,haty(x112,x212) t "" w p ls 7


#######################################
print "plotting weglaenge_f2_hatbeta1_hatbeta2.eps"
set out "weglaenge_f2_hatbeta1_hatbeta2.eps"
#######################################

set term post eps enhanced color solid "Helvetica" 18

set param
set view map
set multiplot
set contour surface
set cntrparam bspline
set cntrparam levels 10 
set isosamples 30,30
set palette defined ( 0 "white", 5 "yellow", 30 "orange",\
  70 "red",  100 "#99003") 

xmin=beta1-3*sigbeta1
#xmin=0
#xmax=beta1+4*sigbeta1
xmax=0

#ymin=beta2-4*sigbeta2
ymin=0
ymax=beta2+3*sigbeta2
#ymax=0

set xlabel "{/Symbol b}_1 - Sch\344tzer"
set xrange [xmin:xmax]

set ylabel "{/Symbol b}_2 - Sch\344tzer"
set yrange [ymin:ymax]

set pm3d
unset clabel
#set title "2d-Dichte der {/Symbol b}-Sch \344tzer unter H_0: {/Symbol s} und {/Symbol b}_j wie gemessen"

set nokey
splot[x=xmin:xmax][y=ymin:ymax]\
  x,y,gauss2d(x-beta1,y-beta2,sigbeta1,sigbeta2,r_beta1beta2)\
   t "f_2(hat({/Symbol b})_1, hat({/Symbol b})_2)" w l ls 98

set cntrparam levels discrete 0.25  # eine Contour!
set key at screen 0.90,0.96
unset pm3d; unset surface
splot[x=xmin:xmax][y=ymin:ymax]\
  x,y,gauss2d(x-beta1,y-beta2,sigbeta1,sigbeta2,r_beta1beta2)\
   t "Konfidenzregion F-Test" w l ls 17

# Falls gruenes Rechteck statt vier gruene schneidende Linien, diese
# Zeilen in folg. Plotbefehl einsetzten
#  beta1-dbeta1,beta2-dbeta2+2*t*dbeta2,100+0.01*t\
#       t "Konfidenzregion t-Test" w l ls 15,\
#  beta1+dbeta1,beta2-dbeta2+2*t*dbeta2,100+0.01*t t "" w l ls 15,\
#  beta1-dbeta1+2*t*dbeta1,beta2-dbeta2,100+0.01*t t "" w l ls 15,\
#  beta1-dbeta1+2*t*dbeta1,beta2+dbeta2,100+0.01*t t "" w l ls 15,\

set surface
set key at screen 0.46,0.96
splot[t=0:1]\
 beta1-dbeta1, ymin+t*(ymax-ymin),100+0.01*t\
      t "Konfidenzintervalle T-Test {/Symbol b_1, b_2}" w l ls 15,\
 beta1+dbeta1, ymin+t*(ymax-ymin),100+0.01*t t "" w l ls 15,\
 xmin+t*(xmax-xmin), beta2-dbeta2,100+0.01*t t "" w l ls 15,\
 xmin+t*(xmax-xmin), beta2+dbeta2,100+0.01*t t "" w l ls 15
#  beta10,beta20,0.01 t "Nullhypothese 1" w p ls 1
set nomultiplot
set notitle

#######################################
print "plotting weglaenge_scatter3d_1.eps"
set out "weglaenge_scatter3d_1.eps"
#######################################
set nogrid
unset colorbox

set param
set multiplot
set contour surface
set cntrparam bspline
unset cntrparam; set cntrparam levels 20 
unset clabel
set isosamples 15,15
set palette defined ( 0 "white", 5 "yellow", 50 "orange",\
  80 "#FF6666",  99 "#FF44BB", 100 "white") 



unset surface
set pm3d  
set contour surface

set pm3d hidden3d 98

set view 50,160

set xlabel "x_1"
set xrange [x1datamin:x1datamax] reverse
set ylabel "x_2"
set yrange [x2datamin:x2datamax]
#set zlabel "y" offset 6,4
set label 1 "y" at x1datamin, x2datamin, ydatamax+0.2*(ydatamax-ydatamin)
set zrange [ydatamin:ydatamax]

set key at screen 0.78,0.97
splot[x1=x1datamin:x1datamax][x2=x2datamin:x2datamax]\
  x1,x2,hatymod(x1,x2)\
   t "{/Symbol b}_0+{/Symbol b}_1 x_1+{/Symbol b}_2 x_2" w l ls 98


set key at screen 0.78,0.92
unset pm3d; 
set surface
splot\
 x11,x21,y1 t "Daten" w p ls 1,\
 x12,x22,y2 t "" w p ls 1,\
 x13,x23,y3 t "" w p ls 1,\
 x14,x24,y4 t "" w p ls 1,\
 x15,x25,y5 t "" w p ls 1,\
 x16,x26,y6 t "" w p ls 1,\
 x17,x27,y7 t "" w p ls 1,\
 x18,x28,y8 t "" w p ls 1,\
 x19,x29,y9 t "" w p ls 1,\
 x110,x210,y10 t "" w p ls 1,\
 x111,x211,y11 t "" w p ls 1,\
 x112,x212,y12 t "" w p ls 1


unset pm3d; 
set surface
set key at screen 0.78,0.87
splot[t=0:1]\
 x11,x21,haty(x11,x21)+t*(y1-haty(x11,x21)) t "{/Symbol e}_i" w l ls 17,\
 x12,x22,haty(x12,x22)+t*(y2-haty(x12,x22)) t "" w l ls 17,\
 x13,x23,haty(x13,x23)+t*(y3-haty(x13,x23)) t "" w l ls 17,\
 x14,x24,haty(x14,x24)+t*(y4-haty(x14,x24)) t "" w l ls 17,\
 x15,x25,haty(x15,x25)+t*(y5-haty(x15,x25)) t "" w l ls 17,\
 x16,x26,haty(x16,x26)+t*(y6-haty(x16,x26)) t "" w l ls 17,\
 x17,x27,haty(x17,x27)+t*(y7-haty(x17,x27)) t "" w l ls 17,\
 x18,x28,haty(x18,x28)+t*(y8-haty(x18,x28)) t "" w l ls 17,\
 x19,x29,haty(x19,x29)+t*(y9-haty(x19,x29)) t "" w l ls 17,\
 x110,x210,haty(x110,x210)+t*(y10-haty(x110,x210)) t "" w l ls 17,\
 x111,x211,haty(x111,x211)+t*(y11-haty(x111,x211)) t "" w l ls 17,\
 x112,x212,haty(x112,x212)+t*(y12-haty(x112,x212)) t "" w l ls 17

set nomultiplot

#######################################
print "plotting weglaenge_scatter3d_2.eps"
set out "weglaenge_scatter3d_2.eps"
#######################################

set multiplot

set view 31,65
set xrange [x1datamin:x1datamax] noreverse


set key at screen 0.78,0.97
splot[x1=x1datamin:x1datamax][x2=x2datamin:x2datamax]\
  x1,x2,hatymod2(x1,x2)\
   t "{/Symbol b}_0+{/Symbol b}_1 x_1+{/Symbol b}_2 x_2" w l ls 98


set key at screen 0.78,0.92
unset pm3d; 
set surface
splot\
 x11,x21,y1 t "Daten" w p ls 1,\
 x12,x22,y2 t "" w p ls 1,\
 x13,x23,y3 t "" w p ls 1,\
 x14,x24,y4 t "" w p ls 1,\
 x15,x25,y5 t "" w p ls 1,\
 x16,x26,y6 t "" w p ls 1,\
 x17,x27,y7 t "" w p ls 1,\
 x18,x28,y8 t "" w p ls 1,\
 x19,x29,y9 t "" w p ls 1,\
 x110,x210,y10 t "" w p ls 1,\
 x111,x211,y11 t "" w p ls 1,\
 x112,x212,y12 t "" w p ls 1


unset pm3d; 
set surface
set key at screen 0.78,0.87
splot[t=0:1]\
 x11,x21,haty(x11,x21)+t*(y1-haty(x11,x21)) t "{/Symbol e}_i" w l ls 17,\
 x12,x22,haty(x12,x22)+t*(y2-haty(x12,x22)) t "" w l ls 17,\
 x13,x23,haty(x13,x23)+t*(y3-haty(x13,x23)) t "" w l ls 17,\
 x14,x24,haty(x14,x24)+t*(y4-haty(x14,x24)) t "" w l ls 17,\
 x15,x25,haty(x15,x25)+t*(y5-haty(x15,x25)) t "" w l ls 17,\
 x16,x26,haty(x16,x26)+t*(y6-haty(x16,x26)) t "" w l ls 17,\
 x17,x27,haty(x17,x27)+t*(y7-haty(x17,x27)) t "" w l ls 17,\
 x18,x28,haty(x18,x28)+t*(y8-haty(x18,x28)) t "" w l ls 17,\
 x19,x29,haty(x19,x29)+t*(y9-haty(x19,x29)) t "" w l ls 17,\
 x110,x210,haty(x110,x210)+t*(y10-haty(x110,x210)) t "" w l ls 17,\
 x111,x211,haty(x111,x211)+t*(y11-haty(x111,x211)) t "" w l ls 17,\
 x112,x212,haty(x112,x212)+t*(y12-haty(x112,x212)) t "" w l ls 17

set nomultiplot
unset label 1


print "\nKlausuraufgaben"
print "==========="

print"\n(b)"
beta_costUSA=-3;
beta_costD=-0.5;

ybarD=(y1+y2+y3+y4+y5)/5.
ybarUSA=(y6+y7+y8+y9+y10)/5.
kmCentD=20.
kmCentUSA=6.

epsilonUSA=kmCentUSA/ybarUSA*beta_costUSA
epsilonD=kmCentD/ybarD*beta_costD
print "Reiseweite D=",ybarD," Reiseweite USA=",ybarUSA
print "kmCentD=",kmCentD," kmCentUSA=",kmCentUSA
print "epsilonD=",epsilonD
print "epsilonUSA=",epsilonUSA


print"\n(c)"
print "\nKovarianzmatrix exogene Var:"
print "s11=",s11," s12=",s12," s22=",s22

print "\nKovarianzvektor exogen-endogene Var:"
print "s1y=",s1y," s2y=",s2y

detS=s11*s22-s12**2
beta1=(s1y*s22-s2y*s12)/detS# mit klass Methode
beta2=(s2y*s11-s1y*s21)/detS
beta0=ybar-beta1*x1bar-beta2*x2bar
print "detS=",detS
print "beta1=",beta1
print "beta2=",beta2
print "beta0=",beta0

