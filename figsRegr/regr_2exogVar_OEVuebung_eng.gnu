
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
# (Line- and Punkttypen Oct. 2009)
#geordnet nach hue (ps)
# BUG: Manchmal
#     obskure Abhaengigkeit of the Symbole von Reihenfolge bzw linestyles
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
set style line 13 lt 8 lw 6 pt 3 ps 1.2 #blassrot, offener star
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

betaFun(x,y)=gamma(x)*gamma(y)/(gamma(x+y))
fisher(x,d1,d2)=sqrt( (d1*x)**d1 * d2**d2/( (d1*x+d2)**(d1+d2)))\
 / (x*betaFun(0.5*d1, 0.5*d2))
fisherCum(x,d1,d2)=ibeta(0.5*d1, 0.5*d2, d1*x/(d1*x+d2))



max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x


set term post eps enhanced color solid "Helvetica" 22

####################################################
####################################################
# Data Hotel-Beispiel
# x1=Zahl of the stars
# x2=Preis (Euro/Nacht
# y=Auslastung (%)

n=5.
J=2.
alpha=0.05

t_nm1_0975=tQuantil(1-0.5*alpha,n-1)
t_nm2_0975=tQuantil(1-0.5*alpha,n-2)
t_nm3_0975=tQuantil(1-0.5*alpha,n-3)
t_nm3_095=tQuantil(1-alpha,n-3)


x11=1.8;        x21=20; y1=220;
x12=3.5;        x22=30; y2=180;
x13=2.5;        x23=25; y3=200;
x14=1.7;        x24=25; y4=240;
x15=2.5;        x25=20; y5=160;

x16=0.;         x26=0;  y6=0;
x17=0.;         x27=0;  y7=0;
x18=0.;         x28=0;  y8=0;
x19=0.;         x29=0;  y9=0;
x110=0.;        x210=0;  y10=0;
x111=0.;        x211=0;  y11=0;
x112=0.;        x212=0;  y12=0;


x1datamin=1.4
x1datamax=3.7
x2datamin=18
x2datamax=32
ydatamin=140
ydatamax=260
ydatamax3d=310

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

# Errorquadratsumme S(beta)=y'y-2 beta' (X'y)+beta'(X'X)beta
S(beta0,beta1,beta2)=sumyy\
  -2*(beta0*sum0y+beta1*sum1y+beta2*sum2y)\
  + beta0**2*sum00 + beta1**2*sum11 + beta2**2*sum22\
  + 2*(beta0*beta1*sum01 + beta0*beta2*sum02 + beta1*beta2*sum12)



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
print "Data: Mittelwerte and Matrizen"
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

#########################################################
print "\nParametersch\"atzer and stat. Eigenschaften"
print "=================================="
#########################################################

beta0=inv00*sum0y+inv01*sum1y+inv02*sum2y
beta1=inv10*sum0y+inv11*sum1y+inv12*sum2y
beta2=inv20*sum0y+inv21*sum1y+inv22*sum2y

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

print "\nDeskriptive Residualvarianz s_epseps=",s_epseps
print "Deskriptive erklaerte Variance s_hatyhaty=",s_hatyhaty
print "Deskriptive Gesamtvarianz s_yy=",s_yy
print "Deskriptives Best.-ma3 B=",B

print "\nInduktiver Residualvarianzschaetzer hatsig_epseps=",hatsig_epseps
print "Induktiver Gesamtvarianzschaetzer hatsig_yy=",hatsig_yy
print "Induktives Best.-ma3 Bquer=",Bquer(B,n,J)

print "\nKovarianzmatrix  hatsig_epseps*(X'X)^{-1} of the beta-Schaetzer"
print "(",covbeta00,"\t",covbeta01,"\t",covbeta02,")"
print "(",covbeta10,"\t",covbeta11,"\t",covbeta12,")"
print "(",covbeta20,"\t",covbeta21,"\t",covbeta22,")"

print "\nStandard Deviationen of the beta-Schaetzer:"
print "sqrt(covbeta00)=",sqrt(covbeta00)
print "sqrt(covbeta11)=",sqrt(covbeta11)
print "sqrt(covbeta22)=",sqrt(covbeta22)

print "\nKorrelationen of the beta-Schaetzer:"
print "r_beta0beta1=",r_beta0beta1," r_beta0beta2=",r_beta0beta2
print "r_beta1beta2=",r_beta1beta2

#########################################################
print "\nConfidence Intervals for alpha=0.05"
print "======================================"
#########################################################

dbeta0=t_nm3_0975*sqrt(covbeta00)
dbeta1=t_nm3_0975*sqrt(covbeta11)
dbeta2=t_nm3_0975*sqrt(covbeta22)

print "Quantil t^(",n-3,")_{0.975}=",t_nm3_0975
print "Quantil t^(",n-3,")_{0.95}=",t_nm3_095
print "Halbe Breiten of the KI and KI selbst:"
print "beta0=",beta0
print "     dbeta0=",dbeta0,"  KI(beta_0)=[",beta0-dbeta0,",",beta0+dbeta0,"]"
print "beta1=",beta1
print "     dbeta1=",dbeta1,"  KI(beta_1)=[",beta1-dbeta1,",",beta1+dbeta1,"]"
print "beta2=",beta2
print "     dbeta2=",dbeta2,"  KI(beta_2)=[",beta2-dbeta2,",",beta2+dbeta2,"]"



#########################################################
print "\nSymmetrische Tests auf Parameterwerte=0: t- and p-Valuee"
print "======================================"
t0=abs(beta0/sqrt(covbeta00))
t1=abs(beta1/sqrt(covbeta11))
t2=abs(beta2/sqrt(covbeta22))
p0=2*(1-studentCum(t0,n-3))
p1=2*(1-studentCum(t1,n-3))
p2=2*(1-studentCum(t2,n-3))
print "beta0: t-Value=",t0," p-Value=",p0
print "beta1: t-Value=",t1," p-Value=",p1
print "beta2: t-Value=",t2," p-Value=",p2
#########################################################
print "\nAsymmetrischer Test auf beta2<-1.5: t- and p-Valuee"
beta20=-1.5
t2asym=(beta2-beta20)/sqrt(covbeta22)
p2asym=1*(1-studentCum(t2asym,n-3))
print "beta2asym: t-Value=",t2asym," p-Value=",p2asym


#########################################################
print "\nConfidence Intervals of the Linearkombination U=beta2+c beta1"
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
print "Halbe Breiten von KI(U1/2) and KI(U1/2) selbst"
print "dU1=",dU1,"   KI(U1)=[",erwU1-dU1,",",erwU1+dU1,"]"
print "dU2=",dU2,"   KI(U2)=[",erwU2-dU2,",",erwU2+dU2,"]"

#########################################################
print "\nTest of the verbundenen Nullhypothese beta1=beta10, beta2=beta20"
print "==================================================="
#########################################################

beta101=30.
beta20=-1.

# Nullmodell yNull=betaNull definiert durch
# y (x1,x2;betaNull)=betaNull+beta10*x1+beta20*x2

#y-"Messwerte" of the effektiven Trivialmodells yNull=betaNull
yNull1=y1-beta101*x11-beta20*x21
yNull2=y2-beta101*x12-beta20*x22
yNull3=y3-beta101*x13-beta20*x23
yNull4=y4-beta101*x14-beta20*x24
yNull5=y5-beta101*x15-beta20*x25
yNull6=y6-beta101*x16-beta20*x26
yNull7=y7-beta101*x17-beta20*x27
yNull8=y8-beta101*x18-beta20*x28
yNull9=y9-beta101*x19-beta20*x29
yNull10=y10-beta101*x110-beta20*x210
yNull11=y11-beta101*x111-beta20*x211
yNull12=y12-beta101*x112-beta20*x212

sumyNullyNull=yNull1**2+yNull2**2+yNull3**2+yNull4**2+yNull5**2+yNull6**2+yNull7**2+yNull8**2+yNull9**2+yNull10**2+yNull11**2+yNull12**2

betaNull=ybar-beta101*x1bar-beta20*x2bar

FminNull=sumyNullyNull-n*betaNull**2
Fmin=n*s_epseps
f_realis1=(FminNull-Fmin)/Fmin * (n-3.)/(3.-1.)

print  "\nH01: beta101=",beta101," beta20=",beta20," f_realis1=",f_realis1
print "Errorquadratsumme volles Modell: Fmin=",Fmin
print "Kalibrierter Parameter of the Nullmodells: betaNull=",betaNull
print "Errorquadratsumme Nullmodell: FminNull=",FminNull

print  "f_realis1=",f_realis1," sim F(3-1,n-3) distributed"
print "fisherCum(f_realis1,2,n-3)=",fisherCum(f_realis1,2,n-3)
print "p-Value=1-fisherCum(f_realis1,2,n-3)=",1-fisherCum(f_realis1,2,n-3)
print "Veransch: (FminNull-Fmin)=",FminNull-Fmin," zwei DOF entsorgt"
print " => ",0.5*(FminNull-Fmin)," pro DOF"
print " Es verbleibt Fmin/(n-3)=",Fmin/(n-3)," Errorquadratsumme pro DOF"

##############################
beta102=34.
beta20=-1.
yNull1=y1-beta102*x11-beta20*x21
yNull2=y2-beta102*x12-beta20*x22
yNull3=y3-beta102*x13-beta20*x23
yNull4=y4-beta102*x14-beta20*x24
yNull5=y5-beta102*x15-beta20*x25
yNull6=y6-beta102*x16-beta20*x26
yNull7=y7-beta102*x17-beta20*x27
yNull8=y8-beta102*x18-beta20*x28
yNull9=y9-beta102*x19-beta20*x29
yNull10=y10-beta102*x110-beta20*x210
yNull11=y11-beta102*x111-beta20*x211
yNull12=y12-beta102*x112-beta20*x212

sumyNullyNull=yNull1**2+yNull2**2+yNull3**2+yNull4**2+yNull5**2+yNull6**2+yNull7**2+yNull8**2+yNull9**2+yNull10**2+yNull11**2+yNull12**2

betaNull=ybar-beta102*x1bar-beta20*x2bar

FminNull=sumyNullyNull-n*betaNull**2
Fmin=n*s_epseps
f_realis2=(FminNull-Fmin)/Fmin * (n-3.)/(3.-1.)
print  "\nH02: beta102=",beta102," beta20=",beta20,":"
print " Fmin=", Fmin," FminNull=",FminNull," f_realis2=",f_realis2
print " p-Value=1-fisherCum(f_realis2,2,n-3)=",1-fisherCum(f_realis2,2,n-3)

##############################
beta103=30.
beta203=-0.5
yNull1=y1-beta103*x11-beta203*x21
yNull2=y2-beta103*x12-beta203*x22
yNull3=y3-beta103*x13-beta203*x23
yNull4=y4-beta103*x14-beta203*x24
yNull5=y5-beta103*x15-beta203*x25
yNull6=y6-beta103*x16-beta203*x26
yNull7=y7-beta103*x17-beta203*x27
yNull8=y8-beta103*x18-beta203*x28
yNull9=y9-beta103*x19-beta203*x29
yNull10=y10-beta103*x110-beta203*x210
yNull11=y11-beta103*x111-beta203*x211
yNull12=y12-beta103*x112-beta203*x212

sumyNullyNull=yNull1**2+yNull2**2+yNull3**2+yNull4**2+yNull5**2+yNull6**2+yNull7**2+yNull8**2+yNull9**2+yNull10**2+yNull11**2+yNull12**2

betaNull=ybar-beta103*x1bar-beta203*x2bar

FminNull=sumyNullyNull-n*betaNull**2
Fmin=n*s_epseps
f_realis3=(FminNull-Fmin)/Fmin * (n-3.)/(3.-1.)
print  "\nH03: beta103=",beta103," beta203=",beta203,":"
print " Fmin=", Fmin," FminNull=",FminNull," f_realis3=",f_realis3
print " p-Value=1-fisherCum(f_realis3,2,n-3)=",1-fisherCum(f_realis3,2,n-3)




#########################################################
print ""
print "Vergleich: Modelle nur einer exog. Var x1 bzw. x2"
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
print "Modell M1: stars allein:"
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
print "plotting PTuebung_Ftest_eng.eps"
set out "PTuebung_Ftest_eng.eps"
#######################################
set param
set key bottom left
set xlabel "f"
set ylabel "Cumulated Fisher-F-Distribution F^{2,n-3}(f)"
xmax=14.
set xrange [0:xmax]
plot[t=0:1]\
  xmax*t, fisherCum(xmax*t,2,n-3) t "Cumulated Fisher-Distribution" w l ls 12,\
  f_realis1, t\
   t "Realized f-Value bei {/Symbol b}_{10}=30,{/Symbol b}_{20}=-1" w l ls 17,\
  f_realis2, t\
   t "Realized f-Value bei {/Symbol b}_{10}=34,{/Symbol b}_{20}=-1" w l ls 15,\
  t*xmax, 0.95 t "F=0.95" w l ls 1

#######################################
print "plotting PTuebung_f_hatbeta1_eng.eps"
set out "PTuebung_f_hatbeta1_eng.eps"
#######################################

set title "Confidence Interval ({/Symbol b}_1) f\374r {/Symbol a}=0.05"
tmin=-1.5*t_nm3_0975
tmax=2.1*t_nm3_0975
sigbeta1=sqrt(covbeta11)
sigbeta2=sqrt(covbeta22)
beta1min=beta1+tmin*sigbeta1
beta1max=beta1+tmax*sigbeta1
beta2min=beta2+tmin*sigbeta2
beta2max=beta2+tmax*sigbeta2

studentNoNorm(x,mux,hatsigx,nu)=1./hatsigx*student( (x-mux)/hatsigx,nu)

set xlabel "{/Symbol b}_1 - Estimator"
set xrange [beta1min:beta1max]
set ylabel "Density f"
set param
set key at screen 0.94,0.78



plot [t=tmin:tmax]\
  beta1+t*sigbeta1,studentNoNorm(beta1+t*sigbeta1,beta1,sigbeta1,n-3)\
    t "Density" w l ls 12,\
  beta1-dbeta1 + (t-tmin)/(tmax-tmin)*2*dbeta1,0 t "KI" w l ls 29


#######################################
print "plotting PTuebung_F_hatbeta1_eng.eps"
set out "PTuebung_F_hatbeta1_eng.eps"
#######################################

set ylabel "Distribution Function F"

plot [t=tmin:tmax]\
  beta1+t*sigbeta1,studentCum(t,n-3) t "F" w l ls 12,\
  beta1-dbeta1 + (t-tmin)/(tmax-tmin)*2*dbeta1,0 t "KI" w l ls 29,\
  beta1-dbeta1, (t-tmin)/(tmax-tmin) t "" w l ls 1,\
  beta1+dbeta1, (t-tmin)/(tmax-tmin) t "" w l ls 1,\
  beta1+t*sigbeta1,0.025 t "F=0.025" w l ls 16,\
  beta1+t*sigbeta1,0.975 t "F=0.975" w l ls 17

#######################################
print "plotting PTuebung_f_hatbeta2_eng.eps"
set out "PTuebung_f_hatbeta2_eng.eps"
#######################################

set title "Confidence Interval for H_0: {/Symbol b}_2 and {/Symbol s} as estimated and {/Symbol a}=0.05"

set xlabel "{/Symbol b}_2 - Estimator"
set xrange [beta2min:beta2max]
set ylabel "Density f "
plot [t=tmin:tmax]\
  beta2+t*sigbeta2,studentNoNorm(beta2+t*sigbeta2,beta2,sigbeta2,n-3)\
    t "Density" w l ls 12,\
  beta2-dbeta2 + (t-tmin)/(tmax-tmin)*2*dbeta2,0 t "KI" w l ls 29

#######################################
print "plotting PTuebung_F_hatbeta2_eng.eps"
set out "PTuebung_F_hatbeta2_eng.eps"
#######################################

set ylabel "Distribution Function F"

plot [t=tmin:tmax]\
  beta2+t*sigbeta2,studentCum(t,n-3) t "F" w l ls 12,\
  beta2-dbeta2 + (t-tmin)/(tmax-tmin)*2*dbeta2,0 t "KI" w l ls 29,\
  beta2-dbeta2, (t-tmin)/(tmax-tmin) t "" w l ls 1,\
  beta2+dbeta2, (t-tmin)/(tmax-tmin) t "" w l ls 1,\
  beta2+t*sigbeta2,0.025 t "F=0.025" w l ls 16,\
  beta2+t*sigbeta2,0.975 t "F=0.975" w l ls 17


#######################################
print "plotting PTuebung_scatterplot_x1x2_eng.eps"
set out "PTuebung_scatterplot_x1x2_eng.eps"
#######################################

set xlabel "Exogenous Variable x_1"
set ylabel "Exogenous Variable x_2"
set nokey
set xrange [x1datamin:x1datamax]
set yrange [x2datamin:x2datamax]
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
print "plotting PTuebung_scatterplot_x1y_eng.eps"
set out "PTuebung_scatterplot_x1y_eng.eps"
#######################################

set xlabel "Exogenous Variable x_1"
set ylabel "Endogenous Variable y"
set key at screen 0.95,0.95

set yrange [ydatamin:ydatamax]
set notitle
haty(x1,x2)=beta0+beta1*x1+beta2*x2

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
  x110,y10 t "" w p ls 1,\
  x111,y11 t "" w p ls 1,\
  x112,y12 t "" w p ls 1,\
  x1, ybar+beta1_M1*(x1-x1bar) t "Simple Regression (x_1)" w l ls 12,\
  x11,haty(x11,x21) t "Multiple Regression (x_1,x_2)" w p ls 7,\
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
print "plotting PTuebung_scatterplot_x2y_eng.eps"
set out "PTuebung_scatterplot_x2y_eng.eps"
#######################################

set key at screen 0.60,0.95

set xlabel "Exogenous Variable x_2"

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
  x210,y10 t "" w p ls 1,\
  x211,y11 t "" w p ls 1,\
  x212,y12 t "" w p ls 1,\
  x2, ybar+beta2_M2*(x2-x2bar) t "Simple Regression (x_2)" w l ls 12,\
  x21,haty(x11,x21) t "Multiple Regression (x_1,x_2)" w p ls 7,\
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
print "plotting PTuebung_scatterplot_x2ycondx1_eng.eps"
set out "PTuebung_scatterplot_x2ycondx1_eng.eps"
#######################################
# Nur sinnvoll bei Hotelbeispiel ?!
min1=17
max1=42
min2=30
max2=60
min3=55
max3=100
min4=80
max4=x2datamax

plot[t=0:1]\
  x21,y1 t "" w p ls 1,\
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
  x2datamin+t*(x2datamax-x2datamin),\
    ybar+beta2_M2*(x2datamin+t*(x2datamax-x2datamin)-x2bar) \
    t "Simple Regression (x_2)" w l ls 11,\
  min1+t*(max1-min1),haty(1, min1+t*(max1-min1) ) t "y est(x_1=1 star, x_2)" w l ls 12,\
  min2+t*(max2-min2),haty(2,min2+t*(max2-min2)) t "y est (x_1=2 stars, x_2)" w l ls 13,\
  min3+t*(max3-min3),haty(3,min3+t*(max3-min3)) t "y est (x_1=3 stars, x_2)" w l ls 14,\
  min4+t*(max4-min4),haty(4,min4+t*(max4-min4)) t "y est (x_1=4 stars, x_2)" w l ls 15

#######################################
print "plotting PTuebung_scatterplot_x2ycondx1_notCal1_eng.eps"
set out "PTuebung_scatterplot_x2ycondx1_notCal1_eng.eps"
#######################################
# Nur sinnvoll bei Hotelbeispiel ?!

beta1_notCal1=beta1+dbeta1
beta2_notCal1=beta2+dbeta2
haty_notCal1(x1,x2)=beta0+beta1_notCal1*x1+beta2_notCal1*x2
set title "{/Symbol b}_1 and {/Symbol b}_2 um {/Symbol Db_1} bzw. {/Symbol Db_2} verschoben"
set key at screen 0.43,0.20
plot[t=0:1]\
  x21,y1 t "" w p ls 2,\
  x22,y2 t "" w p ls 2,\
  x23,y3 t "" w p ls 2,\
  x24,y4 t "" w p ls 3,\
  x25,y5 t "" w p ls 3,\
  x26,y6 t "" w p ls 3,\
  x27,y7 t "" w p ls 4,\
  x28,y8 t "" w p ls 4,\
  x29,y9 t "" w p ls 4,\
  x210,y10 t "" w p ls 5,\
  x211,y11 t "" w p ls 5,\
  x212,y12 t "" w p ls 5,\
  min1+t*(max1-min1),haty_notCal1(1, min1+t*(max1-min1) )\
      t "y est(x_1=1 star, x_2)" w l ls 12,\
  min2+t*(max2-min2),haty_notCal1(2,min2+t*(max2-min2))\
      t "y est (x_1=2 stars, x_2)" w l ls 13,\
  min3+t*(max3-min3),haty_notCal1(3,min3+t*(max3-min3))\
      t "y est (x_1=3 stars, x_2)" w l ls 14,\
  min4+t*(max4-min4),haty_notCal1(4,min4+t*(max4-min4))\
      t "y est (x_1=4 stars, x_2)" w l ls 15

#######################################
print "plotting PTuebung_scatterplot_x2ycondx1_notCal2_eng.eps"
set out "PTuebung_scatterplot_x2ycondx1_notCal2_eng.eps"
#######################################
# Nur sinnvoll bei Hotelbeispiel ?!

set title "{/Symbol b}_1 and {/Symbol b}_2 um {/Symbol Db_1} bzw. - {/Symbol Db_2} verschoben"
beta1_notCal2=beta1+dbeta1
beta2_notCal2=beta2-dbeta2
haty_notCal2(x1,x2)=beta0+beta1_notCal2*x1+beta2_notCal2*x2

plot[t=0:1]\
  x21,y1 t "" w p ls 2,\
  x22,y2 t "" w p ls 2,\
  x23,y3 t "" w p ls 2,\
  x24,y4 t "" w p ls 3,\
  x25,y5 t "" w p ls 3,\
  x26,y6 t "" w p ls 3,\
  x27,y7 t "" w p ls 4,\
  x28,y8 t "" w p ls 4,\
  x29,y9 t "" w p ls 4,\
  x210,y10 t "" w p ls 5,\
  x211,y11 t "" w p ls 5,\
  x212,y12 t "" w p ls 5,\
  min1+t*(max1-min1),haty_notCal2(1, min1+t*(max1-min1) )\
      t "y est(x_1=1 star, x_2)" w l ls 12,\
  min2+t*(max2-min2),haty_notCal2(2,min2+t*(max2-min2))\
      t "y est (x_1=2 stars, x_2)" w l ls 13,\
  min3+t*(max3-min3),haty_notCal2(3,min3+t*(max3-min3))\
      t "y est (x_1=3 stars, x_2)" w l ls 14,\
  min4+t*(max4-min4),haty_notCal2(4,min4+t*(max4-min4))\
      t "y est (x_1=4 stars, x_2)" w l ls 15

#######################################
print "plotting PTuebung_scatterplot_x2ycondx1_notCal3_eng.eps"
set out "PTuebung_scatterplot_x2ycondx1_notCal3_eng.eps"
#######################################
# Nur sinnvoll bei Hotelbeispiel ?!

set title "{/Symbol b}_1 and {/Symbol b}_2 um - {/Symbol Db_1} bzw. - {/Symbol Db_2} verschoben"
beta1_notCal3=beta1-dbeta1
beta2_notCal3=beta2-dbeta2
haty_notCal3(x1,x2)=beta0+beta1_notCal3*x1+beta2_notCal3*x2

plot[t=0:1]\
  x21,y1 t "" w p ls 2,\
  x22,y2 t "" w p ls 2,\
  x23,y3 t "" w p ls 2,\
  x24,y4 t "" w p ls 3,\
  x25,y5 t "" w p ls 3,\
  x26,y6 t "" w p ls 3,\
  x27,y7 t "" w p ls 4,\
  x28,y8 t "" w p ls 4,\
  x29,y9 t "" w p ls 4,\
  x210,y10 t "" w p ls 5,\
  x211,y11 t "" w p ls 5,\
  x212,y12 t "" w p ls 5,\
  min1+t*(max1-min1),haty_notCal3(1, min1+t*(max1-min1) )\
      t "y est(x_1=1 star, x_2)" w l ls 12,\
  min2+t*(max2-min2),haty_notCal3(2,min2+t*(max2-min2))\
      t "y est (x_1=2 stars, x_2)" w l ls 13,\
  min3+t*(max3-min3),haty_notCal3(3,min3+t*(max3-min3))\
      t "y est (x_1=3 stars, x_2)" w l ls 14,\
  min4+t*(max4-min4),haty_notCal3(4,min4+t*(max4-min4))\
      t "y est (x_1=4 stars, x_2)" w l ls 15

#######################################
print "plotting PTuebung_scatterplot_x2ycondx1_notCal4_eng.eps"
set out "PTuebung_scatterplot_x2ycondx1_notCal4_eng.eps"
#######################################
# Nur sinnvoll bei Hotelbeispiel ?!

set title "{/Symbol b}_1 and {/Symbol b}_2 um - {/Symbol Db_1} bzw. +{/Symbol Db_2} verschoben"

beta1_notCal4=beta1-dbeta1
beta2_notCal4=beta2+dbeta2
haty_notCal4(x1,x2)=beta0+beta1_notCal4*x1+beta2_notCal4*x2

plot[t=0:1]\
  x21,y1 t "" w p ls 2,\
  x22,y2 t "" w p ls 2,\
  x23,y3 t "" w p ls 2,\
  x24,y4 t "" w p ls 3,\
  x25,y5 t "" w p ls 3,\
  x26,y6 t "" w p ls 3,\
  x27,y7 t "" w p ls 4,\
  x28,y8 t "" w p ls 4,\
  x29,y9 t "" w p ls 4,\
  x210,y10 t "" w p ls 5,\
  x211,y11 t "" w p ls 5,\
  x212,y12 t "" w p ls 5,\
  min1+t*(max1-min1),haty_notCal4(1, min1+t*(max1-min1) )\
      t "y est(x_1=1 star, x_2)" w l ls 12,\
  min2+t*(max2-min2),haty_notCal4(2,min2+t*(max2-min2))\
      t "y est (x_1=2 stars, x_2)" w l ls 13,\
  min3+t*(max3-min3),haty_notCal4(3,min3+t*(max3-min3))\
      t "y est (x_1=3 stars, x_2)" w l ls 14,\
  min4+t*(max4-min4),haty_notCal4(4,min4+t*(max4-min4))\
      t "y est (x_1=4 stars, x_2)" w l ls 15
set notitle


#######################################
print "plotting PTuebung_f2_hatbeta1_hatbeta2_eng.eps"
set out "PTuebung_f2_hatbeta1_hatbeta2_eng.eps"
#######################################

unset label

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

xmin=beta1-6*sigbeta1
xmax=beta1+6*sigbeta1

ymin=beta2-6*sigbeta2
ymax=beta2+6*sigbeta2

set xlabel "{/Symbol b}_1 - Estimator"
set xrange [xmin:xmax]

set ylabel "{/Symbol b}_2 - Estimator"
set yrange [ymin:ymax]

set pm3d
unset clabel
set nokey
splot[x=xmin:xmax][y=ymin:ymax]\
  x,y,gauss2d(x-beta1,y-beta2,sigbeta1,sigbeta2,r_beta1beta2)\
   t "f_2(hat({/Symbol b})_1, hat({/Symbol b})_2)" w l ls 98

# eine contour fuer Confidence Region

set label "Confidence Region F-Test" at screen 0.5,0.4 #bug key
unset pm3d; unset surface
set cntrparam levels discrete 0.000001
splot[x=xmin:xmax][y=ymin:ymax]\
  x,y,gauss2d(x-beta1,y-beta2,sigbeta1,sigbeta2,r_beta1beta2)\
   t "" w l ls 17


set surface
#bug falls y-Value > 0.96
set key at screen 0.52,0.80
splot[t=0:1]\
 beta1-dbeta1, ymin+t*(ymax-ymin),100+0.01*t\
      t "Confidence Intervals T-Test {/Symbol b_1, b_2}" w l ls 15,\
 beta1+dbeta1, ymin+t*(ymax-ymin),100+0.01*t t "" w l ls 15,\
 xmin+t*(xmax-xmin), beta2-dbeta2,100+0.01*t t "" w l ls 15,\
 xmin+t*(xmax-xmin), beta2+dbeta2,100+0.01*t t "" w l ls 15

set nomultiplot
set notitle

#######################################
print "plotting PTuebung_S_hatbeta1_hatbeta2_eng.eps"
set out "PTuebung_S_hatbeta1_hatbeta2_eng.eps"
#######################################
#!!! ACHTUNG; hier nicht Rekurs auf 1d-Modell moegl, da fuer festes
# beta0 (im Ggs zur Gaussverteilung of the Wahrsch, wo alle Valuee beta0
# als multil Faktoren enthalten sind)

set term post eps enhanced color solid "Helvetica" 18

set param
set view map
set multiplot
set contour surface
unset cntrparam
set cntrparam bspline
set cntrparam levels 10 
set isosamples 30,30
set palette defined ( 0 "white", 5 "yellow", 30 "orange",\
  70 "red",  100 "#99003") 

xmin=beta1-12*sigbeta1
xmax=beta1+12*sigbeta1

ymin=beta2-8*sigbeta2
ymax=beta2+8*sigbeta2

set xlabel "{/Symbol b}_1 - Estimator"
set xrange [xmin:xmax]

set ylabel "{/Symbol b}_2 - Estimator"
set yrange [ymin:ymax]

set zrange [0:10*S(beta0,0,0)]
set cbrange [0:10*S(beta0,0,0)]
set cntrparam levels incr 0,10000,100000
set pm3d
unset clabel

set nokey
splot[x=xmin:xmax][y=ymin:ymax]\
  x,y,S(beta0,x,y)\
   t "f_2(hat({/Symbol b})_1, hat({/Symbol b})_2)" w l ls 98

unset pm3d; unset surface


set surface
set key at screen 0.42,0.96
splot[t=0:1]\
 beta1-dbeta1, ymin+t*(ymax-ymin),100+0.01*t\
      t "Confidence Intervals T-Test {/Symbol b_1, b_2}" w l ls 15,\
 beta1+dbeta1, ymin+t*(ymax-ymin),100+0.01*t t "" w l ls 15,\
 xmin+t*(xmax-xmin), beta2-dbeta2,100+0.01*t t "" w l ls 15,\
 xmin+t*(xmax-xmin), beta2+dbeta2,100+0.01*t t "" w l ls 15

set nomultiplot
set notitle

#######################################
print "plotting PTuebung_fisherCum_eng.eps"
set out "PTuebung_fisherCum_eng.eps"
#######################################
set term post eps enhanced color solid "Helvetica" 20
set size 1,0.7
fmin=1
fmax=5
set param
set nokey
set xlabel "f"
set xrange [fmin:fmax]
set ylabel "F_{Fisher(2,n-3)}"
set auto y
plot[f=fmin:fmax]\
 f, fisherCum(f,2,n-3) w l ls 12
set size 1,1


#######################################
print "plotting PTuebung_studentCum_eng.eps"
set out "PTuebung_studentCum_eng.eps"
#######################################
set term post eps enhanced color solid "Helvetica" 20
set size 1,0.7
tmin=1
tmax=8
set param
set nokey
set xlabel "t"
set xrange [tmin:tmax]
set ylabel "F_{T(n-3)}"
set yrange [0.8:1]
plot[t=tmin:tmax]\
 t, studentCum(t,n-3) w l ls 12

#######################################
print "plotting PTuebung_studentCumSolved_eng.eps"
set out "PTuebung_studentCumSolved_eng.eps"
#######################################
set key at screen 0.5,0.2
plot[t=tmin:tmax]\
 t, studentCum(t,n-3) t "" w l ls 12,\
 t1, 0.8+0.2*(t-tmin)/(tmax-tmin) t "t_1" w l ls 15,\
 t2asym, 0.8+0.2*(t-tmin)/(tmax-tmin) t "t_2" w l ls 17



#######################################
print "plotting PTuebung_scatter3d_1_eng.eps"
set out "PTuebung_scatter3d_1_eng.eps"
#######################################
set term post eps enhanced color solid "Helvetica" 18

hatymod(x1,x2)=min(ydatamax3d,max(ydatamin,beta0+beta1*x1+beta2*x2))
hatymod2(x1,x2)=max(ydatamin,beta0+beta1*x1+beta2*x2)

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
set cbrange [ydatamin:ydatamax3d]


unset surface
set pm3d  
set contour surface
set pm3d hidden3d 98

set view 50,150

set xlabel "x_1"
set xrange [x1datamin:x1datamax] reverse
set ylabel "x_2"
set yrange [x2datamin:x2datamax]
#set zlabel "y" offset 6,4
set label 1 "y" at x1datamin, x2datamin, ydatamax3d+0.2*(ydatamax3d-ydatamin)
set zrange [ydatamin:ydatamax3d]

set key at screen 0.78,0.67
splot[x1=x1datamin:x1datamax][x2=x2datamin:x2datamax]\
  x1,x2,hatymod(x1,x2)\
   t "" w l ls 98


set key at screen 0.78,0.62
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
 x110,x210,y10 t "" w p ls 1,\
 x111,x211,y11 t "" w p ls 1,\
 x112,x212,y12 t "" w p ls 1


unset pm3d; 
set surface
set key at screen 0.78,0.57
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
print "plotting PTuebung_scatter3d_2_eng.eps"
set out "PTuebung_scatter3d_2_eng.eps"
#######################################

set multiplot

set view 50,60
set xrange [x1datamin:x1datamax] noreverse
set cbrange [ydatamin:ydatamax3d]

unset surface
set pm3d  
set contour surface
set pm3d hidden3d 98

set key at screen 0.78,0.67
splot[x1=x1datamin:x1datamax][x2=x2datamin:x2datamax]\
  x1,x2,hatymod2(x1,x2)\
   t "" w l ls 98


set key at screen 0.78,0.57
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
 x110,x210,y10 t "" w p ls 1,\
 x111,x211,y11 t "" w p ls 1,\
 x112,x212,y12 t "" w p ls 1


unset pm3d; 
set surface
splot[t=0:1]\
 x11,x21,haty(x11,x21)+t*(y1-haty(x11,x21)) t "" w l ls 17,\
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
