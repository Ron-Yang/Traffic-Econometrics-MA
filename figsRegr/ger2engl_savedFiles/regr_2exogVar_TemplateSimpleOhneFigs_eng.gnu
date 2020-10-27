
# (feb12) aus 
# ~/vorlesungen/Verkehrsoekonometrie_Ma/skript/figsRegr/regr_2exogVarSimpleOhneFigs.gnu

set encoding iso_8859_1   # dann äöüßÄÖÜ durch woertl. Eingabe korrekt
#\304     "A
#\326     "O
#\334     "U
#\344     "a
#\366     "o
#\374     "u
#\337     &szlig;

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



max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x


 ############################
 # Gauss-Distributionen und Derivate
 ############################


gauss(mu,sigma,x) = exp(-(x-mu)**2/(2* sigma**2)) / (sigma * sqrt(2*pi))
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
 # Student-t-Distribution  und Derivate, nu=Zahl der FG
 ############################

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
 # Chi-Quadrat-Distribution und Derivate, nu=Zahl der FG
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

#Test
p=0.999; nu=2; t0=nu
print "\nTest tQuantil und chi2Quantil fuer p=",p," nu=",nu
print "\ntQuantil(p,nu)=",tQuantil(p,nu)
print "studentCum(tQuantil(p,nu), nu)=",studentCum(tQuantil(p,nu),nu)
print "\nchi2Quantil(p,nu)=",chi2Quantil(p,nu)
print "chi2Cum(chi2Quantil(p,nu), nu)=",chi2Cum(chi2Quantil(p,nu),nu)


##################
# Fisher's F Distribution
##################

betaFun(x,y)=gamma(x)*gamma(y)/(gamma(x+y))
fisher(x,d1,d2)=sqrt( (d1*x)**d1 * d2**d2/( (d1*x+d2)**(d1+d2)))\
 / (x*betaFun(0.5*d1, 0.5*d2))
fisherCum(x,d1,d2)=ibeta(0.5*d1, 0.5*d2, d1*x/(d1*x+d2))



set term post eps enhanced color solid "Helvetica" 22

####################################################
####################################################
# Daten 
# y=Car-Weglaenge pro Tag
# x1=Einwohner/km^2
# x2=Dummy D=0, USA=1

n=8.
J=2.

t_nm1_0975=tQuantil(0.975,n-1)
t_nm2_0975=tQuantil(0.975,n-2)
t_nm3_0975=tQuantil(0.975,n-3)

x11=65.;		x21=0;		y1=3;
x12=75.;		x22=6;		y2=7;
x13=210.;	x23=2;		y3=10;
x14=190.;	x24=5;		y4=12;
x15=34.;		x25=13;		y5=9;
x16=49;		x26=7;		y6=5;
x17=95.;		x27=6;		y7=6;
x18=90.;		x28=12;		y8=11;

# nur Pseudo, um danach Summen nicht aendern zu muessen

x19=0.;		x29=0;		y9=0;
x110=0.;		x210=0;		y10=0;
x111=0.;		x211=0;		y11=0;
x112=0.;		x212=0;		y12=0;

x1datamin=30
x1datamax=220
x2datamin=0
x2datamax=14
ydatamin=3
ydatamax=15

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
print "beta0=",beta0
print "beta1=",beta1
print "beta2=",beta2
print "Test klassisch: beta1Test=",beta1Test," beta2Test=",beta2Test

print "\nInduktiver Gesamtvarianzschaetzer hatsig_yy=",hatsig_yy
print "Induktiver Residualvarianzschaetzer hatsig_epseps=",hatsig_epseps
print "\nDeskriptives Best.-ma3 B=",B
print "Induktives Best.-ma3 Bquer=",Bquer(B,n,J)
print "\nCovarianzmatrix  hatsig_epseps*(X'X)^{-1} der beta-Schaetzer"
print "(",covbeta00,"\t",covbeta01,"\t",covbeta02,")"
print "(",covbeta10,"\t",covbeta11,"\t",covbeta12,")"
print "(",covbeta20,"\t",covbeta21,"\t",covbeta22,")"

print "Korrelationen beta"

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


print "\nSymmetrische Tests auf Parameterwerte=0: t- und p-Werte"
print "======================================"
t0=abs(beta0/sqrt(covbeta00))
t1=abs(beta1/sqrt(covbeta11))
t2=abs(beta2/sqrt(covbeta22))
p0=2*(1-studentCum(t0,n-3))
p1=2*(1-studentCum(t1,n-3))
p2=2*(1-studentCum(t2,n-3))
print "beta0=",beta0,": t-Wert=",t0," p-Wert=",p0
print "beta1=",beta1,": t-Wert=",t1," p-Wert=",p1
print "beta2=",beta2,": t-Wert=",t2," p-Wert=",p2



print ""
print "Vergleich: Modelle mit nur einer exog. Var x1 bzw. x2"
print "====================================="


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

print "Modell M1: x1 allein:"
print "beta1_M1=",beta1_M1, " KI=[",beta1_M1-dbeta1_M1,",",beta1_M1+dbeta1_M1,"]"
print "Modell M2: x2 allein:"
print "beta2_M2=",beta2_M2, " KI=[",beta2_M2-dbeta2_M2,",",beta2_M2+dbeta2_M2,"]"


print "\nKlausuraufgaben"
print "==========="

print"\n(a)"
print "\nKovarianzmatrix exogene Var:"
print "s11=",s11," s12=",s12," s22=",s22

print "\nKovarianzvektor exogen-endogene Var:"
print "s1y=",s1y," s2y=",s2y

print "\nx1bar=",x1bar
print "x2bar=",x2bar
print "ybar=",ybar
detS=s11*s22-s12**2
beta1=(s1y*s22-s2y*s12)/detS# mit klass Methode
beta2=(s2y*s11-s1y*s21)/detS
beta0=ybar-beta1*x1bar-beta2*x2bar
print "detS=",detS
print "beta1=",beta1
print "beta2=",beta2
print "beta0=",beta0

print"\n(b)"
epsilon1=x1bar/ybar*beta1
print "epsilon1 mit echtem Wert: epsilon1=",epsilon1
print "epsilon1 mit Wert der Aufgabenstellung: x1bar/ybar*0.05=",x1bar/ybar*0.05

print"\n(c)"
print "gegeben: Induktiver Residualvarianzschaetzer hatsig_epseps=",hatsig_epseps
r12=s12/sqrt(s11*s22)
V11=hatsig_epseps/(n*s11*(1-r12**2))
t1Test=beta1/sqrt(V11)
t1Test_betaAufg=0.05/sqrt(V11)
V22=hatsig_epseps/(n*s22*(1-r12**2))
t2Test=beta2/sqrt(V22)
t2Test_betaAufg=0.6/sqrt(V22)

print "r12=",r12
print "V11=",V11
print "V22=",V22
print "t1Test=",t1Test," mit beta-Wert aus Aufgabenstellung: t1_aufg=",t1Test_betaAufg
print "t2Test=",t2Test," mit beta-Wert aus Aufgabenstellung: t2_aufg=",t2Test_betaAufg
print "beta1=",beta1,": t-Wert=",t1," p-Wert=",p1
print "beta2=",beta2,": t-Wert=",t2," p-Wert=",p2


