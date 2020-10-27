
###########################
# Template:
# von ~/vorlesungen/Verkehrsoekonometrie_Ma/skript/figsRegr/regr_2exogVarMatrix.gnu
############################

set encoding iso_8859_1   # dann ������� durch woertl. Eingabe korrekt
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

# if independently color, dash-style, point-style gnuplot 4.2 and  later
# !!! BUG as of oct 2017: DOS in contours at first splot takes lt
# instead of ls need to assume lt 6 (black) at ls 98

set style line 98 lt 6 lw 0 linecolor rgb "#000000" # beliebige Farben:Schwarz
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
set style line 12 lt 1 lw 6 pt 2 ps 0.8  lc rgb "#CC0022" #rot, dash, Kreuz
set style line 13 lt 8 lw 6 pt 3 ps 1.2 #blassrot, offener Stern
set style line 14 lt 6 lw 6 pt 4 ps 1.5  lc rgb "#FFAA00"  #gelb, offenes Quadrat
set style line 15 lt 1 lw 6 pt 5 ps 0.8  lc rgb "#00DD22"  #gruen,solid,ClosedBox
set style line 16 lt 5 lw 6 pt 7 ps 0.8  lc rgb "#00AAAA" #offener Kreis
set style line 17 lt 1 lw 9 pt 7 ps 1.5  lc rgb "#1100AA"  #blau,solid,Bullet
set style line 18 lt 4 lw 6 pt 8 ps 1.5  lc rgb "#6600AA"  #lila, aufrechtes geschloss. Dreieck
set style line 19 lt 7 lw 6 pt 9 ps 1.5  lc rgb "#999999"  #grau, aufr. gschl. Dreieck


set style line 21 lt 1 lw 6 pt 7 ps 0.8  lc rgb "#000000" #schwarz,solid,bullet
set style line 22 lt 1 lw 20 pt 2 ps 1.5  lc rgb "#CC0022" #rot, dash, Kreuz
set style line 29 lt 7 lw 20 pt 9 ps 1.5  lc rgb "#999999"  #grau, aufr. gschl. Dreieck



############### Beispiele fuer Funktionen ####################



xlimited(x)       =(x<-5) ? -5 : ( (x>5) ? 5 : x)

  # Normierte Dichte der Standardnormalverteilung

gauss(mu,sigma,x) = exp(-(x-mu)**2/(2* sigma**2)) / (sigma*sqrt(2*pi))

  # (kum) Verteilungsfunktion der Standardnormalverteilung

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

# 2d-Dichten (bei Student einfach jede Koordinate wie student(1d)+Korrel.
#  transformiert)

gauss2d(x,y,sx,sy,r) \
 = exp(-(0.5/(1-r*r))*((x/sx)**2+(y/sy)**2-2*r*x*y/(sx*sy)))\
  / (2*pi*sx*sy*sqrt(1-r*r))
gauss2d_invCovMatrix(x1,x2,invs11, invs12,invs22)\
 = sqrt(invs11*invs22-invs12**2) /  (2*pi)\
 * exp(-0.5*(x1*invs11*x1+2*x1*invs12*x2+x2*invs22*x2))

# auf Max=1 normierte 2d-Studentverteilung

xgauss2xstudent(x,nu)=tQuantil(phi(x) ,nu)

student2d(x,y,r,nu) \
 =student(sqrt((1./(1-r*r))*(x**2+y**2-2*r*x*y)),nu)/student(0,nu)


max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x


set term post eps enhanced color solid "Helvetica" 20


####################################################
####################################################
# Daten Hotel-Beispiel
# x1=Zahl der Sterne
# x2=Preis (Euro/Nacht
# y=Auslastung (%)

n=12.
J=2.
alpha=0.05

t_nm1_0975=tQuantil(1-0.5*alpha,n-1)
t_nm2_0975=tQuantil(1-0.5*alpha,n-2)
t_nm3_0975=tQuantil(1-0.5*alpha,n-3)
t_nm3_095=tQuantil(1-alpha,n-3)

x11=1.;		x21=15;		y1=42;
x12=1.;		x22=31;		y2=38;
x13=1.;		x23=40;		y3=24;
x14=2.;		x24=34;		y4=76;
x15=2.;		x25=50;		y5=52;
x16=2.;		x26=58;		y6=40;
x17=3.;		x27=67;		y7=90;
x18=3.;		x28=72;		y8=77;
x19=3.;		x29=84;		y9=62;
x110=4.;	x210=82;	y10=90;
x111=4.;	x211=98;	y11=82;
x112=4.;	x212=116;	y12=68;

x1datamin=0.8
x1datamax=4.2
x2datamin=10
x2datamax=120
ydatamin=20
ydatamax=120

#Fehlerquadratsumme als f(beta)

SE(x1,x2,y,beta0,beta1,beta2)=(beta0+beta1*x1+beta2*x2-y)**2

SSE(beta0,beta1,beta2)\
 = SE(x11,x21,y1,beta0,beta1,beta2)\
 + SE(x12,x22,y2,beta0,beta1,beta2)\
 + SE(x13,x23,y3,beta0,beta1,beta2)\
 + SE(x14,x24,y4,beta0,beta1,beta2)\
 + SE(x15,x25,y5,beta0,beta1,beta2)\
 + SE(x16,x26,y6,beta0,beta1,beta2)\
 + SE(x17,x27,y7,beta0,beta1,beta2)\
 + SE(x18,x28,y8,beta0,beta1,beta2)\
 + SE(x19,x29,y9,beta0,beta1,beta2)\
 + SE(x110,x210,y10,beta0,beta1,beta2)\
 + SE(x111,x211,y11,beta0,beta1,beta2)\
 + SE(x112,x212,y12,beta0,beta1,beta2)


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

# Klassische Varianz-Kovarianz matrix exog-exog und Kov-Vektor exog-endog
 
s11=sum11/n - x1bar**2
s22=sum22/n - x2bar**2
s12=sum12/n - x1bar*x2bar
s21=s12
s1y=sum1y/n - x1bar*ybar
s2y=sum2y/n - x2bar*ybar
syy=sumyy/n - ybar**2   # Deskriptive Gesamtvarianz


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

print "\nKlassisch"
print "==========="

print "\nx1bar=",x1bar
print "x2bar=",x2bar
print "ybar=",ybar


print "\nKovarianzmatrix exogene Var:"
print "s11=",s11," s12=",s12," s22=",s22

print "\nKovarianzvektor exogen-endogene Var:"
print "s1y=",s1y," s2y=",s2y

print "\nVarianz endogene Var:"
print "syy=",syy

print "\nMit Matrizenformulierung
print "=========================="

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
print "\nParametersch\"atzer und stat. Eigenschaften"
print "=================================="
#########################################################

beta0=inv00*sum0y+inv01*sum1y+inv02*sum2y
beta1=inv10*sum0y+inv11*sum1y+inv12*sum2y
beta2=inv20*sum0y+inv21*sum1y+inv22*sum2y


SSEmin=SSE(beta0,beta1,beta2)

hatsig_yy=n/(n-1) *syy    # Induktive Gesamtvarianz

s_hatyhaty=beta1*s1y+beta2*s2y       # Deskr. erkl. Var (OHNE beta0, da s0y=0!)
s_epseps=syy-s_hatyhaty                       # Deskr. Residualvarianz
hatsig_epseps=n/(n-1-J) * s_epseps       # Indukt. Residualvarianz

B=s_hatyhaty/syy                               # Deskriptives Best.-ma3
Bquer(B,n,J)=1-(n-1)/(n-1-J) * (1-B)     # Induktives Best.-ma3

cov00=hatsig_epseps*inv00
cov01=hatsig_epseps*inv01
cov02=hatsig_epseps*inv02
cov11=hatsig_epseps*inv11
cov12=hatsig_epseps*inv12
cov22=hatsig_epseps*inv22
cov10=cov01
cov20=cov02
cov21=cov12

r_01=cov01/sqrt(cov00*cov11)
r_02=cov02/sqrt(cov00*cov22)
r_12=cov12/sqrt(cov11*cov22)

###########

print "\nVektor beta"
print "(",beta0,"\t",beta1,"\t",beta2,")^T"

print "\nMinimale SSE: SSEmin=",SSEmin
print "Deskriptive Gesamtvarianz syy=",syy
print "Deskriptive Residualvarianz s_epseps=",s_epseps
print "Deskriptive erklaerte Varianz s_hatyhaty=",s_hatyhaty
print "Deskriptives Best.-ma3 B=",B

print "\nInduktiver Gesamtvarianzschaetzer hatsig_yy=",hatsig_yy
print "Induktiver Residualvarianzschaetzer hatsig_epseps=",hatsig_epseps
print "Induktives Best.-ma3 Bquer=",Bquer(B,n,J)
print "\ngeschaetzte Covarianzmatrix  hatsig_epseps*(X'X)^{-1} der beta-Schaetzer"
print "(",cov00,"\t",cov01,"\t",cov02,")"
print "(",cov10,"\t",cov11,"\t",cov12,")"
print "(",cov20,"\t",cov21,"\t",cov22,")"

print "\nKorrelationen der beta-Schaetzer"

print "r_01=",r_01," r_02=",r_02
print "r_12=",r_12


#########################################################
print "\nKonfidenzintervalle zu alpha=0.05"
print "======================================"
#########################################################

dbeta0=t_nm3_0975*sqrt(cov00)
dbeta1=t_nm3_0975*sqrt(cov11)
dbeta2=t_nm3_0975*sqrt(cov22)

print "Quantil t^(",n-3,")_{0.975}=",t_nm3_0975
print "Quantil t^(",n-3,")_{0.95}=",t_nm3_095
print "Halbe Breiten der KI und KI selbst:"
print "beta0=",beta0
print "     dbeta0=",dbeta0,"  KI(beta_0)=[",beta0-dbeta0,",",beta0+dbeta0,"]"
print "beta1=",beta1
print "     dbeta1=",dbeta1,"  KI(beta_1)=[",beta1-dbeta1,",",beta1+dbeta1,"]"
print "beta2=",beta2
print "     dbeta2=",dbeta2,"  KI(beta_2)=[",beta2-dbeta2,",",beta2+dbeta2,"]"



#########################################################
print "\nSymmetrische Tests auf Parameterwerte=0: t- und p-Werte"
print "======================================"
t0=beta0/sqrt(cov00)
t1=beta1/sqrt(cov11)
t2=beta2/sqrt(cov22)
p0=2*(1-studentCum(abs(t0),n-3))
p1=2*(1-studentCum(abs(t1),n-3))
p2=2*(1-studentCum(abs(t2),n-3))

beta202=-1.5
t2asym=(beta2-beta202)/sqrt(cov22)
p2asym=1*(1-studentCum(t2asym,n-3))

print "beta0: t-Wert=",t0," p-Wert=",p0
print "beta1: t-Wert=",t1," p-Wert=",p1
print "beta2: t-Wert=",t2," p-Wert=",p2

print "\nAsymmetrische Tests auf Parameterwerte< oder >0: t- und p-Werte"
print "======================================"
print "H0: beta1<0: t=",t1," p=1-studentCum(t,n-3)=", 1-studentCum(t1,n-3)
print "H0: beta2>0: t=",t2," p=studentCum(t,n-3)=", studentCum(t2,n-3)
print "H0: beta2<",beta202,": t=",t2asym," p=",1-studentCum(t2asym,n-3)


#########################################################

c=30.
print "\nKombi-Test auf beta1<-",c,"*beta2 bzw. GamH0=beta1+",c,"*beta2<0"
print"=================================================================
erwGamH0=beta1+c*beta2
varGamH0=cov11+2*c*cov12+c**2*cov22
tGamH0=erwGamH0/sqrt(varGamH0)
print " erwGamH0=",erwGamH0," varGamH0=",varGamH0
print " tGamH0=erwGamH0/sqrt(varGamH0)=",tGamH0
print " pGamH0=1-studentCum(tGamH0,n-3)=", 1-studentCum(tGamH0,n-3)





#########################################################
c1=dbeta2/dbeta1        #dbeta1/2 sind KI von beta1/2, hier nur als Leitlinie
c2=-dbeta2/dbeta1
print "\nKonfidenzintervalle der Linearkombinationen"
print "   gamma1=beta2+",c1,"*beta1
print "   gamma2=beta2+",c2,"*beta1
print "================================================="

erwGam1=beta2+c1*beta1
erwGam2=beta2+c2*beta1
varGam1=cov22 + c1**2*cov11+2*c1*cov12
varGam2=cov22 + c2**2*cov11+2*c2*cov12
dGam1=t_nm3_0975*sqrt(varGam1)
dGam2=t_nm3_0975*sqrt(varGam2)
print "  erwGam1=",erwGam1," varGam1=",varGam1
print "  erwGam2=",erwGam2," varGam2=",varGam2
print "  Halbe Breiten von KI(Gam1/2) und KI(Gam1/2) selbst"
print "    dGam1=",dGam1,"   KI(Gam1)=[",erwGam1-dGam1,",",erwGam1+dGam1,"]"
print "    dGam2=",dGam2,"   KI(Gam2)=[",erwGam2-dGam2,",",erwGam2+dGam2,"]"



#########################################################
print "\nTest der verbundenen Nullhypothese beta1=beta10, beta2=beta20"
print "==================================================="
#########################################################

#H01

beta101=30.
beta201=-1.0
str_H01=sprintf("{/Symbol b}_{10}=%1.1f, {/Symbol b}_{20}=%1.1f",\
  beta101, beta201)
str_realH01=sprintf("Realisierter f-Wert bei %s",str_H01)


#H03

beta103=30.
beta203=-0.6
str_H03=sprintf("{/Symbol b}_{10}=%1.1f, {/Symbol b}_{20}=%1.1f",\
  beta103, beta203)
str_realH03=sprintf("Realisierter f-Wert bei %s",str_H03)
print "str_realH03=",str_realH03

# Nullmodell yNull=betaNull definiert durch
# y (x1,x2;betaNull)=betaNull+beta10*x1+beta20*x2

#y-"Messwerte" des effektiven Trivialmodells yNull=betaNull
yNull1=y1-beta101*x11-beta201*x21
yNull2=y2-beta101*x12-beta201*x22
yNull3=y3-beta101*x13-beta201*x23
yNull4=y4-beta101*x14-beta201*x24
yNull5=y5-beta101*x15-beta201*x25
yNull6=y6-beta101*x16-beta201*x26
yNull7=y7-beta101*x17-beta201*x27
yNull8=y8-beta101*x18-beta201*x28
yNull9=y9-beta101*x19-beta201*x29
yNull10=y10-beta101*x110-beta201*x210
yNull11=y11-beta101*x111-beta201*x211
yNull12=y12-beta101*x112-beta201*x212

sumyNullyNull=yNull1**2+yNull2**2+yNull3**2+yNull4**2+yNull5**2+yNull6**2+yNull7**2+yNull8**2+yNull9**2+yNull10**2+yNull11**2+yNull12**2

betaNull=ybar-beta101*x1bar-beta201*x2bar

FminNull=sumyNullyNull-n*betaNull**2
Fmin=n*s_epseps
f_realis1=(FminNull-Fmin)/Fmin * (n-3.)/(3.-1.)

print  "\nH01: beta101=",beta101," beta201=",beta201
print "  Fehlerquadratsumme volles Modell: Fmin=",Fmin
print "  Kalibrierter Parameter des Nullmodells: betaNull=",betaNull
print "  Fehlerquadratsumme Nullmodell: FminNull=",FminNull

print  "  f_realis1=",f_realis1," at H01 drawn from F(3-1,n-3) distribution"
print "  fisherCum(f_realis1,2,n-3)=",fisherCum(f_realis1,2,n-3)
print "  p=1-fisherCum(f_realis1,2,n-3)=",1-fisherCum(f_realis1,2,n-3)
print "  Veransch: (FminNull-Fmin)=",FminNull-Fmin," mit zwei FG entsorgt"
print "             => ",0.5*(FminNull-Fmin)," pro FG"
print "             rest Fmin/(n-3)=",Fmin/(n-3)," Fehlerquadratsumme pro FG"



# H03 (def siehe oben)

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

print  "\n\nH03: beta103=",beta103," beta203=",beta203
print "  Fehlerquadratsumme volles Modell: Fmin=",Fmin
print "  Kalibrierter Parameter des Nullmodells: betaNull=",betaNull
print "  Fehlerquadratsumme Nullmodell: FminNull=",FminNull

print  "  f_realis3=",f_realis3," at H03 drawn from F(3-1,n-3) distribution"
print "  fisherCum(f_realis3,2,n-3)=",fisherCum(f_realis3,2,n-3)
print "  p=1-fisherCum(f_realis3,2,n-3)=",1-fisherCum(f_realis3,2,n-3)
print "  Veransch: (FminNull-Fmin)=",FminNull-Fmin," mit zwei FG entsorgt"
print "             => ",0.5*(FminNull-Fmin)," pro FG"
print "             rest Fmin/(n-3)=",Fmin/(n-3)," Fehlerquadratsumme pro FG"





#########################################################
print ""
print "Vergleich: Modelle mit nur einer exog. Var x1 bzw. x2"
print "==================================================="
#########################################################


beta1_M1=s1y/s11
beta2_M2=s2y/s22

s_hatyhaty_M1=beta1_M1*s1y
s_hatyhaty_M2=beta2_M2*s2y

sig2_est_M1=(syy-s_hatyhaty_M1)*n/(n-2)
sig2_est_M2=(syy-s_hatyhaty_M2)*n/(n-2)

cov11_M1 = 1./n * sig2_est_M1/s11
cov22_M2 = 1./n * sig2_est_M2/s22

dbeta1_M1=t_nm2_0975*sqrt(cov11_M1)
dbeta2_M2=t_nm2_0975*sqrt(cov22_M2)

F_M1=(syy-s_hatyhaty_M1)*n
F_M2=(syy-s_hatyhaty_M2)*n

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
print "\n\nplotting hotel_Ftest.eps"
set out "hotel_Ftest.eps"
#######################################

set param
set key bottom left
set xlabel "f"
set ylabel "Kumulierte Fisher-F-Verteilung F^{2,n-3}(f)"
xmax=14.
set xrange [0:xmax]
plot[t=0:1]\
  xmax*t, fisherCum(xmax*t,2,n-3) t "" w l ls 12,\
  f_realis1, t t str_realH01 w l ls 17,\
  f_realis3, t t str_realH03 w l ls 15,\
  t*xmax, 0.95 t "F=0.95" w l ls 1


#######################################
print "plotting hotel_f_hatbeta1.eps"
set out "hotel_f_hatbeta1.eps"
#######################################

set title "Konfidenzintervall ({/Symbol b}_1) f\374r {/Symbol a}=0.05"
tmin=-1.5*t_nm3_0975
tmax=2.1*t_nm3_0975
sigbeta1=sqrt(cov11)
sigbeta2=sqrt(cov22)
beta1min=beta1+tmin*sigbeta1
beta1max=beta1+tmax*sigbeta1
beta2min=beta2+tmin*sigbeta2
beta2max=beta2+tmax*sigbeta2

studentNoNorm(x,mux,hatsigx,nu)=1./hatsigx*student( (x-mux)/hatsigx,nu)

set xlabel "{/Symbol b}_1 - Sch\344tzer"
set xrange [beta1min:beta1max]
set ylabel "Dichtefunktion f"
set param
set key center right



plot [t=tmin:tmax]\
  beta1+t*sigbeta1,studentNoNorm(beta1+t*sigbeta1,beta1,sigbeta1,n-3)\
    t "Dichte" w l ls 12,\
  beta1-dbeta1 + (t-tmin)/(tmax-tmin)*2*dbeta1,0 t "KI" w l ls 29


#######################################
print "plotting hotel_F_hatbeta1.eps"
set out "hotel_F_hatbeta1.eps"
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
print "plotting hotel_f_hatbeta2.eps"
set out "hotel_f_hatbeta2.eps"
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
print "plotting hotel_F_hatbeta2.eps"
set out "hotel_F_hatbeta2.eps"
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
print "plotting hotel_scatterplot_x1x2.eps"
set out "hotel_scatterplot_x1x2.eps"
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
print "plotting hotel_scatterplot_x1y.eps"
set out "hotel_scatterplot_x1y.eps"
#######################################

set xlabel "Exogene Variable x_1"
set ylabel "Endogene Variable y"
unset key
set key at screen 0.62,0.86



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
  x11,haty(x11,x21) t "Mehrfachregr. (x_1,x_2)" w p ls 7,\
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
print "plotting hotel_scatterplot_x2y.eps"
set out "hotel_scatterplot_x2y.eps"
#######################################



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
  x21,haty(x11,x21) t "Mehrfachregr. (x_1,x_2)" w p ls 7,\
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
print "plotting hotel_scatterplot_x2ycondx1.eps"
set out "hotel_scatterplot_x2ycondx1.eps"
#######################################
min1=17
max1=42
min2=30
max2=60
min3=55
max3=100
min4=80
max4=x2datamax
unset key
set key at screen 0.6,0.84

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
    t "Einfachregression (x_2)" w l ls 11,\
  min1+t*(max1-min1),haty(1, min1+t*(max1-min1) ) t "ydach(x_1=1 Stern, x_2)" w l ls 12,\
  min2+t*(max2-min2),haty(2,min2+t*(max2-min2)) t "ydach (x_1=2 Sterne, x_2)" w l ls 13,\
  min3+t*(max3-min3),haty(3,min3+t*(max3-min3)) t "ydach (x_1=3 Sterne, x_2)" w l ls 14,\
  min4+t*(max4-min4),haty(4,min4+t*(max4-min4)) t "ydach (x_1=4 Sterne, x_2)" w l ls 15

#######################################
print "plotting hotel_scatterplot_x2ycondx1_notCal1.eps"
set out "hotel_scatterplot_x2ycondx1_notCal1.eps"
#######################################

set key at screen 0.6,0.79


beta1_notCal1=beta1+dbeta1
beta2_notCal1=beta2+dbeta2
haty_notCal1(x1,x2)=beta0+beta1_notCal1*x1+beta2_notCal1*x2
set title "{/Symbol b}_1 und {/Symbol b}_2 um {/Symbol Db_1} bzw. {/Symbol Db_2} verschoben"

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
      t "ydach(x_1=1 Stern, x_2)" w l ls 12,\
  min2+t*(max2-min2),haty_notCal1(2,min2+t*(max2-min2))\
      t "ydach (x_1=2 Sterne, x_2)" w l ls 13,\
  min3+t*(max3-min3),haty_notCal1(3,min3+t*(max3-min3))\
      t "ydach (x_1=3 Sterne, x_2)" w l ls 14,\
  min4+t*(max4-min4),haty_notCal1(4,min4+t*(max4-min4))\
      t "ydach (x_1=4 Sterne, x_2)" w l ls 15

#######################################
print "plotting hotel_scatterplot_x2ycondx1_notCal2.eps"
set out "hotel_scatterplot_x2ycondx1_notCal2.eps"
#######################################

set title "{/Symbol b}_1 und {/Symbol b}_2 um {/Symbol Db_1} bzw. - {/Symbol Db_2} verschoben"
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
      t "ydach(x_1=1 Stern, x_2)" w l ls 12,\
  min2+t*(max2-min2),haty_notCal2(2,min2+t*(max2-min2))\
      t "ydach (x_1=2 Sterne, x_2)" w l ls 13,\
  min3+t*(max3-min3),haty_notCal2(3,min3+t*(max3-min3))\
      t "ydach (x_1=3 Sterne, x_2)" w l ls 14,\
  min4+t*(max4-min4),haty_notCal2(4,min4+t*(max4-min4))\
      t "ydach (x_1=4 Sterne, x_2)" w l ls 15

#######################################
print "plotting hotel_scatterplot_x2ycondx1_notCal3.eps"
set out "hotel_scatterplot_x2ycondx1_notCal3.eps"
#######################################

set title "{/Symbol b}_1 und {/Symbol b}_2 um - {/Symbol Db_1} bzw. - {/Symbol Db_2} verschoben"
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
      t "ydach(x_1=1 Stern, x_2)" w l ls 12,\
  min2+t*(max2-min2),haty_notCal3(2,min2+t*(max2-min2))\
      t "ydach (x_1=2 Sterne, x_2)" w l ls 13,\
  min3+t*(max3-min3),haty_notCal3(3,min3+t*(max3-min3))\
      t "ydach (x_1=3 Sterne, x_2)" w l ls 14,\
  min4+t*(max4-min4),haty_notCal3(4,min4+t*(max4-min4))\
      t "ydach (x_1=4 Sterne, x_2)" w l ls 15

#######################################
print "plotting hotel_scatterplot_x2ycondx1_notCal4.eps"
set out "hotel_scatterplot_x2ycondx1_notCal4.eps"
#######################################

set title "{/Symbol b}_1 und {/Symbol b}_2 um - {/Symbol Db_1} bzw. +{/Symbol Db_2} verschoben"

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
      t "ydach(x_1=1 Stern, x_2)" w l ls 12,\
  min2+t*(max2-min2),haty_notCal4(2,min2+t*(max2-min2))\
      t "ydach (x_1=2 Sterne, x_2)" w l ls 13,\
  min3+t*(max3-min3),haty_notCal4(3,min3+t*(max3-min3))\
      t "ydach (x_1=3 Sterne, x_2)" w l ls 14,\
  min4+t*(max4-min4),haty_notCal4(4,min4+t*(max4-min4))\
      t "ydach (x_1=4 Sterne, x_2)" w l ls 15
set notitle




#######################################
print "plotting hotel_f2_hatbeta1_hatbeta2.eps"
set out "hotel_f2_hatbeta1_hatbeta2.eps"
#######################################

set term post eps enhanced color solid "Helvetica" 18
set multiplot
set param
set view map
set contour surface
set cntrparam bspline
set cntrparam levels 10 
set isosamples 30,30
set palette defined ( 0 "white", 5 "yellow", 50 "orange",\
  80 "#FF6666",  99 "#FF44BB", 100 "white") 

xmin=20
xminGamH0=0.
xmax=beta1+3*sigbeta1

ymin=beta2-4*sigbeta2
ymax=0

set xlabel "{/Symbol b}_1 - Sch\344tzer"
set xrange [xmin:xmax]

set ylabel "{/Symbol b}_2 - Sch\344tzer"
set yrange [ymin:ymax]


# (1) Main density plot w/ contour lines

set pm3d
unset clabel
unset key
splot[x=xmin:xmax][y=ymin:ymax]\
  x,y,student2d((x-beta1)/sigbeta1,(y-beta2)/sigbeta2,r_12,n-3)\
   t "f_2(hat({/Symbol b})_1, hat({/Symbol b})_2)" w l ls 98

# (2) Ellipsoid-shaped confidence region

unset pm3d; unset surface # here unset surface necessary!
set cntrparam levels discrete 0.05  # eine Contour!
#set key
set key at screen -0.12,0.96
# BUG doppelt; workaround left label left unvisible

splot[x=xmin:xmax][y=ymin:ymax]\
  x,y,student2d((x-beta1)/sigbeta1,(y-beta2)/sigbeta2,r_12,n-3)\
   t "Konfidenzregion F-Test" w l ls 17

#(3) one-parameter confidence limits

set surface; # here set surface necessary!
#bug falls y-Wert > 0.96
set key at screen 0.92,0.96

#BUGS in TUD gnuplot but not at home (aber Ubuntu 12 LTS!)

# BUG 3d necessary because of scaling bug otherwise!!!
#     lines produce unreproducible/illogical bugs => points OK (strange)
# BUG surface necessary but produces huge files for reasonable sampling=>png
# set isosamples 100  # because of points ("w l" <-> "w l" etc)

splot[t=0:1]\
 beta1-dbeta1, ymin+t*(ymax-ymin),10\
      t "Konfidenzintervalle T-Test {/Symbol b_1, b_2}" w l ls 16,\
 beta1+dbeta1, ymin+t*(ymax-ymin),0 t "" w l ls 16,\
 xmin+t*(xmax-xmin), beta2-dbeta2,0 t "" w l ls 16,\
 xmin+t*(xmax-xmin), beta2+dbeta2,0 t "" w l ls 16,\
 xminGamH0+t*(xmax-xminGamH0), erwGam1-c1*(xminGamH0+t*(xmax-xminGamH0))-dGam1,0\
   t "Konfidenzintervalle T-Test {/Symbol g}_1, {/Symbol g}_2" w l ls 12,\
 xminGamH0+t*(xmax-xminGamH0), erwGam1-c1*(xminGamH0+t*(xmax-xminGamH0))+dGam1,0\
   t "" w l ls 12,\
 xminGamH0+t*(xmax-xminGamH0), erwGam2-c2*(xminGamH0+t*(xmax-xminGamH0))-dGam2,0\
   t "" w l ls 12,\
 xminGamH0+t*(xmax-xminGamH0), erwGam2-c2*(xminGamH0+t*(xmax-xminGamH0))+dGam2,0\
   t "" w l ls 12,\
 20+30*1.2*t,-2./3. -1.2*t,0 t "Grenze Test {/Symbol b_1 + 30 b_2}<0" w l ls 21


#(4) verbundene H0 (always w p, not w l !!)

set key at screen 0.40,0.92

splot[t=0:1]\
  beta101,beta201,0.01 t "Verbundene Nullhypothese H_{01}" w p ls 1,\
  beta103,beta203,0.01 t "Verbundene Nullhypothese H_{03}" w p ls 11


set nomultiplot
set notitle
set isosamples 30




#######################################
print "plotting hotel_f2_hatbeta1_hatbeta2_simple.eps"
set out "hotel_f2_hatbeta1_hatbeta2_simple.eps"
#######################################

# gauss2d, not student2d as the "not simple" plot!

set nokey
set title "f(hat({/Symbol b})_1, hat({/Symbol b})_2)"
unset cntrparam
set cntrparam levels 10 
set pm3d
xmin=beta1-2.2*sigbeta1
xmax=beta1+2.2*sigbeta1
ymin=beta2-2.2*sigbeta2
ymax=beta2+2.2*sigbeta2

set auto x
set auto y
splot[x=xmin:xmax][y=ymin:ymax]\
  x,y,gauss2d(x-beta1,y-beta2,sigbeta1,sigbeta2,r_12)\
  t "f_2(hat({/Symbol b})_1, hat({/Symbol b})_2)" w l ls 98
unset title

#######################################
print "plotting hotel_fisherCum.eps"
set out "hotel_fisherCum.eps"
#######################################
unset colorbox
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
 f, fisherCum(f,2,n-3) t "" w l ls 12

#######################################
print "plotting hotel_fisherCumSolved.eps"
set out "hotel_fisherCumSolved.eps"
#######################################
set size 1,0.7
set key at screen 0.6,0.2
plot[f=fmin:fmax]\
 f, fisherCum(f,2,n-3) t "" w l ls 12,\
 f_realis1, 0.55+0.45*(f-fmin)/(fmax-fmin) t str_realH01 w l ls 17



#######################################
print "plotting hotel_studentCum.eps"
set out "hotel_studentCum.eps"
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
print "plotting hotel_studentCumSolved.eps"
set out "hotel_studentCumSolved.eps"
#######################################
set key at screen 0.8,0.2
plot[t=tmin:tmax]\
 t, studentCum(t,n-3) t "" w l ls 12,\
 t1, 0.8+0.2*(t-tmin)/(tmax-tmin) t "Realisierung symm. Test" w l ls 15,\
 t2asym, 0.8+0.2*(t-tmin)/(tmax-tmin) t "Realisierung asymm. Test" w l ls 17



#######################################
print "plotting hotel_scatter3d_1.eps"
set out "hotel_scatter3d_1.eps"
#######################################
set term post eps enhanced color solid "Helvetica" 18
set nogrid
unset colorbox

set param
set multiplot
set contour surface
set cntrparam bspline
unset cntrparam; set cntrparam levels 20 
unset clabel
set isosamples 15,15




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

set key at screen 0.82,0.67
splot[x1=x1datamin:x1datamax][x2=x2datamin:x2datamax]\
  x1,x2,hatymod(x1,x2)\
   t "{/Symbol b}_0+{/Symbol b}_1 x_1+{/Symbol b}_2 x_2" w l ls 98


set key at screen 0.82,0.62
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
set key at screen 0.82,0.57
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
print "plotting hotel_scatter3d_2.eps"
set out "hotel_scatter3d_2.eps"
#######################################

set multiplot

set view 31,65
set xrange [x1datamin:x1datamax] noreverse


set key at screen 0.82,0.67
splot[x1=x1datamin:x1datamax][x2=x2datamin:x2datamax]\
  x1,x2,hatymod2(x1,x2)\
   t "{/Symbol b}_0+{/Symbol b}_1 x_1+{/Symbol b}_2 x_2" w l ls 98


set key at screen 0.82,0.62
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
set key at screen 0.82,0.57
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
unset key

#######################################
print "plotting hotel_SSE.eps"
set out "hotel_SSE.eps"
#######################################

set palette defined ( 0 "#99003", 5  "red", 20 "orange",\
  50 "yellow",  100 "white") 


unset contour
set pm3d; set pm3d map  
set contour surface

set colorbox

set isosample 50,50
set auto
set xlabel "{/Symbol b}_1"
set ylabel "{/Symbol b}_2"
set title "SSE(hat({/Symbol b}), {/Symbol b}_1, {/Symbol b}_2)"

#set xrange[x1datamin:x1datamax]
#set yrange[x2datamin:x2datamax]
set zrange[SSEmin:100*SSEmin]
set cbrange [SSEmin:100*SSEmin]
set cntrparam levels incr 1.01*SSEmin, 2*SSEmin, 10*SSEmin

splot[b1=beta1-4*dbeta1:beta1+4*dbeta1][b2=beta2-4*dbeta2:beta2+4*dbeta2]\
 b1,b2,SSE(beta0,b1,b2) w l ls 98
