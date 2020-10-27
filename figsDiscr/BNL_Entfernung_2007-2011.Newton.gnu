


#########################################################
# Minimalbeispiel: 2 Modi, 2 Modellparameter, Studentenbefragung 
# in VL zu Wegeentfernung und Verkehrsmittelwahl 2008
# k=1 => Fuss+Rad
# k=2 => \"OPNV und MIV
#########################################################

#########################################################
# Mittl. Entfernung Entf.-klasse l; 
# BNL_kalib2.gnu unterscheidet sich nur in anderen Klassenmittelwerten!
r1=0.5
r2=1.5   
r3=3.5
r4=7.5
r5=15.

#########################################################
#########################################################
# Abs. Haeufigkeiten y_{ni}

y11=33.; y12=10.;
y21=32.; y22=22.;
y31=29.; y32=59.;
y41=7.;  y42=49.;
y51=0.;  y52=25.;


y1=y11+y12   # Zahl aller Befragten in Entf. Klasse 1 
y2=y21+y22 
y3=y31+y32 
y4=y41+y42 
y5=y51+y52 

# Beobachtete Gesamtkilometerzahl in Verkehrsmodus k=1(Fuss+Rad)
R1=y11*r1+y21*r2+y31*r3+y41*r4+y51*r5
R2=y12*r1+y22*r2+y32*r3+y42*r4+y52*r5
R=R1+R2
print ""
print "Merkmalssummen (Alt 1 relevant):\n"
print "R1=",R1
print "R2=",R2

# Beobachtete Gesamtzahl an Entscheidungen fuer k=1
N1=y11+y21+y31+y41+y51
N2=y12+y22+y32+y42+y52
N=N1+N2
print "N1=",N1
print "N2=",N2
print "N=",N

#V_{ni}=V_i(r_n)

V1(r,beta1,beta2)=beta1*r+beta2
V2(r,beta1,beta2)=0.

# logit probabilities

P1(r,beta1,beta2)=1./(1.+exp(-V1(r,beta1,beta2)))
P2(r,beta1,beta2)=1./(1.+exp(V1(r,beta1,beta2)))

#P_{ki}=Pk(r_i)

P11(beta1,beta2)=P1(r1,beta1,beta2)
P21(beta1,beta2)=P1(r2,beta1,beta2)
P31(beta1,beta2)=P1(r3,beta1,beta2)
P41(beta1,beta2)=P1(r4,beta1,beta2)
P51(beta1,beta2)=P1(r5,beta1,beta2)

P12(beta1,beta2)=P2(r1,beta1,beta2)
P22(beta1,beta2)=P2(r2,beta1,beta2)
P32(beta1,beta2)=P2(r3,beta1,beta2)
P42(beta1,beta2)=P2(r4,beta1,beta2)
P52(beta1,beta2)=P2(r5,beta1,beta2)


# Log-Likelihood

lnL(b1,b2)=\
   y11*log(P11(b1,b2))+y12*log(P12(b1,b2))\
 + y21*log(P21(b1,b2))+y22*log(P22(b1,b2))\
 + y31*log(P31(b1,b2))+y32*log(P32(b1,b2))\
 + y41*log(P41(b1,b2))+y42*log(P42(b1,b2))\
 + y51*log(P51(b1,b2))+y52*log(P52(b1,b2))



# Kalibrierungsbedingung 1: F1=<R(k=1)>_theo - R1=0
# Kalibrierungsbedingung 2: F2=<n(k=1)>_theo - N1=0

F1(b1,b2)=y1*r1*P11(b1,b2)+y2*r2*P21(b1,b2)\
   +y3*r3*P31(b1,b2)+y4*r4*P41(b1,b2)+y5*r5*P51(b1,b2) - R1
F2(b1,b2)=y1*P11(b1,b2)+y2*P21(b1,b2)\
   +y3*P31(b1,b2)+y4*P41(b1,b2)+y5*P51(b1,b2) - N1

R1mod(b1,b2)=F1(b1,b2)+R1
N1mod(b1,b2)=F2(b1,b2)+N1

# PPi=P1i*P2i
PP1(b1,b2)=P11(b1,b2)*P12(b1,b2)
PP2(b1,b2)=P21(b1,b2)*P22(b1,b2)
PP3(b1,b2)=P31(b1,b2)*P32(b1,b2)
PP4(b1,b2)=P41(b1,b2)*P42(b1,b2)
PP5(b1,b2)=P51(b1,b2)*P52(b1,b2)

#Functional matrix d(F_k)/d(beta_k')


J11(b1,b2)=y1*r1**2*PP1(b1,b2)+y2*r2**2*PP2(b1,b2)\
                  +y3*r3**2*PP3(b1,b2)+y4*r4**2*PP4(b1,b2)+y5*r5**2*PP5(b1,b2)
J12(b1,b2)=y1*r1*PP1(b1,b2)+y2*r2*PP2(b1,b2)+y3*r3*PP3(b1,b2)\
       +y4*r4*PP4(b1,b2)+y5*r5*PP5(b1,b2)
J21(b1,b2)=J12(b1,b2)
J22(b1,b2)=y1*PP1(b1,b2)+y2*PP2(b1,b2)+y3*PP3(b1,b2)+y4*PP4(b1,b2)+y5*PP5(b1,b2)

# Inverse of functional matrix

detJ(b1,b2)=J11(b1,b2)*J22(b1,b2)-J12(b1,b2)**2

invJ11(b1,b2)=J22(b1,b2)/detJ(b1,b2)
invJ12(b1,b2)=-J12(b1,b2)/detJ(b1,b2)
invJ21(b1,b2)=-J21(b1,b2)/detJ(b1,b2)
invJ22(b1,b2)=J11(b1,b2)/detJ(b1,b2)

# Newton-update step - invJ . F

dbeta1(b1,b2)= -1* (invJ11(b1,b2)*F1(b1,b2) + invJ12(b1,b2)*F2(b1,b2))
dbeta2(b1,b2)= -1* (invJ21(b1,b2)*F1(b1,b2) + invJ22(b1,b2)*F2(b1,b2))


#############################################
# Perform Newton iteration
#############################################

# Zielfunktion zum lotten u. Auswerten; 2*J =Hesse-Matrix davon
errorSQR(b1,b2)=F1(b1,b2)**2+F2(b1,b2)**2

#############################################
b1=0.
b2=0.
it=0
#############################################


print ""
print "Newton-Iteration"
print "============"
print ""
print "Initialwerte: b1=",b1,"  b2=",b2," errorSQR=",errorSQR(b1,b2)
print "R1=",R1," R1mod(b1,b2)=", R1mod(b1,b2), " grad1=F1(b1,b2)=",F1(b1,b2)
print "N1=",N1," N1mod(b1,b2)=", N1mod(b1,b2), " grad2=F2(b1,b2)=",F2(b1,b2)
it=it+1

print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," errorSQR=",errorSQR(b1,b2)
print "R1=",R1," R1mod(b1,b2)=", R1mod(b1,b2), " grad1=F1(b1,b2)=",F1(b1,b2)
print "N1=",N1," N1mod(b1,b2)=", N1mod(b1,b2), " grad2=F2(b1,b2)=",F2(b1,b2)
it=it+1

print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," errorSQR=",errorSQR(b1,b2)
print "R1=",R1," R1mod(b1,b2)=", R1mod(b1,b2), " grad1=F1(b1,b2)=",F1(b1,b2)
print "N1=",N1," N1mod(b1,b2)=", N1mod(b1,b2), " grad2=F2(b1,b2)=",F2(b1,b2)
it=it+1

print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," errorSQR=",errorSQR(b1,b2)
print "R1=",R1," R1mod(b1,b2)=", R1mod(b1,b2), " grad1=F1(b1,b2)=",F1(b1,b2)
print "N1=",N1," N1mod(b1,b2)=", N1mod(b1,b2), " grad2=F2(b1,b2)=",F2(b1,b2)
it=it+1

print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," errorSQR=",errorSQR(b1,b2)
print "R1=",R1," R1mod(b1,b2)=", R1mod(b1,b2), " grad1=F1(b1,b2)=",F1(b1,b2)
print "N1=",N1," N1mod(b1,b2)=", N1mod(b1,b2), " grad2=F2(b1,b2)=",F2(b1,b2)
it=it+1

print ""
print "iteration ",it,":"
b1new=b1+dbeta1(b1,b2)
b2new=b2+dbeta2(b1,b2)
b1=b1new
b2=b2new
print "b1=",b1,"  b2=",b2," errorSQR=",errorSQR(b1,b2)
print "R1=",R1," R1mod(b1,b2)=", R1mod(b1,b2), " grad1=F1(b1,b2)=",F1(b1,b2)
print "N1=",N1," N1mod(b1,b2)=", N1mod(b1,b2), " grad2=F2(b1,b2)=",F2(b1,b2)
it=it+1

print "Parameter Covariance matrix H^{-1}:\n"
print "V11=invJ11(b1,b2)=",invJ11(b1,b2)
print "V12=invJ12(b1,b2)=",invJ12(b1,b2)
print "V22=invJ22(b1,b2)=",invJ22(b1,b2)
print ""
print "sqrt(V11)=",sqrt(invJ11(b1,b2))
print "sqrt(V22)=",sqrt(invJ22(b1,b2))
print "rbeta12=",invJ12(b1,b2)/(sqrt(invJ11(b1,b2))*sqrt(invJ22(b1,b2)))

print "lnLmax=lnL(b1,b2)=",lnL(b1,b2)
print "lnLmax=lnL(-0.866677,-0.313785)=",lnL(-0.866677,-0.313785)
