

min(x,y)=(x<y) ? x : y


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

h11=33.  # Beobachtete absolute Hauef. Alternative k, Klasse i
h12=32.
h13=29.
h14=7.
h15=0.

h21=10.
h22=22.
h23=59.
h24=49.
h25=25.


n1=h11+h21   # Zahl aller Befragten in Entf. Klasse 1 
n2=h12+h22 
n3=h13+h23 
n4=h14+h24 
n5=h15+h25 

# Beobachtete Gesamtkilometerzahl in Verkehrsmodus k=1(Fuss+Rad)
R1=h11*r1+h12*r2+h13*r3+h14*r4+h15*r5
R2=h21*r1+h22*r2+h23*r3+h24*r4+h25*r5
R=R1+R2
print ""
print "Merkmalssummen (Alt 1 relevant):\n"
print "R1=",R1
print "R2=",R2

# Beobachtete Gesamtzahl an Entscheidungen fuer k=1
N1=h11+h12+h13+h14+h15
N2=h21+h22+h23+h24+h25
n=N1+N2
print "N1=",N1
print "N2=",N2
print "n=",n

#V_{ki}=Vk(r_i)

V1(r,beta1,beta2)=beta1*r+beta2
V2(r,beta1,beta2)=0.

#P_{ki}=Pk(r_i)

P1(r,beta1,beta2)=1./(1.+exp(-V1(r,beta1,beta2)))
P2(r,beta1,beta2)=1./(1.+exp(V1(r,beta1,beta2)))

P11(beta1,beta2)=P1(r1,beta1,beta2)
P12(beta1,beta2)=P1(r2,beta1,beta2)
P13(beta1,beta2)=P1(r3,beta1,beta2)
P14(beta1,beta2)=P1(r4,beta1,beta2)
P15(beta1,beta2)=P1(r5,beta1,beta2)

P21(beta1,beta2)=P2(r1,beta1,beta2)
P22(beta1,beta2)=P2(r2,beta1,beta2)
P23(beta1,beta2)=P2(r3,beta1,beta2)
P24(beta1,beta2)=P2(r4,beta1,beta2)
P25(beta1,beta2)=P2(r5,beta1,beta2)

# Log-Likelihood

lnL(b1,b2)=\
   h11*log(P11(b1,b2))+h21*log(P21(b1,b2))\
 + h12*log(P12(b1,b2))+h22*log(P22(b1,b2))\
 + h13*log(P13(b1,b2))+h23*log(P23(b1,b2))\
 + h14*log(P14(b1,b2))+h24*log(P24(b1,b2))\
 + h15*log(P15(b1,b2))+h25*log(P25(b1,b2))



# Kalibrierungsbedingung 1: F1=<R(k=1)>_theo - R1=0
# Kalibrierungsbedingung 2: F2=<n(k=1)>_theo - N1=0

F1(b1,b2)=n1*r1*P11(b1,b2)+n2*r2*P12(b1,b2)\
   +n3*r3*P13(b1,b2)+n4*r4*P14(b1,b2)+n5*r5*P15(b1,b2) - R1
F2(b1,b2)=n1*P11(b1,b2)+n2*P12(b1,b2)\
   +n3*P13(b1,b2)+n4*P14(b1,b2)+n5*P15(b1,b2) - N1

R1mod(b1,b2)=F1(b1,b2)+R1
N1mod(b1,b2)=F2(b1,b2)+N1

# PPi=P1i*P2i
PP1(b1,b2)=P11(b1,b2)*P21(b1,b2)
PP2(b1,b2)=P12(b1,b2)*P22(b1,b2)
PP3(b1,b2)=P13(b1,b2)*P23(b1,b2)
PP4(b1,b2)=P14(b1,b2)*P24(b1,b2)
PP5(b1,b2)=P15(b1,b2)*P25(b1,b2)

#Functional matrix d(F_k)/d(beta_k')


J11(b1,b2)=n1*r1**2*PP1(b1,b2)+n2*r2**2*PP2(b1,b2)\
                  +n3*r3**2*PP3(b1,b2)+n4*r4**2*PP4(b1,b2)+n5*r5**2*PP5(b1,b2)
J12(b1,b2)=n1*r1*PP1(b1,b2)+n2*r2*PP2(b1,b2)+n3*r3*PP3(b1,b2)\
       +n4*r4*PP4(b1,b2)+n5*r5*PP5(b1,b2)
J21(b1,b2)=J12(b1,b2)
J22(b1,b2)=n1*PP1(b1,b2)+n2*PP2(b1,b2)+n3*PP3(b1,b2)+n4*PP4(b1,b2)+n5*PP5(b1,b2)

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
