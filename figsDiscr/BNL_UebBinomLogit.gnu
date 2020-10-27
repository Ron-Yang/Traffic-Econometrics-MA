



#########################################################
# Minimalbeispiel: 2 Modi, 2 Modellparameter, Studentenbefragung 
# in VL zu Wegeentfernung und Verkehrsmittelwahl 2007-2011
# k=1 => Fuss+Rad
# k=2 => \"OPNV und MIV
#########################################################

set encoding iso_8859_1 

set style line 1 lt 7 lw 6 pt 1 ps 1.2 #schwarz, plus sign
set style line 2 lt 1 lw 8 pt 4 ps 1.2 #rot, open box
set style line 7 lt 3 lw 4 pt 10 ps 1.2 #blau, upside-down open triangle

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x



#########################################################
# Zeiten der beiden Alternativen 1=Fuss/Rad, 2=OEV
T11=20;	 T12=25;
T21=60;	 T22=35;
T31=40;	 T32=30;
T41=80;	 T42=45;
T51=30;	 T52=20;
T61=40;	 T62=40;

#########################################################
# Ad-Hoc Kosten OEV
K12=2;
K22=0;
K32=2;
K42=3.5;
K52=0;
K62=0;

#########################################################
# Abs. Haeufigkeiten y_{ni}

y11=3.;  y12=2.;
y21=0.;  y22=5.;
y31=2.;  y32=3.;
y41=1.;  y42=4.;
y51=0.;  y52=5.;
y61=1.;  y62=4.;


y1=y11+y12   # Zahl aller Befragten in Entf. Klasse 1 
y2=y21+y22 
y3=y31+y32 
y4=y41+y42 
y5=y51+y52 
y6=y61+y62 

#V_{ni}=V_i(r_n)

V1(T,beta1,beta3)=beta1*T+beta3
V2(T,K,beta1,beta2)=beta1*T+beta2*K

#P_{ki}=Pk(r_i)

denom(T1,T2,K2,beta1,beta2,beta3)=exp(V1(T1,beta1,beta3))+exp(V2(T2,K2,beta1,beta2))

P11(b1,b2,b3)=exp(V1(T11,b1,b3))/denom(T11,T12,K12,b1,b2,b3)
P21(b1,b2,b3)=exp(V1(T21,b1,b3))/denom(T21,T22,K22,b1,b2,b3)
P31(b1,b2,b3)=exp(V1(T31,b1,b3))/denom(T31,T32,K32,b1,b2,b3)
P41(b1,b2,b3)=exp(V1(T41,b1,b3))/denom(T41,T42,K42,b1,b2,b3)
P51(b1,b2,b3)=exp(V1(T51,b1,b3))/denom(T51,T52,K52,b1,b2,b3)
P61(b1,b2,b3)=exp(V1(T61,b1,b3))/denom(T61,T62,K62,b1,b2,b3)

P12(b1,b2,b3)=exp(V2(T12,K12,b1,b2))/denom(T11,T12,K12,b1,b2,b3)
P22(b1,b2,b3)=exp(V2(T22,K22,b1,b2))/denom(T21,T22,K22,b1,b2,b3)
P32(b1,b2,b3)=exp(V2(T32,K32,b1,b2))/denom(T31,T32,K32,b1,b2,b3)
P42(b1,b2,b3)=exp(V2(T42,K42,b1,b2))/denom(T41,T42,K42,b1,b2,b3)
P52(b1,b2,b3)=exp(V2(T52,K52,b1,b2))/denom(T51,T52,K52,b1,b2,b3)
P62(b1,b2,b3)=exp(V2(T62,K62,b1,b2))/denom(T61,T62,K62,b1,b2,b3)



logL(b1,b2,b3)=\
     y11*log(P11(b1,b2,b3))+y12*log(P12(b1,b2,b3))\
   + y21*log(P21(b1,b2,b3))+y22*log(P22(b1,b2,b3))\
   + y31*log(P31(b1,b2,b3))+y32*log(P32(b1,b2,b3))\
   + y41*log(P41(b1,b2,b3))+y42*log(P42(b1,b2,b3))\
   + y51*log(P51(b1,b2,b3))+y52*log(P52(b1,b2,b3))\
   + y61*log(P61(b1,b2,b3))+y62*log(P62(b1,b2,b3))


########################################################
print ""
print "Newton-Iteration, Start"
print "======================="
print ""

beta1=0.
beta2=0.
beta3=0.


# Beobachtete und modellierte Gesamtkilometerzahl

T=y11*T11+y21*T21+y31*T31+y41*T41+y51*T51+y61*T61\
 +y12*T12+y22*T22+y32*T32+y42*T42+y52*T52+y62*T62

Tmod(b1,b2,b3)=y1*P11(b1,b2,b3)*T11+y2*P21(b1,b2,b3)*T21+y3*P31(b1,b2,b3)*T31\
              +y4*P41(b1,b2,b3)*T41+y5*P51(b1,b2,b3)*T51+y6*P61(b1,b2,b3)*T61\
              +y1*P12(b1,b2,b3)*T12+y2*P22(b1,b2,b3)*T22+y3*P32(b1,b2,b3)*T32\
              +y4*P42(b1,b2,b3)*T42+y5*P52(b1,b2,b3)*T52+y6*P62(b1,b2,b3)*T62

# Beobachtete und modellierte Gesamtkosten

K=y12*K12+y22*K22+y32*K32+y42*K42+y52*K52+y62*K62

Kmod(b1,b2,b3)=y1*P12(b1,b2,b3)*K12+y2*P22(b1,b2,b3)*K22+y3*P32(b1,b2,b3)*K32\
              +y4*P42(b1,b2,b3)*K42+y5*P52(b1,b2,b3)*K52+y6*P62(b1,b2,b3)*K62

# Beobachtete und modellierte Gesamtentscheidungszahl Alt. 1

N1=y11+y21+y31+y41+y51+y61

N1mod(b1,b2,b3)=y1*P11(b1,b2,b3)+y2*P21(b1,b2,b3)+y3*P31(b1,b2,b3)\
               +y4*P41(b1,b2,b3)+y5*P51(b1,b2,b3)+y6*P61(b1,b2,b3)

print "beta1=",beta1,"  beta2=",beta2,"  beta3=",beta3
print "Gesamtreisezeit: T=",T," Tmod(beta1,beta2,beta3)=",Tmod(beta1,beta2,beta3)
print "Gesamtkosten: K=",K," Kmod(beta1,beta2,beta3)=",Kmod(beta1,beta2,beta3)
print "Gesamtentscheidungszahl Alt. 1: N1=",N1," N1mod(beta1,beta2,beta3)=",N1mod(beta1,beta2,beta3)


print "\nErgebnis Kalibrierung:
beta1=-0.09138
beta2=-1.064
beta3=-1.9043
print "beta1=",beta1,"  beta2=",beta2,"  beta3=",beta3
print "Gesamtreisezeit: T=",T," Tmod(beta1,beta2,beta3)=",Tmod(beta1,beta2,beta3)
print "Gesamtkosten: K=",K," Kmod(beta1,beta2,beta3)=",Kmod(beta1,beta2,beta3)
print "Gesamtentscheidungszahl Alt. 1: N1=",N1," N1mod(beta1,beta2,beta3)=",N1mod(beta1,beta2,beta3)

print "\nWahlhwahrsch. der ersten Person nach Kalibrierung:"
print "P11(beta1,beta2,beta3)=",P11(beta1,beta2,beta3)
print "P12(beta1,beta2,beta3)=",P12(beta1,beta2,beta3)
print "\nPrognostizierte Wahlhaeufigkeiten der ersten Person nach Kalibrierung:"
print "y1*P11(beta1,beta2,beta3)=",y1*P11(beta1,beta2,beta3)
print "y1*P12(beta1,beta2,beta3)=",y1*P12(beta1,beta2,beta3)

quit
