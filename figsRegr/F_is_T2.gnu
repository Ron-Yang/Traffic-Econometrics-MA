# copied from regr_2exogVar_OEV.gnu

x11=1;          x21=30;         y1=280;
x12=1.1;        x22=22;         y2=230;
x13=1.5;        x23=25;         y3=200;
x14=1.6;        x24=21;         y4=140;
x15=2.0;        x25=29;         y5=190;
x16=2.2;        x26=37;         y6=200;
x17=2.4;        x27=35;         y7=200;
x18=2.5;        x28=32;         y8=180;
x19=3.0;        x29=41;         y9=210;
x10=3.2;        x20=28;         y0=110;

n=10.
J=2.


x0bar=1;
x1bar=(x11+x12+x13+x14+x15+x16+x17+x18+x19+x10)/n
x2bar=(x21+x22+x23+x24+x25+x26+x27+x28+x29+x20)/n
ybar=(y1+y2+y3+y4+y5+y6+y7+y8+y9+y0)/n

#################################################
# 3 x 3 - Matrix (\m{X}\tr\cdot\m{X})_{jk}=sumjk
#################################################

sum00=n
sum11=x11**2+x12**2+x13**2+x14**2+x15**2+x16**2+x17**2+x18**2+x19**2+x10**2
sum22=x21**2+x22**2+x23**2+x24**2+x25**2+x26**2+x27**2+x28**2+x29**2+x20**2

sum01=n*x1bar; sum10=sum01;
sum02=n*x2bar; sum20=sum02;
sum12=x11*x21+x12*x22+x13*x23+x14*x24+x15*x25+x16*x26+x17*x27+x18*x28+x19*x29+x10*x20
sum21=sum12;

#################################################
# 3-Vektor (\m{X}\tr\cdot \vec{y})_j=sumjy
#################################################

sum0y=n*ybar
sum1y=x11*y1+x12*y2+x13*y3+x14*y4+x15*y5+x16*y6+x17*y7+x18*y8+x19*y9+x10*y0
sum2y=x21*y1+x22*y2+x23*y3+x24*y4+x25*y5+x26*y6+x27*y7+x28*y8+x29*y9+x20*y0

#################################################
# endogenous var square sum 
#################################################

sumyy=y1**2+y2**2+y3**2+y4**2+y5**2+y6**2+y7**2+y8**2+y9**2+y0**2

#################################################
# exog-endogen covariances
#################################################

s0y=sum0y/n - x0bar*ybar  # =0 if xi0=1
s1y=sum1y/n - x1bar*ybar
s2y=sum2y/n - x2bar*ybar
syy=sumyy/n - ybar**2   # Descr total variance


#################################################
# Matrixinverse (X'X)^{-1}
#################################################

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
# parameter full model
#########################################################

beta0=inv00*sum0y+inv01*sum1y+inv02*sum2y
beta1=inv10*sum0y+inv11*sum1y+inv12*sum2y
beta2=inv20*sum0y+inv21*sum1y+inv22*sum2y

#########################################################
# SSE
#########################################################

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
 + SE(x10,x20,y0,beta0,beta1,beta2)

#########################################################
# inductive residual variance in two ways
#########################################################

# (i) with descriptive covariances (s0y=0 if xi0=1)

s_hatyhaty    = beta0*s0y+beta1*s1y+beta2*s2y     # descr explained var
s_epseps      = syy-s_hatyhaty                    # descr. residual var
hatsig_epseps = n/(n-1-J) * s_epseps              # induct. residual var

# (ii) with SSE

hatsig_epseps_SSE=SSE(beta0,beta1,beta2)/(n-1-J)

print "full model:"
print "  hatsig_epseps=",hatsig_epseps
print "  hatsig_epseps_SSE=",hatsig_epseps_SSE


#########################################################
# estimation variance covariance matrix
#########################################################

hatV00=hatsig_epseps*inv00
hatV11=hatsig_epseps*inv11
hatV22=hatsig_epseps*inv22
hatV01=hatsig_epseps*inv01; hatV10=hatV01
hatV02=hatsig_epseps*inv02; hatV20=hatV02
hatV12=hatsig_epseps*inv12; hatV21=hatV12




#########################################################
# restricted model beta2=0
#########################################################


detXXr=sum00*sum11-sum01**2
inv00r=sum11/detXXr;
inv11r=sum00/detXXr;
inv01r=-sum01/detXXr;
inv10r=inv01r

beta0r=inv00r*sum0y+inv01r*sum1y
beta1r=inv10r*sum0y+inv11r*sum1y


#########################################################
# t-test and F-test for beta2=0
#########################################################

tdata_tTest=beta2/sqrt(hatV22)
tdata_FTest=(n-J-1)*(SSE(beta0r,beta1r,0)/SSE(beta0,beta1,beta2)-1)

print "SSE(beta0r,beta1r,0)=",SSE(beta0r,beta1r,0)
print "SSE(beta0,beta1,beta2)=",SSE(beta0,beta1,beta2)


print "H0: beta2=0"
print "tdata_FTest=",tdata_FTest
print "tdata_tTest**2=",tdata_tTest**2





