gauss(mu,sigma,x) = exp(-(x-mu)**2/(2* sigma**2)) / (sigma * sqrt(2*pi))
gaussNorm(x)=exp(-0.5*x**2) /sqrt(2*pi)
xlimited(x)       =(x<-5) ? -5 : ( (x>5) ? 5 : x)
phi(x)            =0.5*(erf(xlimited(x/sqrt(2)))+1) 
                  # integral der Standardnormalverteilung
pi(x)=norm(x)
lorenz(x)         = 1/(pi*(1+x**2))
expo(lambda,x)    = lambda*exp(-lambda*x)
gleich(a,b,x)     = ((x>=a)&&(x<=b)) ? 1/(b-a) : 0

studnorm(nu)  = gamma(0.5*(nu+1)) / (sqrt(nu*pi)*gamma(0.5*nu))
student(x,nu) =  studnorm(nu) / (1.+x**2/nu)**(0.5*(nu+1))

chinorm(nu)   = 1./(2**(0.5*nu)*gamma(0.5*nu))
chi2(x,nu)    =  chinorm(nu) *exp(-0.5*x) * x**(0.5*nu-1)

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x

eulerKonst=0.577216

#gamma(n+1)=n! for integer n>=0

#erlang(1,x)=standard Poisson distribution

b_erlang(n)=gamma(n+1) /gamma(n)
a_erlang(n)=(b_erlang(n))**(n)/gamma(n)
erlang(n,x)=a_erlang(n)*x**(n-1)*exp(-b_erlang(n)*x)

###########################################################

# Konstanter Binomialkoeffizient-Term (Bsp. I=6 Gruppen)

lnfac(x)=(x<50) ? log(gamma(x+1)) : x*log(x)-x

logB(ni, y1i)=lnfac(ni)-lnfac(y1i)-lnfac(ni-y1i)  

n1=5
n2=5
n3=5
n4=5
n5=5
n6=5
#logL_BinomKoeff=logB(n1, y11)+logB(n2, y12)+logB(n3, y13)\
#  +logB(n4, y14)+logB(n5, y15)+logB(n6, y16)

#print "logL_BinomKoeff=",logL_BinomKoeff



# Konstanter Multinomialkoeffizient-Term (4 Alternativen, Bsp. I=2 Gruppen)

logMul4(ni, y1i, y2i, y3i)= lnfac(ni)-lnfac(y1i)-lnfac(y2i)-lnfac(y3i)-lnfac(ni-y1i-y2i-y3i)  

# logL_MulKoeff4=logMul4(n1, y11,y21,y31)+logMul4(n2, y12,y22,y32)

# print "logL_MulKoeff4=",logL_MulKoeff4
