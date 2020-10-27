
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

  # integral of the Standard Normal Distribution
phi(x)            =norm(x)

                
 # Quantil of the Standard Normal Distribution
phiQuantil(q)=invnorm(q)

                


max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x

set term post eps enhanced color solid "Helvetica" 18
#set term post eps enhanced color dashed "Helvetica" 20

####################################################
####################################################
# Data 
# x=kilometer
# f=Anteil PT/Car bei Summe of the RC-Umfragen WS 2014 and 2015
# n=Zahl of the Personen in Gruppe

# y* = ln(f1/(1-f1)), #f1=Anteil Alternative 1 (hier PT/Car)
ystar(f1)=log(f1/(1.-f1))

# Umkehrfunktion=logist. Distribution
f1(ystar)=exp(ystar)/(exp(ystar)+1.)


n1=10.; x1=0.5; f1=0.6;   ys1=ystar(f1)
n2=8.;  x2=1.5; f2=0.375;  ys2=ystar(f2)
n3=11.; x3=2.5; f3=9./11.; ys3=ystar(f3)
n4=16.; x4=4.0; f4=15./16.; ys4=ystar(f4)
n5=12.; x5=7.5; f5=0.999999;  ys5=ystar(f5)  # n=5 als Demo fuer Ungeeignetheit

#######################################
print "plotting regr_logistic_WS1516cum_eng.eps"
set out "regr_logistic_WS1516cum_eng.eps"
#######################################

n=n1+n2+n3+n4

xbar=(n1*x1+n2*x2+n3*x3+n4*x4)/n
ysbar=(n1*ys1+n2*ys2+n3*ys3+n4*ys4)/n

sxx=1/n*( n1*(x1-xbar)**2+n2*(x2-xbar)**2+n3*(x3-xbar)**2+n4*(x4-xbar)**2)
sxys=1/n*( n1*(x1-xbar)*(ys1-ysbar)+n2*(x2-xbar)*(ys2-ysbar)\
     +n3*(x3-xbar)*(ys3-ysbar)+n4*(x4-xbar)*(ys4-ysbar))
print "sxx=",sxx
beta1=sxys/sxx
beta0=ysbar-beta1*xbar
haty(x)=beta0+beta1*x

print "  beta0=",beta0
print "  beta1=",beta1


set xlabel "Distance x_1 [km]"
set ylabel "y^*=ln(f/(1-f))"
set key right bottom

set notitle

plot[x=x1-0.5:x4+1]\
  x1,ys1 t "Data" w p ls 1,\
  x2,ys2 t "" w p ls 1,\
  x3,ys3 t "" w p ls 1,\
  x4,ys4 t "" w p ls 1,\
  x, haty(x) t "Logistic Regression" w l ls 2

#######################################
print "plotting regr_logistic_WS1516cum_f_eng.eps"
set out "regr_logistic_WS1516cum_f_eng.eps"
#######################################

set ylabel "Modal Split PT/Car together [%]"
plot[x=x1-0.5:x4+1]\
  x1,f1(ys1) t "Data" w p ls 1,\
  x2,f1(ys2) t "" w p ls 1,\
  x3,f1(ys3) t "" w p ls 1,\
  x4,f1(ys4) t "" w p ls 1,\
  x, f1(haty(x)) t "Logistic Regression" w l ls 2



#######################################
print "plotting regr_logistic_WS1516cum_alt_eng.eps"
set out "regr_logistic_WS1516cum_alt_eng.eps"
#######################################

n=n1+n2+n3+n4+n5

xbar=(n1*x1+n2*x2+n3*x3+n4*x4+n5*x5)/n
ysbar=(n1*ys1+n2*ys2+n3*ys3+n4*ys4+n5*ys5)/n

sxx=1/n*( n1*(x1-xbar)**2+n2*(x2-xbar)**2+n3*(x3-xbar)**2+n4*(x4-xbar)**2+n5*(x5-xbar)**2)
sxys=1/n*( n1*(x1-xbar)*(ys1-ysbar)+n2*(x2-xbar)*(ys2-ysbar)\
     +n3*(x3-xbar)*(ys3-ysbar)+n4*(x4-xbar)*(ys4-ysbar)+n5*(x5-xbar)*(ys5-ysbar))
print "sxx=",sxx
beta1=sxys/sxx
beta0=ysbar-beta1*xbar
haty(x)=beta0+beta1*x

print "  beta0=",beta0
print "  beta1=",beta1


set xlabel "Distance x_1 [km]"
set ylabel "y^*=ln(f/(1-f))"
set key right bottom

set notitle

plot[x=x1-0.5:x5+1]\
  x1,ys1 t "Data" w p ls 1,\
  x2,ys2 t "" w p ls 1,\
  x3,ys3 t "" w p ls 1,\
  x4,ys4 t "" w p ls 1,\
  x5,ys5 t "" w p ls 1,\
  x, haty(x) t "Logistic Regression" w l ls 2

#######################################
print "plotting regr_logistic_WS1516cum_alt_f_eng.eps"
set out "regr_logistic_WS1516cum_alt_f_eng.eps"
#######################################

set ylabel "Modal Split PT/Car together [%]"
plot[x=x1-0.5:x5+1]\
  x1,f1(ys1) t "Data" w p ls 1,\
  x2,f1(ys2) t "" w p ls 1,\
  x3,f1(ys3) t "" w p ls 1,\
  x4,f1(ys4) t "" w p ls 1,\
  x5,f1(ys5) t "" w p ls 1,\
  x, f1(haty(x)) t "Logistic Regression" w l ls 2



