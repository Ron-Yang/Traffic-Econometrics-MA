


# siehe ~/info/gnuplot, ~/info/gnuTemplate.gnu

##########################################################
#geordnet nach hue (ps)
# set style line: linetype lt, point type pt

set encoding iso_8859_1 # dann -bäöüßÄÖÜ durch woertl. Eingabe korrekt-A

set style line 1 lt 7 lw 6 pt 1 ps 1.2 #schwarz, plus sign
set style line 2 lt 1 lw 6 pt 4 ps 0.8 #rot, open box
set style line 3 lt 8 lw 4 pt 5 ps 1.2 #blassrot, closed square
set style line 4 lt 6 lw 4 pt 6 ps 1.2 #gelb, open circle
set style line 5 lt 2 lw 4 pt 8  ps 1.5 #gruen, open triangle
set style line 6 lt 5 lw 4 pt 9  ps 1.2 #blasstuerkisblau, closed triangle
set style line 7 lt 3 lw 4 pt 10 ps 1.8 #blau, upside-down open triangle
set style line 8 lt 4 lw 4 pt 11 ps 1.2 #lila, upside-down closed triangle


############### Beispiele fuer Funktionen ####################

gauss(mu,sigma,x) = exp(-(x-mu)**2/(2* sigma**2)) / (sigma * sqrt(2*pi))

xlimited(x)       =(x<-5) ? -5 : ( (x>5) ? 5 : x)
phi(x)            =0.5*(erf(xlimited(x/sqrt(2)))+1) 
                  # integral der Standardnormalverteilung

max(x,y)    = (x>y) ? x : y
min(x,y)    = (x>y) ? y : x


#########################################################
# 4-Alternativen-Beisoiel aus dem Skript
# in VL zu Wegeentfernung und Verkehrsmittelwahl 2008
# k=0 => Fuss
# k=1 => Rad
# k=2 => OEV
# k=3 => MIV
#########################################################

# Mittl. Entfernung Entf.-klasse l; 
# MNL_kalib2.gnu unterscheidet sich nur in anderen Klassenmittelwerten!
r0=1.   
r1=3.5
r2=7.5
r3=15.

r(i)=(i<0.5) ? r0 : (i<1.5) ? r1:(i<2.5) ? r2 : r3

select(k,i)=((i>k-exp(-10))&&(i<k+exp(-10))) ? 1 : -1

delta(k,i)=((i>k-exp(-10))&&(i<k+exp(-10))) ? 1 : 0


# betas aus dem unten zu plottenden Datenfile als Comment am Anfang
beta0=-1.041204
beta1=-0.502038
beta2=-0.267425
beta3=3.199936
beta4=3.249662
beta5=2.715138

V(k,r)=(beta0*r+beta3)*delta(k,0)\
  +(beta1*r+beta4)*delta(k,1) +(beta2*r+beta5)*delta(k,2)

denom(r)=exp(V(0,r)) + exp(V(1,r)) + exp(V(2,r)) + exp(V(3,r))
pMNL(k,r)=exp(V(k,r))/denom(r)

print "denom(4)=",denom(4)
print "pMNL(2,4)=",pMNL(2,4)

#########################################################
# Ergebnis plotten
#########################################################

set term post eps enhanced color solid "Helvetica" 24


set size 1,1
set key at screen 0.55,0.9
set param
set xlabel "Reiseweite (km)"
xmax=20.
set xrange [0:xmax]

set ylabel "Anteilswerte"

set out "MNL_theo1.eps"
print "plotting MNL_theo1.eps"

plot[x=0:xmax]\
 x,pMNL(0,x) w l ls 7,\
 x,pMNL(1,x) w l ls 5,\
 x,pMNL(2,x) w l ls 3,\
 x,pMNL(3,x) w l ls 2

set out "MNL_result1.eps"
print "plotting MNL_result1.eps"

plot[x=0:xmax]\
 x,pMNL(0,x) t "Modell, Fu\337" w l ls 7,\
 "discrChoice_results.dat" u (select($2,0)*r($1)):($5/$7)\
          t "Daten, Fu\337" w p ls 7,\
 x,pMNL(1,x) t "" w l ls 5,\
 "discrChoice_results.dat" u (select($2,1)*r($1)):($5/$7)\
          t "Rad" w p ls 5,\
 x,pMNL(2,x)  t "" w l ls 3,\
 "discrChoice_results.dat" u (select($2,2)*r($1)):($5/$7)\
          t "\326PNV" w p ls 3,\
 x,pMNL(3,x)  t "" w l ls 2,\
 "discrChoice_results.dat" u (select($2,3)*r($1)):($5/$7)\
          t "MIV" w p ls 2

### Achtung! Kum Plotten geht nur "zu Fuss!"

  h00=3.;       h10=5.;	h20=4.; h30=1.;
  h01=1.;       h11=7.;	h21=8.; h31=1.;
  h02=0.;       h12=2.;	h22=8.; h32=3.;
  h03=0.;       h13=0.;	h23=1.; h33=5.;
  n0=h00+h10+h20+h30
  n1=h01+h11+h21+h31
  n2=h02+h12+h22+h32
  n3=h03+h13+h23+h33
  h10kum=h00+h10;   h20kum=h10kum+h20
  h11kum=h01+h11;   h21kum=h11kum+h21
  h12kum=h02+h12;   h22kum=h12kum+h22
  h13kum=h03+h13;   h23kum=h13kum+h23

pMNLkum1(x)=pMNL(0,x) +pMNL(1,x) 
pMNLkum2(x)=pMNLkum1(x) +pMNL(2,x) 

set out "MNL_result1Kum.eps"
print "plotting MNL_result1Kum.eps"

set label 1 "Fu\337" at 0.5,0.05
set label 2 "Rad" at 3,0.2
set label 3 "\326PNV" at 6,0.4
set label 4 "MIV" at 12,0.8
set nokey

plot[x=0:xmax]\
 x,pMNL(0,x) t "Modell, Fu\337" w l ls 7,\
 "discrChoice_results.dat" u (select($2,0)*r($1)):($5/$7) t "" w p ls 7,\
 x,pMNLkum1(x) t "Rad" w l ls 5,\
   r0,h10kum/n0 t "" w p ls 5,\
   r1,h11kum/n1 t "" w p ls 5,\
   r2,h12kum/n2 t "" w p ls 5,\
   r3,h13kum/n3 t "" w p ls 5,\
 x,pMNLkum2(x) t "" w l ls 3,\
   r0,h20kum/n0 t "" w p ls 3,\
   r1,h21kum/n1 t "" w p ls 3,\
   r2,h22kum/n2 t "" w p ls 3,\
   r3,h23kum/n3 t "" w p ls 3

quit
############################################################

