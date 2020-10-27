
# Allgemein

gauss(mu,sigma,x) = exp(-(x-mu)**2/(2* sigma**2)) / (sigma * sqrt(2*pi))
phi(x)            =0.5*(erf(xlimited(x/sqrt(2)))+1)
                  # integral of the Standard Normal Distribution

max(x,y)=(x>y) ? x : y
min(x,y)=(x>y) ? y : x

n=20.

ybar=4.
y_sigRes=3.
set title "n=20;  Standard Deviation of the Residual Errors:  {/Symbol s}_{/Symbol e}=3"
xbar=3.
sigx=3.

b=0.5
a=ybar-b*xbar

regr_erw(x)=a+b*x
regr_sig2(x)=y_sigRes**2/n*(1.+((x-xbar)/sigx)**2)
regr_sig(x)=sqrt(regr_sig2(x))
fdens_regr(x,y)=gauss(regr_erw(x),regr_sig(x),y)

#fdens_regr(x,y)=(abs(y-regr_erw(x))>0.2*regr_sig(xbar))\
# ? gauss(regr_erw(x),regr_sig(x),y)\
# : 1*gauss(regr_erw(xbar),regr_sig(xbar),ybar)


#################################
# Aufloesung analytischer Plots
#################################

set isosample 100,100
#set isosample 30,30
#set samples 0,50
#set isosamples 20


#################################
# Definition of the Contourlinien  (falls set contour aktiv, s.u.)
#################################


set cntrparam bspline 

set cntrparam levels 12                          # n equidistant levels from min-max
#set cntrparam levels incr -10,1,-2          # min, increment, max
#set cntrparam levels discrete -8,-1.6,0,1 # freely set lines
#set cntrparam levels discrete -8,-7,-6,-5,-4,-3,-2,-1,-0.01,0,0.01,1,2

unset clabel  # dann lauter gleiche Kontourlinien; 
                     # Farbe and Typ "w l ls" beim splot-Kommando




#################################
# Color-Coding Oberflaeche
#################################


#set palette color            # is default
#set palette model RGB  # is default

# set palette defined (z1, farbe1, ...), wobei zi immer automat. auf
# min/max normiert werden
# (unabhaengig sowohl vom Bereich of the zi als auch of the tats. z-Valuee!!)
# tatsaechl. min/max. of the Farben unabhaengig von zrange 
# "set cbrange" steuerbar

valGreen=0.1 # geht auch Variablen and Float's 

set palette defined ( 0 "white", 0.8"yellow",  1.6 "orange",\
      2.4 "red", 4 "#8800ff") 

#set cbrange [-8:2]  #[min:max] of color-coded  range (default zrange)


#################################
# Aussehen of the Plots: 
# Farbige Oberflaechen, Contourlinien, Gitternetz ...
# Reihenfolge of the sets and unsets egal!
#################################


####### (1) Farbcodierte Oberflaeche+optionales Gitternetz ######

# Achtung!! "pm3d map" setzt "set view" um danach 3d => neues "set view" noetig!
# Achtung!! "pm3d map" laesst setting "hidden3d" unveraendert (auch nach unset!)
# => "set pm3d", DANN "set pm3d map" noetig fuer Wegbringen of the Linien!

set style line 99 lt 7 lw 1        # Farbe+Dicke fuer optionales Gitternetz
unset pm3d                            # no color coding
#set pm3d                                # color coded 3d surface
#set pm3d  hidden3d 99         # color coded 3d surface with grid, ls 99
#set pm3d  map                       # color coded xy surface
#set pm3d map hidden3d 99  # color coded xy surface with grid


####### (2) Contourlinien (Linienstil immer bei  splot definieren!) ######

unset contour              # no contour lines
set contour surface      # Aktiviert Kontourlinien auf 3D-Flaeche (egal bei pm3d map)
#set contour base          # Aktiviert Kontourlinien auf xy-Ebene (egal bei pm3d map)
#set contour both           # Aktiviert Kontourlinien auf beidem


####### (3) Gitternetz (Linienstil: pm3d on => dort def.; off=> splot) ######


#set surface         # pm3d off: Erzeugt ueberhaupt Gitternetz
                            # pm3d on: Realisiert Verdeckte Flaechen and KS-Linien
unset surface     # pm3d off: kein Gitternetz gemacht
                            # pm3d on: Schnelles Plotten mit/ohne Netz
                            # (je nach pm3d Einstellung) Artefakten



#################################
# eigentliches Plotten
#################################

# Allgemeines

set noparam
set term post eps color solid enhanced "Helvetica" 36
set nogrid  #ansonsten 2D-Gitter auf xy-Ebene
set size 2,2
set nokey

set grid

set xlabel "Independent Variable x"
set ylabel "Dependent Variable y"
set zlabel "f(y|x)" offset 0,1
xmin=xbar-3*sigx
xmax=xbar+3*sigx
set xrange [xmin:xmax]
set yrange [a+b*xmin-regr_sig(xmin): a+b*xmax+regr_sig(xmax)]
#set xrange [5:70] reverse
#set yrange [10:35] 
#set xtics 10,10


#set auto
set view 20,250


#################################

set out "regression3d_eng.eps"
print "plotting regression3d_eng.eps ..."


set style line 99 lt 6 lw 2        # Farbe+Dicke fuer optionales Gitternetz+Contourlinien


set pm3d; set pm3d map 
splot\
 fdens_regr(x,y) w l ls 99


#################################

quit

set out "color2dContour2d_eng.eps"
print "plotting color2dContour2d_eng.eps ..."

unset surface
set contour surface
set pm3d
set pm3d; set pm3d map 
set yrange [10:35] 
set grid
splot aIDMdv0bmax(y,x) w l lt 6 lw 4 
# Rueckgaengigmachen of the spezif. 2D-Einstellungen
set nogrid
set yrange [10:35] reverse
set view 10,350

#################################

set out "color3dContour3d_eng.eps"
print "plotting color3dContour3d_eng.eps ..."
unset surface
set pm3d  
set contour surface
splot aIDMdv0bmax(y,x) w l lt 6 lw 4 

#################################

set out "color3dContour3dGrid3d_eng.eps"
print "plotting color3dContour3dGrid3d_eng.eps ..."
unset surface
set style line 99 lt 7 lw 1        # Farbe+Dicke fuer optionales Gitternetz
set pm3d  hidden3d 99 
set contour surface
splot aIDMdv0bmax(y,x) w l lt 6 lw 4 

#################################

set out           "color3dContour3dGrid3dHidden3d_eng.eps"
print "plotting color3dContour3dGrid3dHidden3d_eng.eps ..."
set surface
set style line 99 lt 7 lw 1        # Farbe+Dicke fuer optionales Gitternetz
set pm3d  hidden3d 99 
set contour surface
splot aIDMdv0bmax(y,x) w l lt 6 lw 4 

#################################

set out           "color3dContour2d3dGrid3dHidden3d_eng.eps"
print "plotting color3dContour2d3dGrid3dHidden3d_eng.eps ..."
set surface
set pm3d  hidden3d 99 
set contour both
splot aIDMdv0bmax(y,x) w l lt 6 lw 4 

#################################

set out "onlyContour3dGrid3d_eng.eps"
print "plotting onlyContour3dGrid3d_eng.eps ..."
unset pm3d
set surface
set contour surface
splot aIDMdv0bmax(y,x) w l lt 2 lw 4 

quit


#################################
# schnitt a(dv,s)
#################################

set out "accIDM.dvs_eng.eps"
print "plotting accIDM.dvs_eng.eps ..."
set pm3d  
set contour surface

set xlabel "s (m)" offset -1
set ylabel "{/Symbol D}v (m/s)"
set zlabel "a (m/s^2)" offset 0,1
set xrange [0:70] reverse
set yrange [-2:10] reverse
set zrange [-9:3]  #reverse
#set ztics -8,4
set ztics ("-8" -8, "0" 0)
set auto
set view 10,60

splot aIDMv20bmax(y,x) w l ls 99



set view map
set nozlabel
unset surface
set xtics 0,10
set xrange [0:70]
set yrange [-2:10]
set zrange [-9:3]
set out "accIDM.dvs_cont_eng.eps"
print "plotting accIDM.dvs_cont_eng.eps ..."

splot aIDMv20bmax(y,x) w l ls 99

