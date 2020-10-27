#!/bin/bash

# makes some standard replacements; further manual edits necessary!

if (($#<1)); 
  then echo "calling sequence: gnuplot_ger2engl.sh <engl gnuplot files>"; 
       exit -1; 
fi

savedir="ger2engl_savedFiles"
if test -d $savedir; then echo "saving unchanged files in $savedir"; 
  else mkdir $savedir;echo "saving unchanged files in $savedir"; 
fi
cp $@ $savedir

for f in "$@"; do
    echo "treating $f ...";
    perl -i -p -e 's/\.eps/_eng\.eps/g' $f
    perl -i -p -e 's/eng_eng_eng_eng/eng/g' $f
    perl -i -p -e 's/eng_eng_eng/eng/g' $f
    perl -i -p -e 's/eng_eng/eng/g' $f
    
    perl -i -p -e 's/Wahrscheinlichkeit/Probability/g'  $f
    perl -i -p -e 's/Freiheitsgrade/DOF/g'  $f
    perl -i -p -e 's/Kumulierte/Cumulated/g'  $f
    perl -i -p -e 's/Dichtefunktion/Density/g'  $f
    perl -i -p -e 's/Dichte/Density/g'  $f
    perl -i -p -e 's/Standardnormalverteilung/Standard Normal Distribution/g'  $f
    perl -i -p -e 's/Verteilungsfunktion/Distribution Function/g'  $f
    perl -i -p -e 's/Verteilungen/Distributions/g'  $f
    perl -i -p -e 's/Verteilung/Distribution/g'  $f
    perl -i -p -e 's/Quantilsfunktion/Quantile Function/g'  $f
    perl -i -p -e 's/Exogene/Exogenous/g'  $f
    perl -i -p -e 's/Endogene/Endogenous/g'  $f
    perl -i -p -e 's/Sch\\344tzer/Estimator/g'  $f
    perl -i -p -e 's/gesch\\344tzt/estimated/g'  $f
    perl -i -p -e 's/Fahrpreis/Price/g'  $f
    perl -i -p -e 's/Geschwindigkeit/Speed/g'  $f
    perl -i -p -e 's/Fahrgastzahlen/Demand/g'  $f
    perl -i -p -e 's/Fahrten/Trips/g'  $f
    perl -i -p -e 's/Jahr/Year/g'  $f
    perl -i -p -e 's/unabhaengige/Independent/g'  $f
    perl -i -p -e 's/abhaengige/Dependent/g'  $f
    perl -i -p -e 's/OEV/PT/g'  $f
    perl -i -p -e 's/PT_/OEV_/g'  $f  # correct side effect file name
    perl -i -p -e 's/MIV/Car/g'  $f
    perl -i -p -e 's/zusammen/together/g'  $f
    perl -i -p -e 's/Entfernung/Distance/g'  $f

    perl -i -p -e 's/G\\374tefunktion/Power Function/g'  $f
    perl -i -p -e 's/Realisierter/Realized/g'  $f
    perl -i -p -e 's/Konfidenzintervalle/Confidence Intervals/g'  $f
    perl -i -p -e 's/Konfidenzintervall/Confidence Interval/g'  $f
    perl -i -p -e 's/Konfidenzregion/Confidence Region/g'  $f
    perl -i -p -e 's/Testvariable/Test Variable/g'  $f
    perl -i -p -e 's/Verbundene Nullhypothese/Compound Null Hypothesis/g'  $f
    perl -i -p -e 's/Einfachregression/Simple Regression/g'  $f
    perl -i -p -e 's/Mehrfachregression/Multiple Regression/g'  $f
    perl -i -p -e 's/Standardabweichung/Standard Deviation/g'  $f
    perl -i -p -e 's/Varianz/Variance/g'  $f
    perl -i -p -e 's/Residualfehler/Residual Error/g'  $f
    perl -i -p -e 's/Logistische/Logistic/g'  $f
    perl -i -p -e 's/Sterne/stars/g'  $f
    perl -i -p -e 's/Grenze/Bounds/g'  $f
    perl -i -p -e 's/Stern/star/g'  $f
    perl -i -p -e 's/dach/ est/g'  $f
    perl -i -p -e 's/Daten/Data/g'  $f
    perl -i -p -e 's/Fehler/Error/g'  $f
    perl -i -p -e 's/Wert/Value/g'  $f
    perl -i -p -e 's/unendlich/infinity/g'  $f
    perl -i -p -e 's/ zu / for /g'  $f
    perl -i -p -e 's/ und / and /g'  $f
    perl -i -p -e 's/ wie / as /g'  $f
    perl -i -p -e 's/ des / of the /g'  $f
    perl -i -p -e 's/ der / of the /g'  $f
    perl -i -p -e 's/ mit / /g'  $f
    perl -i -p -e 's/ FG/ DOF/g'  $f
   
done
 
