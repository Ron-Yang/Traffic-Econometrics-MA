#!/bin/bash

gnuplot RC_skript_4alt_r_eng.gnu
fig2dev -L png RC_skript_4alt_r_eng.corrMatrix.fig RC_skript_4alt_r_eng_corrMatrix.png
#xv RC_skript_4alt_r_eng_fProbKumDist.png &
#xv RC_skript_4alt_r_eng_corrMatrix.png &
pdflatex RC_skript_4alt_r_eng
exit

