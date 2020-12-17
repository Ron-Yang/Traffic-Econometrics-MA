#!/bin/bash

gnuplot RC_skript_4alt_r_probit.gnu
fig2dev -L eps RC_skript_4alt_r_probit.corrMatrix.fig RC_skript_4alt_r_probit.corrMatrix.eps
gv --orientation PORTRAIT RC_skript_4alt_r_probit_fProb.eps &
gv --orientation PORTRAIT RC_skript_4alt_r_probit.corrMatrix.eps &
exit
#tex2psl RC_skript_4alt_r_probit

