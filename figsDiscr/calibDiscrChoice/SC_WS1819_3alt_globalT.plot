#!/bin/bash

gnuplot SC_WS1819_3alt_globalT.gnu
fig2dev -L png SC_WS1819_3alt_globalT.corrMatrix.fig SC_WS1819_3alt_globalT_corrMatrix.png
echo "made SC_WS1819_3alt_globalT_corrMatrix.png"
#gv --orientation PORTRAIT SC_WS1819_3alt_globalT_corrMatrix.png &
pdflatex SC_WS1819_3alt_globalT
echo "produced SC_WS1819_3alt_globalT.pdf"
exit

