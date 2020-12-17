#!/bin/bash

gnuplot SC_WS1819_3alt_globalT_weather.gnu
fig2dev -L png SC_WS1819_3alt_globalT_weather.corrMatrix.fig SC_WS1819_3alt_globalT_weather_corrMatrix.png
echo "made SC_WS1819_3alt_globalT_weather_corrMatrix.png"
#gv --orientation PORTRAIT SC_WS1819_3alt_globalT_weather_corrMatrix.png &
pdflatex SC_WS1819_3alt_globalT_weather
echo "produced SC_WS1819_3alt_globalT_weather.pdf"
exit

