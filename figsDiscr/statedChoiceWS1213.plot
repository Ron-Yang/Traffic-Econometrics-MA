#!/bin/bash

gnuplot statedChoiceWS1213.gnu
fig2dev -L eps statedChoiceWS1213.corrMatrix.fig statedChoiceWS1213.corrMatrix.eps
echo "made statedChoiceWS1213.corrMatrix.eps"


