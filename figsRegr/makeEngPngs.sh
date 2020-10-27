#!/bin/bash

# converts all relevant image files from engl figure files to png
# including intermediate eps for vector reference

myconvert='convert -flatten -units PixelsPerInch -density 150 -trim -matte -type PaletteMatte'

for f in `ls *_eng.fig`; do
    b=`basename $f .fig`
    fig2dev -m 3 -L png $b.fig $b.png  #m3=magnification 3 times OK
    #fig2dev -L eps $b.fig $b.eps
    echo "made $b.png"
done

# gnuplot output

for f in `ls *_eng.eps`; do
    b=`basename $f .eps`
    convert -flatten $b.eps $b.png
done


