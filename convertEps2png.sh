#!/bin/bash

# converts all relevant image files to png
# F... annoying bug: some xfig eps files cannot be converted anymore
# => screenshot e.g., for F_0408.HCT_screen.png
# F... bug 2: black instead of white background 
# eliminate transparency: option -flatten
# (also solves black-bug in inserted xfig figure files)

#grep "\.eps" folien*.tex > epsfiles
#grep "\.eps" folien1_eng_old.tex > epsfiles
grep "\.eps" folien2_eng_old.tex > epsfiles
perl -i -p -e 's/^.+textwidth\}\{//g' epsfiles
perl -i -p -e 's/^.+textwidth\}\{//g' epsfiles
perl -i -p -e 's/^.+\{//g' epsfiles
perl -i -p -e 's/\}//g' epsfiles
# um Einzelzeilen mit 2em etc wegzubekommen
perl -i -p -e 's/^[^f].*\n//g' epsfiles 

myconvert='convert -flatten -units PixelsPerInch -density 150 -trim -matte -type PaletteMatte'

for f in `cat epsfiles`; do
    #echo "converting $f to png";
    path=`dirname "$f"`
    b=`basename "$f" .eps`
    fullb=$path/$b
    echo "converting to $fullb.png"
    $myconvert $fullb.eps $fullb.png
done
