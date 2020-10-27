
#!/bin/bash

 
# Syntax:
# upload_skript.sh <part_nr> <optional -edit oder -e>
####################################################

sourcedir=$PWD
targetdir="$HOME/public_html/professional/Vkoek_Ma_Skript"
htmldir="$HOME/public_html/professional/Vkoek_Ma"

if ( (($#<1)) || (($#>2)) ); then
  echo "uploads lecture skript or transparencies"
  echo "Calling sequence: upload_skriptVkoek.sh <latexname w/o ext><optional -edit oder -e>"
  echo  "ps file must be generated before calling this script!"
  exit -1
fi
 
latexName=$1
psname=$latexName.ps
pdfname=$latexName.pdf

cd ${sourcedir}
if test -r ${latexName}.ps;
  then  echo "generating $pdfname from $psname ...";  ps2pdf $psname;
  else  echo "error: File ${latexName}.ps nonexistent! run latexChapter.sh before this script";
  exit -1;
fi; 

# upload to web page

echo "uploading to $targetdir/$pdfname ..."
echo cp $pdfname $targetdir
cp $pdfname $targetdir
chmod o+r $targetdir/$pdfname

echo "uploaded skript part $1 to $targetdir/$pdfname:"
ls -l $targetdir/$pdfname
if (($#==2)); 
    then emacs ${htmldir}/index.html;
   #(setze neue Links)
fi;

cd $targetdir

# upload using filezilla

echo "enter in filezilla: sftp://mtreiber.de , p537815 , schwerster Onsight mit frz Grad"
echo "or directly filezilla sftp://mtreiber.de"
echo "if updating html as well: switch additionally to html file in $htmldir"
#filezilla --local=$targetdir


