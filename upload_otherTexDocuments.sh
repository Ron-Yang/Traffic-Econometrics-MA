
#!/bin/bash


# Syntax:
# upload_aufgaben.sh <latexname ohne ext> <optional -edit oder -e>
####################################################

if ( (($#<1)) || (($#>2)) ); then
  echo "Calling sequence:"
  echo " upload_otherTexDocuments.sh <latexname ohne ext.> <optional -edit oder -e>"
  exit 0
fi

sourcedir=$PWD
targetdir="$HOME/public_html/professional/Vkoek_Ma_Skript"
htmldir="$HOME/public_html/professional/Vkoek_Ma"
texname="$1"

 pdfname=${texname}.pdf 

 cd ${sourcedir}

 if test -r ${texname}.tex; then 

  echo "pdflatexing ${texname} ..."
  pdflatex ${texname} > logfile

 else 
  echo "error: File ${texname}.tex existiert nicht!!"
 fi; 

# upload to web page


echo "uploading to $targetdir/$pdfname and chmod o+r ..."
echo cp $pdfname $targetdir
cp $pdfname $targetdir
chmod o+r $targetdir/$pdfname

if (($#==2)); 
    then emacs ${htmldir}/index.html;
     #(setze neue Links)
    else echo "potentially edit html file ${htmldir}/index.html"
fi


