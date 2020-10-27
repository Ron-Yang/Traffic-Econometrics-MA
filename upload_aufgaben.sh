
#!/bin/bash


# Syntax:
# upload_aufgaben.sh <latexname ohne ext> <optional -edit oder -e>
####################################################

if ( (($#<1)) || (($#>2)) ); then
  echo "Calling sequence:"
  echo " upload_aufgaben.sh <latexname ohne ext.> <optional -edit oder -e>"
  exit 0
fi

sourcedir=$PWD
targetdir="$HOME/public_html/Vkoek_Ma_Skript"
htmldir="$HOME/public_html/Vkoek_Ma_Skript"
texname="$1"

 psname=${texname}.ps
 pdfname=${texname}.pdf

 cd ${sourcedir}

 if test -r ${texname}.tex; then 

  echo "latexing ${texname} ..."
  latex ${texname} > logfile

  echo "export as pdf ..."
  dvips -Ppdf  -q -o $psname  ${texname}
  ps2pdf $psname;

 else 
  echo "error: File ${texname}.tex existiert nicht!!"
 fi; 

# upload to web page


echo "uploading to $targetdir/$pdfname and chmod o+r ..."
echo cp $pdfname $targetdir
cp $pdfname $targetdir
chmod o+r $targetdir/$pdfname

echo "uploaded skript part $1 to $targetdir/$pdfname:"
ls -l $targetdir/$pdfname
if (($#==2)); 
    then xemacs ${htmldir}/index.html;
   #(setze neue Links)
fi;

cd $targetdir

uploadweblocal.sh  $pdfname
if (($#==2)); 
  then cd $htmldir
  uploadweblocal.sh index.html
fi


