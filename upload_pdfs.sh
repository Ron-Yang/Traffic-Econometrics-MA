
#!/bin/bash


# Syntax:
# upload_pdfs.sh <pdf files>
####################################################

sourcedir=$PWD
targetdir="$HOME/public_html/professional/Vkoek_Ma_Skript"
htmldir="$HOME/public_html/professional/Vkoek_Ma"

if ( (($#<1)) ); then
  echo "Calling sequence:"
  echo " upload_pdfs.sh <pdf files>"
  echo " will be uploaded in $targetdir"
  exit 0
fi

for f in "$@"; do

    if test -r $f;
    then 
	cp $f $targetdir
	chmod o+r $targetdir/$f
	echo "uploaded to $targetdir/$f"
    else
	echo "error: File ${texname}.tex existiert nicht!!"
    fi
done

echo "potentially edit $htmldir/index.html"


