
#! /bin/bash


# Syntax:
# latexChapter.sh <part_nr> 
####################################################


skriptname="Verkehrsoekonometrie_Ma"
skriptname_work="Verkehrsoekonometrie_work"
allParts="vkoek_Teil1,vkoek_Teil2,vkoek_Teil3,vkoek_Teil4,vkoek_Teil5,vkoek_Teil6,vkoek_Teil7"

logfile="logfile"

if (($#<1)) ; then
  echo "Calling sequence: latexChapter.sh <part number or all>"
  exit -1;
fi

chapterName=vkoek_Teil$1
if [ $1 == "all" ]; then chapterName=$skriptname; fi

if test -r ${skriptname}.tex; 
  then  cp -v ${skriptname}.tex ${skriptname_work}.tex;
  else  echo "error: File ${skriptname}.tex existiert nicht!"; exit -1;
fi; 

###################################

if [ $1 != "all" ]; then 
   if test -r $chapterName.tex; then
      echo "select  ${chapterName} to be included and latex  ${skriptname} ..."
      perl -i -p -e "s/\\includeonly.*/\\includeonly\{${chapterName}\}/g" ${skriptname_work}.tex
   else 
      echo "error: File $chapterName.tex existiert nicht!"; exit -1;
   fi;

else
   perl -i -p -e "s/\\includeonly.*/\\includeonly\{${allParts}\}/g" ${skriptname_work}.tex
   cp -v ${skriptname_work}.tex ${skriptname}_all.tex;
fi

echo "latexing (log output in $logfile) ..."
latex ${skriptname_work}
echo "produced ${skriptname_work}.dvi"
exit

echo "export as $chapterName.ps and show it ..."
dvips -Ppdf  -q -o $chapterName.ps  ${skriptname_work}

gv  -orientation=portrait $chapterName.ps 

exit

