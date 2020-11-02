#!/bin/bash

# just simplifies the arg sequence of the mkvmerge command
# e.g. merge.sh tutorial01.mkv tutorial01_*.mkv
# will be expanded to
# mkvmerge -o tutorial01.mkv tutorial01_01.mkv +tutorial01_02.mkv +tutorial01_03.mkv +tutorial01_03.mkv +tutorial01_04.mkv

if [ $# -le 1 ]
then
    echo "Calling sequence: merge.sh <output video><input videos>"
    echo "example: merge.sh tutorial01.mkv tutorial01_*.mkv"
    exit 0;
fi
output=$1
shift
firstinput=$1
shift
remaininginput=$@
commandstring="mkvmerge -o $output $firstinput"
for f in $remaininginput; do
    commandstring="$commandstring +$f"
done
echo "commandstring=$commandstring"
$commandstring
echo "produced $output"

