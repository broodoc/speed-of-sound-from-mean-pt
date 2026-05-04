#!/bin/bash

clear

echo -e "🧮 Performing \033[0;34mcentrality determination\033[0m..."

dir=$1
prepdir=$2
refMult=$3
overwrite=$4

cp -r cleanup.sh ${prepdir}

if [ "$overwrite" == "Y" ]; then
    echo "🧹 Overwriting previous files in $prepdir"
    (cd "$prepdir" && ./cleanup.sh)
fi

cd /Users/calebbroodo/Desktop/sos_analysis/src/ana_code

./percentilization_run.sh ${dir} ${prepdir} ${refMult}
