#!/bin/bash

clear

echo -e "📊 Performing \033[0;33mstatistical bootstrapping\033[0m..."
        
dir=$1
prepdir=$2
refMult=$3
pTMult=$4
InEx=$5
DebugVerbose=$6

cp -r /Users/calebbroodo/Desktop/sos_analysis/src/Merging_Code/* ${dir}

cd ${dir}

./run_merging.sh

cd /Users/calebbroodo/Desktop/sos_analysis/src/ana_code

./run_Bootstrapping.sh ${dir} ${prepdir} ${refMult} ${pTMult} ${InEx} ${DebugVerbose}
