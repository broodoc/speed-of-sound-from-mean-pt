#!/bin/bash

clear

echo -e "📊 Performing \033[0;33mprobability density shading\033[0m..."

dir=$1
prepdir=$2
refMult=$3
#pTMult=$4

cd /Users/calebbroodo/Desktop/sos_analysis/src/ana_code

./run_Multiplicity_Shading.sh ${dir} ${prepdir} ${refMult}
