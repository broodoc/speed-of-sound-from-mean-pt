#!/bin/bash

clear

echo -e "📈 Performing \033[0;31mspeed of sound\033[0m calculation (looks like this 📈)..."

dir=$1
prepdir=$2
refMult=$3
pTMult=$4
InEx=$5
DebugVerbose=$6
DisplayFitRanges=$7
GenerateFigures=$8
MomentAna=$9

cd /Users/calebbroodo/Desktop/sos_analysis/src/ana_code

./analysis_run_without_fp.sh ${dir} ${prepdir} ${refMult} ${pTMult} ${InEx} ${DebugVerbose} ${DisplayFitRanges} ${GenerateFigures} ${MomentAna}
