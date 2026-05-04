#!/bin/bash

dir=$1
prepdir=$2
refMult=$3
pTMult=$4
InEx=$5
DebugVerbose=$6

echo "Checking environment 🖥️ "

ARCH=$(uname -m)

if [ "$ARCH" == "arm64" ]; then
    echo "This script cannot run on an ARM64 environment. Consult README.txt to reconfigure environment in working tab. Exiting..."
    exit 1
fi
echo "System architecture is $ARCH ✅. Proceeding..."

echo "${dir} (DST directory) passing ✅"
echo "${prepdir} (Prepped files directory) passing ✅"
echo "${refMult} (Mult name) passing ✅"
echo "${pTMult} (pTMult name) passing ✅"
echo "${InEx} (Inclusive/Exclusive setting) passing ✅"
echo "${DebugVerbose} (Debug setting) passing ✅"

echo Generating bootstrap statistics 📊...
#root -l -x <<EOF
#.L Bootstrap_Yield.Cl
#Bootstrap_Yield("${dir}","${prepdir}","${refMult}","${pTMult}","${InEx}","${DebugVerbose}")
#.q
#EOF

root -l -x <<EOF
.L Bootstrap_Yield.C+
Bootstrap_Yield("${dir}","${prepdir}","${refMult}","${pTMult}","${InEx}","${DebugVerbose}")
.q
EOF

echo "Statistics successfully produced ✅. See yield_errors/ in ${prepdir} for details"

cp -r "${prepdir}/yield_errors" /Users/calebbroodo/Desktop/sos_analysis/src/ana_code/

for i in {2..1}; do
    echo -ne "🚦 Moving on to \033[1;31mprobability density shading\033[0m in $i   \r"
    sleep 1
done

echo -e "🚦 Starting \033[1;31mprobability density shading\033[0m analysis...          "
sleep 1
