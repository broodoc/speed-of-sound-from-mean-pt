#!/bin/bash

echo "Checking environment 🖥️ "

ARCH=$(uname -m)

if [ "$ARCH" == "arm64" ]; then
    echo "This script cannot run on an ARM64 environment. Consult README.txt to reconfigure environment in working tab. Exiting..."
    exit 1
fi

echo "System architecture is $ARCH ✅. Proceeding..."

# Clean up previous builds
echo ----====----====----====----====----====----====----====----====----====----
echo LINKING DEPENDENCIES 🔧...
make clean

# Compile the code
make
echo DEPENDENCY LINKAGE DONE ✅...
echo ----====----====----====----====----====----====----====----====----====----

root -l -x <<EOF
.L FindPercentiles.cpp
FindPercentiles()
.q
EOF
echo ----====----====----====----====----====----====----====----====----====----
echo UPDATING DEPENDENCY LINKAGE 🔄...
make clean
make
echo DEPENDENCY LINKAGE UPDATED ✅...
echo ----====----====----====----====----====----====----====----====----====----
echo Running analysis 📊...
# Run ROOT script Run_Analysis.cpp
root -l -x <<EOF
.L Run_Analysis.cpp
Run_Analysis()
.q
EOF

echo Analysis completed ✅
echo Note: Rename the output file to the dataset configuration - eta[]_pTlow[]_[Date]_[Run no.]

cd ../Comparison_Plotting/

root -l -n <<EOF
.L MergeGraphsIntoGraphsOutput.C
MergeGraphsIntoGraphsOutput()
.q
EOF

cd ../Analysis_Eta_pT_Comp/

root -l graphs_output.root -e 'new TBrowser'
#root -l graphs_dEta1p0__pTlow0p20_hPtDec_Run21_from_raw_spectrum.root -e 'new TBrowser'
#root -l graphs_dEta1p0__pTlow0p20_hPtDec_Run21_from_raw_spectrum.root -e 'new TBrowser'
#dEta1p0__pTlow0p15_bin_cuts.txt
# Open the specified ROOT file in a new ROOT session
#root -l graphs_Eta0p5__pTlow0p15_hPtOct_08_Run11_from_raw_spectrum.root -e 'new TBrowser'
# root -l graphs_Eta0p5__pTlow0p20_hPtOct_07_from_raw_spectrum.root -e 'new TBrowser'
# root -l graphs_Eta0p5__pTlow0p15_hPtSep_27_from_raw_spectrum.root -e 'new TBrowser'
