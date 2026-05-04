#!/bin/bash

# Clean up previous builds
make clean

# Compile the code
make

root -l -x <<EOF
.L FindPercentiles.cpp
FindPercentiles()
.q
EOF
make clean

make

# Run ROOT script Run_Analysis.cpp
root -l -x <<EOF
.L Run_Analysis.cpp
Run_Analysis()
.q
EOF
#dEta1p0__pTlow0p15_bin_cuts.txt
# Open the specified ROOT file in a new ROOT session
#root -l graphs_Eta0p5__pTlow0p15_hPtOct_08_Run11_from_raw_spectrum.root -e 'new TBrowser'
# root -l graphs_Eta0p5__pTlow0p20_hPtOct_07_from_raw_spectrum.root -e 'new TBrowser'
# root -l graphs_Eta0p5__pTlow0p15_hPtSep_27_from_raw_spectrum.root -e 'new TBrowser'
