#!/bin/bash

dir=$1
prepdir=$2
refMult=$3

rm -f *bin_cuts.txt
rm -f *bin_cutsError.txt

echo "Checking environment 🖥️ "

ARCH=$(uname -m)

if [ "$ARCH" == "arm64" ]; then
    echo "This script cannot run on an ARM64 environment. Consult README.txt to reconfigure environment in working tab. Exiting..."
    exit 1
fi

echo "System architecture is $ARCH ✅. Proceeding..."

# Clean up previous builds
make clean

# Compile the code
make

root -l -x <<EOF
.L FindPercentiles.cpp+
FindPercentiles("${dir}", "${refMult}");
.q
EOF

cp -r *.txt "${prepdir}/"

echo "✅ centrality bin cuts + errors have been copied into ${prepdir}"

for i in {3..1}; do
    echo -ne "🚦 Moving on to \033[0;33mstatistical bootstrapping\033[0m in $i   \r"
    sleep 1
done
