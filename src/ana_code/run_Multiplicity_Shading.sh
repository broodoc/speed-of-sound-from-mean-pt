#!/bin/bash

dir=$1
prepdir=$2
refMult=$3
#pTMult=$4

cd ${prepdir}

if [[ -f pdf_shaded.pdf ]]; then
    echo "⚠️  'pdf_shaded.pdf' already exists. Skipping pdf shading or consider removing it manually."
fi
#else
cd /Users/calebbroodo/Desktop/sos_analysis/src/ana_code/
echo "Checking environment 🖥️ "

ARCH=$(uname -m)

    if [ "$ARCH" == "arm64" ]; then
        echo "This script cannot run on an ARM64 environment. Consult README.txt to reconfigure environment in working tab. Exiting..."
        exit 1
    fi

    echo "System architecture is $ARCH ✅. Proceeding..."

if [[ ! -d "$dir" ]]; then
    echo "❌ Error: DST directory '$dir' does not exist."
    exit 1
else
    echo "$dir (DST directory) passing ✅"
fi

if [[ ! -d "$prepdir" ]]; then
    echo "❌ Error: Prepped files directory '$prepdir' does not exist."
    exit 1
else
    echo "$prepdir (Prepped files directory) passing ✅"
fi

    echo "${refMult} (Mult name) passing ✅"
    #echo "${pTMult} (pTMult name) passing ✅"

echo "Running mult shading code"
root -l -x <<EOF
.L Shade_MultHisto.C
Shade_MultHisto("${dir}","${prepdir}","${refMult}")
.q
EOF

    echo "Shaded probability distribution successfully produced ✅. See Shaded_PDF.root in ${prepdir} for details"
#fi

for i in {2..1}; do
    echo -ne "🚦 Moving on to \033[1;31mspeed of sound\033[0m analysis in $i   \r"
    sleep 1
done

echo -e "🚦 Starting \033[1;31mspeed of sound\033[0m analysis...          "
sleep 2


cd /Users/calebbroodo/Desktop/sos_analysis/src
