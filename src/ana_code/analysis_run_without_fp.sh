#!/bin/bash

dir=$1
prepdir=$2
refMult=$3
pTMult=$4
InEx=$5
DebugVerbose=$6
DisplayFitRanges=$7
GenerateFigures=$8
MomentAna=$9

echo "Checking environment 🖥️ "

ARCH=$(uname -m)

if [ "$ARCH" == "arm64" ]; then
    echo "This script cannot run on an ARM64 environment. Consult README.txt to reconfigure environment in working tab. Exiting..."
    exit 1
fi

echo "System architecture is $ARCH ✅. Proceeding..."

clear

echo "${dir} (DST directory) passing ✅"
echo "${prepdir} (Prepped files directory) passing ✅"
echo "${refMult} (Mult name) passing ✅"
echo "${pTMult} (pTMult name) passing ✅"
echo "${InEx} (Inclusive/Exlsusive Setting) passing ✅"
echo "${DebugVerbose} (Debug verbose setting) passing ✅"
echo "${DisplayFitRanges} (pT fit ranges result display setting) passing ✅"
echo "${GenerateFigures} (Figure generation setting) passing ✅"
echo "${MomentAna} (Moment analysis setting) passing ✅"

echo Running analysis 📊...
root -l -x <<EOF
.L Run_Analysis.cpp
Run_Analysis("${dir}","${prepdir}","${refMult}","${pTMult}","${InEx}","${DebugVerbose}","${DisplayFitRanges}", "${MomentAna}")
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

cd ../ana_code/

prepdir="${prepdir%/}"

if [[ "$GenerateFigures" == "Y" ]]; then
    echo "Generating figures 📊"
    figdir="${prepdir}/Ana_Figures"
    mkdir -p "$figdir"
#    root -l -b -q "view_all_graphs_BestChi2NDF.C(\"graphs_output.root\", \"$figdir\")"
    root -l -b -q "view_all_graphs_pT_dependence.C(\"graphs_output.root\", \"$figdir\")"
fi

newname="${prepdir}/graphs_output_$(date +"%Y-%b-%d_%I-%M%p").root"

echo "Moving graphs_output.root to: $newname"

mv graphs_output.root "$newname"

if [[ "$MomentAna" == "Y" ]]; then

echo ".    .                      .       .              .                  "
echo "|\  /|                     _|_     / \             |           o      "
echo "| \/ | .-. .--.--. .-. .--. |     /___\  .--. .-.  | .  ..--.  .  .--."
echo "|    |(   )|  |  |(.-' |  | |    /     \ |  |(   ) | |  |\\--.  |  \\--."
echo "'    ' \`-' '  '  \`-\`--''  \`-\`-' '       \`'  \`-\`-' \`-\`--|\\--'-' \`-\`--'"
echo "                                                        ;             "
echo "                                                     \`-'              "
echo "  ................................................................................................  "
echo "  ..                                                                                                .."
echo "  ..                                                                                                .."
echo "  ..                                         -:..:.......                                           .."
echo "  ..                                         .............                                          .."
echo "  ..                              .          ....-::......                                          .."
echo "  ..                           .-..-.           ..=+-..           .---.                             .."
echo "  ..                          .-....::.       ..=.....+:         .-....:                            .."
echo "  ..                         .:.......:      .=:.......:#.      ::.....-.                           .."
echo "  ..                         :.........-    .-..........:-.    .:.......=.                          .."
echo "  ..                        .:..........-..:=.............:.  :..........-                          .."
echo "  ..                       .=............-.-..............:+.-...........::                         .."
echo "  ..     .........:........:..............*-...............:*.............:.  .    . .              .."
echo "  ..     . .. .............:.............-::=-...........:=:-:.............=.:::-::.:.::::::::.     .."
echo "  ..          ..........:.-..............=::::--:......:=::::-..............-............           .."
echo "  ..                    .-..............*::::::::=-.:-=:::::::=.............:.                      .."
echo "  ..                    -..............*::::::::::=+-::::::::::#.............=.                     .."
echo "  ..                   ...............:::::::::-=::::::-::::::::=.............-.                    .."
echo "  ..                 ..-.............:::::::-=::::::::::::-::::::*............:-                    .."
echo "  ..                 .-.............=::::=-:::::::::::::::::::-:::+............:.                   .."
echo "  ..                .-............-=:==-:::::::::::::::::::::::::---:...........-.                  .."
echo "  ..              ..-..........:=*=-:::::::::::::::::::::::::::::::::=+..........::.                .."
echo "  ..           ..--:------::..:=-----:-=---:::::::::::::::::::::::::::::+..........:--..            .."
echo "  ..  .........................................................::::::::...........................  .."
echo "  ..                                           ..    ..                                             .."
echo "  ..                                            :--....:..                                          .."
echo "  ..                                                                                                .."
echo "  ...::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::."

if [[ "$GenerateFigures" == "Y" ]]; then
    echo "Generating figures for moment analysis 📊"
    figdir="${prepdir}/Ana_Figures_Moment_Analysis"
    mkdir -p "$figdir"
    root -l -b -q "view_all_graphs_Moment.C(\"Momentgraphs_output.root\", \"$figdir\")"
fi

momentname="${prepdir}/Moment_graphs_output_$(date +"%Y-%b-%d_%I-%M%p").root"

echo "Moving Momentgraphs_output.root to: $momentname"

mv Momentgraphs_output.root "$momentname"
fi

fitname="${prepdir}/Fitlogs_$(date +"%Y-%b-%d_%I-%M%p")/"

mv Fitlogs "${fitname}"

sosfitname="${prepdir}/SOS_Fitlogs_$(date +"%Y-%b-%d_%I-%M%p")/"

mv SOS_Fitlogs "${sosfitname}"

echo -e "Navigate to \033[0;32m${prepdir}/\033[0m"
