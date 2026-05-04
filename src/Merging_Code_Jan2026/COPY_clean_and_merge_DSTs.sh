#!/bin/bash

dir=$1
refMult=$2
pTMult=$3
glauber=$4
Multfromdata=$5

clear

echo "Inside clean_and_merge_DSTs"
    
# Find the first ROOT file in the directory
first_file=$(ls "${dir}"/output.*.root 2>/dev/null | head -n 1)
first_file_cleaned=$(ls "${dir}"/cleaned_output.*.root 2>/dev/null | head -n 1)

if [[ ! -f "$first_file" && ! -f "$first_file_cleaned" ]]; then
    echo "❌ No ROOT files found in ${dir}"
    exit 1
fi

# Check that both histograms exist in the ROOT file

root -l -q -b <<EOF | grep -q "Missing"
TFile *f = TFile::Open("${first_file}");
$( [[ "$glauber" == "N" ]] && echo 'if (!f->Get("'${refMult}'")) { std::cout << "Missing: '${refMult}'" << std::endl; }' )
#if [[ "$glauber" == "N" ]]; then
#if (!f->Get("${refMult}")) { std::cout << "Missing: ${refMult}" << std::endl; }
#fi
if (!f->Get("${pTMult}"))  { std::cout << "Missing: ${pTMult}" << std::endl; }
f->Close();
.q
EOF

if [[ $? -eq 0 ]]; then
    echo "❌ One or both histograms (${refMult}, ${pTMult}) not found in ${first_file}"
    exit 1
fi

echo "${dir} exists ✅"
echo "${refMult} (1D Mult histo) exists ✅"
echo "${pTMult} (2D pT-Mult histo) exists ✅"

echo "🧹 cleaning and merging DST files from ${dir}..."
sleep 2

cp -r delete_RunID.sh clean_runid_histos.C ${dir}
cd ${dir}

find . -type f -size 0 -exec rm -f {} \;

./delete_RunID.sh

if [[ -f temp ]]; then
    echo "⚠️  'temp' already exists. Skipping merge or consider removing it manually."
else
    echo "📦 Merging cleaned DSTs with hadd..."
    hadd temp cleaned_output* > /dev/null 2>&1
    echo "✅ Done. Merged output: temp"
fi

cd ${dir}

if [[ "$glauber" == "Y" ]]; then
    echo "📥 Adding centrality estimator to temp..."
#    root -l -x <<EOF > /dev/null 2>&1
#    root -l -x <<EOF >>
#    root -b -q <<EOF
    root -b -q <<'EOF'
#include <iostream>
#include <cstdio>
    std::cout << "test" << std::endl;
TFile *ftemp = TFile::Open("temp");

TH1D *hRefMult = (TH1D*)ftemp->Get("${Multfromdata}");
//TFile *fin = TFile::Open("bestChi2RootFile"); //Glauber model file
//TH1D *hR = (TH1D*)fin->Get("hRatio");         //Glauber model file
//TH1D *h = (TH1D*)fin->Get("hRefMultSim");     //Glauber model file
TFile *fin = TFile::Open("/Users/calebbroodo/Desktop/DST_summary_files/200GeV/combined/11_v6/temp"); //default 200 GeV Run11 refMultCorr
//TH1D *hR = (TH1D*)fin->Get("hRatio");
TH1D *h = (TH1D*)fin->Get("dEta1p0__pTlow0p20_hRefMult"); //default 200 GeV Run11
//if (h && hR && hRefMult) {
if (h && hRefMult) {
    std::cout << "Updating temp" << std::endl;
    TFile *fout = new TFile("temp", "UPDATE");
    TH1D *h_clone = (TH1D*)h->Clone("hRefMultCorrGlauber2"); //default 200 GeV Run11
//    TH1D *h_clone = (TH1D*)h->Clone("hRefMultCorrGlauber"); //Glauber model file
    TH1D *hRefMult_clone = (TH1D*)hRefMult->Clone("hRefMult_clone");

//    int nbins = h_clone->GetNbinsX();
//    for (int i = 1; i <= nbins; ++i) {
//        double x = h_clone->GetBinCenter(i);
//        if (x > 300) {
//            h_clone->SetBinContent(i, hRefMult_clone->GetBinContent(i));
//            h_clone->SetBinError(i, hRefMult_clone->GetBinError(i));
//        }
//    }
    h_clone->Write("hRefMultCorrGlauber2", TObject::kOverwrite); //default 200 GeV Run11
//    hR->Clone("hRatio")->Write("hRatio", TObject::kOverwrite);
    fout->Close();
}  else {
    std::cout << "❌ Histogram hRefMultSim not found in bestChi2RootFile." << std::endl;
}
fin->Close();
.q
EOF
fi

cd ${dir}

if [[ "$glauber" == "Y" ]]; then
    echo "📥 Adding centrality estimator to merged_*.root files..."
    for i in {1..10}; do
        file="merged_${i}.root"
        if [[ -f "$file" ]]; then
            root -l -x <<EOF > /dev/null 2>&1
TFile *fin = TFile::Open("temp");
TH1D *h = (TH1D*)fin->Get("hRefMultCorrGlauber2");
#TH1D *h = (TH1D*)fin->Get("hRefMultCorrGlauber");
if (h) {
    TFile *fout = new TFile("${file}", "UPDATE");
    TH1D *h_clone = (TH1D*)h->Clone("hRefMultCorrGlauber2");
#    TH1D *h_clone = (TH1D*)h->Clone("hRefMultCorrGlauber");
    h_clone->Write("hRefMultCorrGlauber2",TObject::kOverwrite);
#    h_clone->Write("hRefMultCorrGlauber",TObject::kOverwrite);
    fout->Close();
} else {
    std::cout << "❌ Missing hRefMultSim in bestChi2RootFile." << std::endl;
}
fin->Close();
.q
EOF
        fi
    done
fi

for i in {3..1}; do
    echo -ne "🚦 Moving on to \033[0;34mcentrality determination\033[0m  in $i   \r"
    sleep 1
done

cd /Users/calebbroodo/Desktop/sos_analysis/src
