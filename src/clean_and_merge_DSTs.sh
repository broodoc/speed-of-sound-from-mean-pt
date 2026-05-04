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
first_file_cleaned=$(ls "${dir}"/cleaned_output*.root 2>/dev/null | head -n 1)
#first_file_cleaned=$(ls "${dir}"/cleaned_output*.root 2>/dev/null | head -n 1)

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
    #TFile *ftemp = TFile::Open("temp");
#    TH1D *hRefMult = (TH1D*)ftemp->Get("${Multfromdata}"); using DCA1p5 temp instead for stable centrality selection across all variations
#TFile *ftemp = TFile::Open("/Users/calebbroodo/Desktop/DST_summary_files/200GeV/combined/11_v6_systematics/DCA1p5_default/temp");
#TFile::Open("/Users/calebbroodo/Desktop/DST_summary_files/200GeV/combined/11_v6_systematics/Glauber_fit/bestChi2RootFile");
#54GeV Analysis:
#TFile *ftemp = TFile::Open("/Users/calebbroodo/Desktop/DST_summary_files/54p4GeV/Run17/DCA2p0/17_v2/temp");
#TFile *fin = TFile::Open("/Users/calebbroodo/Desktop/DST_summary_files/54p4GeV/Run17/DCA2p0/Glauber_fit/bestChi2RootFile");
#27GeV Analysis:
#TFile *ftemp = TFile::Open("/Users/calebbroodo/Desktop/DST_summary_files/27GeV/Run18/DCA2p0/18_v2/temp");
#TFile *fin = TFile::Open("/Users/calebbroodo/Desktop/DST_summary_files/27GeV/Run18/DCA2p0/Glauber_fit/bestChi2RootFile");
#TFile *ftemp = TFile::Open("/Users/calebbroodo/Desktop/DST_summary_files/54p4GeV/Run17/DCA2p0/17_v2/temp");
#TFile *fin = TFile::Open("/Users/calebbroodo/Desktop/DST_summary_files/54p4GeV/Run17/DCA2p0/Glauber_fit/bestChi2RootFile");

#        double switch = 250; # 200, 54 GeV
#        double switch = 30; # 27 GeV

ftemp="${dir}/temp"
#fin="${dir}/../../Glauber_fit/bestChi2RootFile"

if [[ -f "${dir}/../Glauber_fit/bestChi2RootFile" ]]; then
    echo "✅ bestChi2RootFile found in ../"
    fin="${dir}/../Glauber_fit/bestChi2RootFile"
    echo "$fin"
elif [[ -f "${dir}/../../Glauber_fit/bestChi2RootFile" ]]; then
    echo "✅ bestChi2RootFile found in ../../"
    fin="${dir}/../../Glauber_fit/bestChi2RootFile"
    echo "$fin"
else
    echo "ERROR: bestChi2RootFile not found in ../ or ../../"
    exit 1
fi

    root -l -x <<EOF > /dev/null 2>&1
TString ftemp_path = "${ftemp}";
TString fin_path   = "${fin}";
TFile *ftemp = TFile::Open(ftemp_path);
TFile *fin   = TFile::Open(fin_path);

TH1D *hRefMult = (TH1D*)ftemp->Get("dEta1p0__pTlow0p20_hRefMult");

TH1D *hR = (TH1D*)fin->Get("hRatio");
TH1D *h = (TH1D*)fin->Get("hRefMultSim");
if (h && hRefMult) {
    std::cout << "Updating temp" << std::endl;
    TFile *fout = new TFile("temp", "UPDATE");
    TH1D *h_clone = (TH1D*)h->Clone("hRefMultCorrGlauber");
    TH1D *hRefMult_clone = (TH1D*)hRefMult->Clone("hRefMult_clone");
    
    int nbins = h_clone->GetNbinsX();
    for (int i = 1; i <= nbins; ++i) {
        double x = h_clone->GetBinCenter(i);
        double sw = 250;
        if (x > sw) {
            h_clone->SetBinContent(i, hRefMult_clone->GetBinContent(i));
            h_clone->SetBinError(i, hRefMult_clone->GetBinError(i));
        }
    }

    h_clone->Write("hRefMultCorrGlauber", TObject::kOverwrite);
    fout->Close();
}  else {
    std::cout << "❌ Histogram hRefMultSim not found in bestChi2RootFile." << std::endl;
}
fin->Close();
finData->Close();
.q
EOF
fi

cd ${dir}

if [[ "$glauber" == "Y" ]]; then
  echo "📥 Adding centrality estimator to bootstrap/*.root files..."

  shopt -s nullglob
  for file in bootstrap/*.root; do
    [[ -f "$file" ]] || continue
    echo "  → updating $file"

    root -l -x <<EOF > /dev/null 2>&1
TFile *fin = TFile::Open("temp");
TH1D *h = (TH1D*)fin->Get("hRefMultCorrGlauber");
if (h) {
    TFile *fout = new TFile("${file}", "UPDATE");
    TH1D *h_clone = (TH1D*)h->Clone("hRefMultCorrGlauber");
    h_clone->SetDirectory(fout); // good hygiene
    h_clone->Write("hRefMultCorrGlauber", TObject::kOverwrite);
    fout->Close();
    delete fout;
} else {
    std::cout << "❌ Missing hRefMultCorrGlauber in temp." << std::endl;
}
fin->Close();
delete fin;
.q
EOF

  done
fi

#if [[ "$glauber" == "Y" ]]; then
#    echo "📥 Adding centrality estimator to merged_*.root files..."
#    for i in {1..10}; do
#        file="merged_${i}.root"
#        if [[ -f "$file" ]]; then
#            root -l -x <<EOF > /dev/null 2>&1
#TFile *fin = TFile::Open("temp");
#TH1D *h = (TH1D*)fin->Get("hRefMultCorrGlauber");
#if (h) {
#    TFile *fout = new TFile("${file}", "UPDATE");
#    TH1D *h_clone = (TH1D*)h->Clone("hRefMultCorrGlauber");
#    h_clone->Write("hRefMultCorrGlauber",TObject::kOverwrite);
#    fout->Close();
#} else {
#    std::cout << "❌ Missing hRefMultSim in bestChi2RootFile." << std::endl;
#}
#fin->Close();
#.q
#EOF
#        fi
#    done
#fi

for i in {3..1}; do
    echo -ne "🚦 Moving on to \033[0;34mcentrality determination\033[0m  in $i   \r"
    sleep 1
done

cd /Users/calebbroodo/Desktop/sos_analysis/src
