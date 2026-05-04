/*
  root -l -b
  gSystem->Load("/Users/nustream/Desktop/Class/Parameterization.so")
  .L FindPercentiles.cpp
  main()
*/

#include <cstdio>
#include <iostream>
#include <string>
#include <iostream>
#include <fstream>

#include "TSystem.h"
#include "Parameterization.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2D.h"

using std::string;

void FindAmplitudePercBinning(string system);
int FindPercentiles(string system);

int FindPercentiles(string system = "AuAu_3") {
    gSystem->Load("/Users/calebbroodo/Desktop/200GeV_Run16/Analysis/Analysis_Eta_pT_Comp/Parameterization.so");
    FindAmplitudePercBinning(system);
    return 0;
}

void FindAmplitudePercBinning(string system = "AuAu_3") {
    Parameterization per{};
    per.SetSystem(system);
    //  string nameInput = per.GetNameInputFile();
    
    cout << '\n' << "\t ---------------" << '\n';
    cout << "\t Sytem: " << per.GetSystem() << '\n';
    cout << "\t ---------------" << '\n';
    std::vector<std::string> fileNames = {
//        "/Users/calebbroodo/Desktop/3GeV/Analysis/Glauber_Histograms/bestChi2RootFile_dEta2p0__pTlow0p15"
//      "/Users/calebbroodo/Desktop/3GeV_Run21/Analysis/DST_Summary_Files/3GeVFXT_Run21_Analysis_SummarydEta1p0__pTlow0p15__Sept13.root"
//        "/Users/calebbroodo/Desktop/200GeV_Run16/Analysis/DST_Summary_Files/200_GeV_Run11_Analysis_SummaryEta0p5__pTlow0p15__June20.root"
        "/Users/calebbroodo/Desktop/200GeV_Run16/Analysis/DST_Summary_Files/Coarse_Study/Dec_8_2024/combined_data_200GeV_Run16_Eta0p5_pTlow0p20_SelfCorr_Dec_8_2024.root"
    };
    
//    string RefMult_label{"hRefMultSim"};
    string RefMult_label{"dEta1p0__pTlow0p20_hNch"};
    
//    "Eta0p5__pTlow0p15_hNch"
    std::vector<std::string> Nch_labels = {
        "dEta1p0__pTlow0p20_hNch"
        };
    int i = 0;
    for (const auto& fileName : fileNames) {
    TFile fIn =
    TFile(fileName.c_str(),"read");
    if (fIn.IsZombie()) {
        std::cout << "Error opening file" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    std::cout << '\n' << "Success opening file: " << fileName << endl;
    
    TH1D *hTPC = static_cast<TH1D *>(fIn.Get(RefMult_label.c_str()));
        hTPC->SetName(Nch_labels[i].c_str());
        cout << '\n' << "Nch label: " << Nch_labels[i].c_str();
        if(!hTPC) {cout << "hTPC does not exist!" << '\n'; std::exit(EXIT_FAILURE);}
    std::cout << '\n' << "Success opening histogram: " << Nch_labels[i].c_str() << endl;
    i++;
    
    bool isCoarseBinning{true};
    per.SetMultEstimator("TPC");
    //  cout << "\t --> COARSE PERCENTILE BINS: " << per.GetMultEstimator() << '\n';
    //  per.FindPercentiles("", hTPC, Parameterization::PercFirstBin::kLast, true,
    //                      isCoarseBinning);
    cout << '\n';
    cout << "\t --> FINE PERCENTILE BINS: " << per.GetMultEstimator() << '\n';
     if (i == 1){ per.FindPercentiles(hTPC, Parameterization::PercFirstBin::kLast, true, !isCoarseBinning); }
    else{ per.FindPercentiles(hTPC, Parameterization::PercFirstBin::kLast, true, !isCoarseBinning); }
        cout << '\n' << "----====----====----====----====----====----====----====----====----====----\n";
}
//    cout << per.GetCoarseAmplitudeBinCutLow(0) << " - " << per.GetCoarseAmplitudeBinCutHigh(0) << '\n';
//    cout << per.GetCoarseAmplitudeBinCutLow(1) << " - " << per.GetCoarseAmplitudeBinCutHigh(1) << '\n';

/*
  //---------------------------

  per.SetMultEstimator("TPC");
  cout << "\t --> COARSE PERCENTILE BINS: " << per.GetMultEstimator() << '\n';
  per.FindPercentiles("", hNchEtaPos, Parameterization::PercFirstBin::kLast,
                      true, isCoarseBinning);
  cout << "\t --> FINE PERCENTILE BINS: " << per.GetMultEstimator() << '\n';
  per.FindPercentiles("", hNchEtaPos, Parameterization::PercFirstBin::kLast,
                      true, !isCoarseBinning);

  //---------------------------

  per.SetMultEstimator("TracksEtaGapTPC");
  cout << "\t --> COARSE PERCENTILE BINS: " << per.GetMultEstimator() << '\n';
  per.FindPercentiles("", hNchEtaGapTPC, Parameterization::PercFirstBin::kLast,
                      true, isCoarseBinning);
  cout << "\t --> FINE PERCENTILE BINS: " << per.GetMultEstimator() << '\n';
  per.FindPercentiles("", hNchEtaGapTPC, Parameterization::PercFirstBin::kLast,
                      true, !isCoarseBinning);

  //---------------------------

  per.SetMultEstimator("TrackletsEtaGap");
  cout << "\t --> COARSE PERCENTILE BINS: " << per.GetMultEstimator() << '\n';
  per.FindPercentiles("", hTrackletsEtaGap,
                      Parameterization::PercFirstBin::kLast, true,
                      isCoarseBinning);
  cout << "\t --> FINE PERCENTILE BINS: " << per.GetMultEstimator() << '\n';
  per.FindPercentiles("", hTrackletsEtaGap,
                      Parameterization::PercFirstBin::kLast, true,
                      !isCoarseBinning);

  //---------------------------

  per.SetMultEstimator("Tracklets14");
  cout << "\t --> COARSE PERCENTILE BINS: " << per.GetMultEstimator() << '\n';
  per.FindPercentiles("", hTracklets14, Parameterization::PercFirstBin::kLast,
                      true, isCoarseBinning);
  cout << "\t --> FINE PERCENTILE BINS: " << per.GetMultEstimator() << '\n';
  per.FindPercentiles("", hTracklets14, Parameterization::PercFirstBin::kLast,
                      true, !isCoarseBinning);
  //---------------------------
*/
  // per.SetMultEstimator("Tracklets10");
  // cout << "\t ----> FINDING COARSE PERCENTILE BINNING Tracklets |Eta| < 1
  // <----"
  //      << '\n';
  // per.FindPercentiles("", hTracklets10,
  // Parameterization::PercFirstBin::kLast,
  //                     true, isCoarseBinning);
  //
  // cout << '\n';
  // cout << "\t ----> FINDING FINE PERCENTILE BINNING Tracklets |Eta| < 1
  // <----"
  //      << '\n';
  //
  // per.FindPercentiles("", hTracklets10,
  // Parameterization::PercFirstBin::kLast,
  //                     true, !isCoarseBinning);
  //
  // //---------------------------
  //
}
