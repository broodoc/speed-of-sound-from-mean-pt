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

void FindAmplitudePercBinning(string dir, string RefMult_label);
int FindPercentiles(string dir, string RefMult_label);

int FindPercentiles(string dir = "", string RefMult_label = "") {
    string system = "AuAu_200";
    gSystem->Load("/Users/calebbroodo/Desktop/sos_analysis/src/ana_code/Parameterization.so");
    FindAmplitudePercBinning(dir, RefMult_label);
    return 0;
}

void FindAmplitudePercBinning(string dir = "", string RefMult_label = "") {
    string system = "AuAu_200";
    Parameterization per{};
    per.SetSystem(system);
    
    cout << '\n' << "\t ---------------" << '\n';
    cout << "\t Sytem: " << per.GetSystem() << '\n';
    cout << "\t ---------------" << '\n';
    
    if (!dir.empty() && dir.back() != '/') dir += '/';
        
    std::string fileName = dir + "temp";
    
        TFile fIn = TFile(fileName.c_str(),"read");
        if (fIn.IsZombie()) {
            std::cout << "Error opening file" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        std::cout << '\n' << "Success opening file: " << fileName << endl;
        
        TH1D *hTPC = static_cast<TH1D *>(fIn.Get(RefMult_label.c_str()));
        hTPC->SetName(RefMult_label.c_str());
        cout << '\n' << "Nch label: " << RefMult_label.c_str();
        if(!hTPC) {cout << "hTPC does not exist!" << '\n'; std::exit(EXIT_FAILURE);}
        std::cout << '\n' << "Success opening histogram: " << RefMult_label.c_str() << endl;
        
        bool isCoarseBinning{true};
        per.SetMultEstimator("TPC");
        cout << '\n';
        cout << "\t --> FINE PERCENTILE BINS: " << per.GetMultEstimator() << '\n';
        per.FindPercentiles(hTPC, Parameterization::PercFirstBin::kLast, true, !isCoarseBinning);
        cout << '\n' << "----====----====----====----====----====----====----====----====----====----\n";
}
