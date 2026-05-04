
#include <TROOT.h>
#include "TFile.h"
#include "TStyle.h"
#include "TString.h"
#include <TH1.h>
#include <TH2.h>
#include "TLatex.h"
#include <TCanvas.h>
#include <TStyle.h>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <random>
#include "TRandom.h"
#include <TError.h>
#include <TGraph.h>
#include <TMath.h>
#include <TTree.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include "Parameterization.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

void Run_Analysis(string system = "AuAu_200");

void Run_Analysis(string system = "AuAu_200"){
    Parameterization per{};
    per.SetSystem(system);
    per.SetMultEstimator("TPC");
    std::vector<std::string> fileNames = {
//        "Eta0p5__pTlow0p15_bin_cuts.txt"
        "dEta1p0__pTlow0p20_bin_cuts.txt"
//        "dEta1p0__pTlow0p20_bin_cuts.txt",
//        "dEta1p0__pTlow0p40_bin_cuts.txt",
//        "dEta1p0__pTlow0p80_bin_cuts.txt",
//        "dEta1p3__pTlow0p15_bin_cuts.txt",
//        "dEta1p3__pTlow0p20_bin_cuts.txt",
//        "dEta1p3__pTlow0p40_bin_cuts.txt",
//        "dEta1p3__pTlow0p80_bin_cuts.txt",
//        "dEta1p6__pTlow0p15_bin_cuts.txt",
//        "dEta1p6__pTlow0p20_bin_cuts.txt",
//        "dEta1p6__pTlow0p40_bin_cuts.txt",
//        "dEta1p6__pTlow0p80_bin_cuts.txt",
//        "dEta2p0__pTlow0p15_bin_cuts.txt",
//        "dEta2p0__pTlow0p20_bin_cuts.txt",
//        "dEta2p0__pTlow0p40_bin_cuts.txt",
//        "dEta2p0__pTlow0p80_bin_cuts.txt"
    };
    for (const auto& filename : fileNames) {
        std::ifstream infile(filename);
        
        if (!infile.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            std::exit(EXIT_FAILURE);
        }
        
        std::string line;
        std::vector<int> BinCuts;
        
        if (std::getline(infile, line)) {
            std::stringstream ss(line);
            std::string item;
            
            while (std::getline(ss, item, ',')) {
                // Remove any leading or trailing whitespace from the item
                item.erase(0, item.find_first_not_of(" \t\n\r"));
                item.erase(item.find_last_not_of(" \t\n\r") + 1);
                
                if (!item.empty()) {
                    // Convert the item to an integer and add to the vector
                    BinCuts.push_back(std::stoi(item));
                }
            }
        }
        infile.close();
        
        // Print the numbers to verify
        int i = 0;
        for (const auto& num : BinCuts) {
            std::cout << num << ", ";
            i++;
        }
        std::cout << std::endl;
        per.SetFineBinCutsTPC(BinCuts.data());
        per.PrintFineBinCuts();
//        std::exit(EXIT_FAILURE);
        bool isCoarseBinning = true;
        std::string modifiedFilename = filename;
        if (modifiedFilename.length() >= 12) {
            // Replace the last four characters with "bin_cuts.txt"
            modifiedFilename.replace(modifiedFilename.length() - 12, 12, "_Oct07.root");
        }
//        /Users/calebbroodo/Desktop/200GeV_Run16/Analysis/DST_Summary_Files/combined_data_200GeV_Run16_Eta_and_pT_Sep_27_2024.root
//        string dirfilename = "/Users/calebbroodo/Desktop/200GeV_Run16/Analysis/DST_Summary_Files" + modifiedFilename;
        string dirfilename = "/Users/calebbroodo/Desktop/200GeV_Run16/Analysis/DST_Summary_Files/Coarse_Study/Dec_8_2024/combined_data_200GeV_Run16_Eta0p5_pTlow0p20_SelfCorr_Dec_8_2024.root";
        TFile *fIn = TFile::Open(dirfilename.c_str(), "READ");
        if (fIn->IsZombie()) {std::cout << "Error opening file" << std::endl; std::exit(EXIT_FAILURE);}
//        if (modifiedFilename.length() >= 10) { modifiedFilename.replace(modifiedFilename.length() - 12, 12, "_hNch");}
        if (modifiedFilename.length() >= 10) { modifiedFilename.replace(modifiedFilename.length() - 12, 12, "_hRefMultCorr");}
        TH1D* Nch_dist = (TH1D*)fIn->Get(modifiedFilename.c_str());
        if (!Nch_dist){cout << "\n\tError pulling Nch histogram from " << fIn->GetName() << '\n';
            cout << "\nModified Filename: " << modifiedFilename << '\n';
            std::exit(EXIT_FAILURE);}
//        if (modifiedFilename.length() >= 4) {modifiedFilename.replace(modifiedFilename.length() - 3, 3, "Pt_Nch");}
        if (modifiedFilename.length() >= 8) {modifiedFilename.replace(modifiedFilename.length() - 11, 11, "Pt_Nch");}
        TH2D* hPt_Nch =  (TH2D*)fIn->Get(modifiedFilename.c_str());
        if (!hPt_Nch){cout << "\nError pulling hPt_Nch histogram from " << fIn->GetName() << '\n';
            cout << "Modified Filename: " << modifiedFilename << '\n';
            std::exit(EXIT_FAILURE);}
        cout << '\n' << "Running analysis on " << dirfilename << "... " << '\n';
        per.Ratio_compute_pt_Nch(hPt_Nch, Nch_dist, !isCoarseBinning, per); ///Fine
//        std::exit(EXIT_FAILURE);
        std::cout << '\n' << "-=-==-=-==-=-==-=-==-=-==-=-==-=-==-=-==-=-==-=-==-=-==-=-==-"<< std::endl;
        gErrorIgnoreLevel = kFatal; //Ignoring useless and nonFatal writing-based errors
        std::cout << '\n' << "END OF Run_Analysis()"<< std::endl;
    }
    gErrorIgnoreLevel = kFatal; //Ignoring useless and nonFatal writing-based errors
}
