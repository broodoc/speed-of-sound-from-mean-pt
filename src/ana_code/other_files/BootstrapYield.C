#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TMath.h>
#include <iostream>
#include <vector>

void BootstrapYield(TH1D* hpt_ref, int numBootstrap = 1000) {
    const int nBins = hpt_ref->GetNbinsX();
    std::vector<std::vector<double>> bootstrapYields(numBootstrap, std::vector<double>(nBins, 0.0));

    // Random generator
    TRandom3 random(0);

    // Perform bootstrapping
    for (int i = 0; i < numBootstrap; ++i) {
        for (int j = 1; j <= nBins; ++j) {
            // Resample bin contents with replacement
            int randomBin = random.Integer(nBins) + 1;
            bootstrapYields[i][j - 1] = hpt_ref->GetBinContent(randomBin);
        }
    }

    // Calculate the mean and standard deviation for each bin
    std::vector<double> meanYields(nBins, 0.0);
    std::vector<double> stdDevYields(nBins, 0.0);

    for (int j = 0; j < nBins; ++j) {
        // Collect yields for the current bin across all bootstrap samples
        std::vector<double> resampledYields(numBootstrap, 0.0);
        for (int i = 0; i < numBootstrap; ++i) {
            resampledYields[i] = bootstrapYields[i][j];
        }

        // Calculate mean and standard deviation for the current bin
        meanYields[j] = TMath::Mean(resampledYields.begin(), resampledYields.end());
        stdDevYields[j] = TMath::RMS(resampledYields.begin(), resampledYields.end());
    }

    // Print or plot results
    for (int j = 0; j < nBins; ++j) {
        std::cout << "Bin " << j + 1 << ": Mean Yield = " << meanYields[j]
                  << ", Std Dev = " << stdDevYields[j] << std::endl;
    }

    // Optionally create a new histogram with bootstrapped errors
    TH1D* hBootstrap = (TH1D*)hpt_ref->Clone("hBootstrap");
    for (int j = 1; j <= nBins; ++j) {
        hBootstrap->SetBinContent(j, meanYields[j - 1]);
        hBootstrap->SetBinError(j, stdDevYields[j - 1]);
    }
    
    TList* lOut = new TList();
    lOut->SetOwner(true);
    lOut->Add(hBootstrap);
    TFile* fOut = new TFile("hBootstrap.root","recreate");
    
    fOut->Write();
    fOut->cd();
    lOut->Write();
    delete fOut;
    delete lOut;
}

int main() {
    // Open the input ROOT file
    TFile* inputFile = TFile::Open("graphs_dEta1p0__pTlow0p20_hPtDec_Run21_from_raw_spectrum.root", "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open file!" << std::endl;
        return 1;
    }

    // Retrieve the hpt_ref histogram
    TH1D* hpt_ref = (TH1D*)inputFile->Get("(ref) 0-5%");
    if (!hpt_ref) {
        std::cerr << "Error: hpt_ref histogram not found in file!" << std::endl;
        inputFile->Close();
        return 1;
    }

    // Enable Sumw2 for proper error handling
    hpt_ref->Sumw2();

    // Perform bootstrapping
    BootstrapYield(hpt_ref);

    // Close the input file
    inputFile->Close();

    return 0;
}

