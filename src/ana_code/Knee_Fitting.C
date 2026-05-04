#include <TF1.h>
#include <TH1.h>
#include <TH1D.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TCanvas.h>
#include <iostream>

double P_n(const double* x, const double* par) {
    double n = x[0];            // n (x-axis value)
    double n_knee = par[0];     // Fit parameter: n_knee
    double a1 = par[1];         // Fit parameter: a1
    double a2 = par[2];         // Fit parameter: a2
    double a3 = par[3];         // Fit parameter: a3
    double sigma = par[4];      // Fit parameter: sigma

    // Integrating from 0 to 1 over c_b
    double integral = 0.0;
    int n_steps = 100;   //n number of integration steps
    double step = 1.0 / n_steps;

    for (int i = 0; i <= n_steps; ++i) {
        double c_b = i * step; // integrating variable
        double exponent = n_knee * TMath::Exp(-a1 * c_b - a2 * c_b * c_b - a3 * c_b * c_b * c_b);
        double arg = (n - exponent) / sigma;
        integral += TMath::Exp(-0.5 * arg * arg) / (sigma * TMath::Sqrt(2 * TMath::Pi()));
    }

    return integral * step; // multiplying step size to approximate the integral
}

void Knee_Fitting() {
    // https://arxiv.org/pdf/1708.00081 using scenario (B) sigma(c_b) = sigma(0) = sigma
    string dirfilename = "/Users/calebbroodo/dst_files/200GeV/combined/ptdep_CE/225MeV/11_v7_systematics/DCA1p8/temp";
    TFile *fIn = TFile::Open(dirfilename.c_str(), "READ");
    if (fIn->IsZombie()) {std::cout << "Error opening file" << std::endl; std::exit(EXIT_FAILURE);}
//    TH1D* h = (TH1D*)fIn->Get("dEta1p0__pTlow0p20_hRefMultCorr");
    TH1D* h = (TH1D*)fIn->Get("dEta1p0__pTlow0p20_hRefMult");
    h->Scale(1./h->Integral());
//    h->Sumw2();   enables error propagation
//    h->SetBinErrorOption(TH1::kPoisson);
    
    h->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT");
    h->GetYaxis()->SetTitle("Event fraction (a.u.)");
    
    //Modulate
//    double n_min = 500; // minimum Ncharged cutoff ( as recommended by https://arxiv.org/pdf/1708.00081 )
    double n_min = 1; // minimum Ncharged cutoff ( as recommended by https://arxiv.org/pdf/1708.00081 )
    double n_max = 800;
    
    cout << "n_min: " << n_min << " | n_max: " << n_max << endl;
    
    gErrorIgnoreLevel = kFatal;
    TList* lOut = new TList();
    lOut->SetOwner(true);
    
    int nBins = h->GetNbinsX();
    TH1D* h_uncorr = (TH1D*)h->Clone("h_uncorr");
    TRandom3 rand;
    //Unccorelated errors using Guassian input:
        for (int i = 1; i <= nBins; ++i) {
            double content = h->GetBinContent(i);
            double error = h->GetBinError(i);
            if (error < 0 || std::isnan(error)) {
                    std::cerr << "Invalid error in bin " << i << std::endl;
                }
            // Generate new bin content with Gaussian fluctuations
            double new_content = rand.Gaus(content, error);
            h_uncorr->SetBinContent(i, new_content);
            h_uncorr->SetBinError(i, error); // Preserve errors
        }
    
    lOut->Add(h);
    lOut->Add(h_uncorr);

    // Define the TF1 with the custom function
    TF1* fitFunc = new TF1("fitFunc", P_n, n_min, n_max, 5); // 5 parameters

    // Set parameter names and initial values
    fitFunc->SetParName(0, "N_{ch,knee}");
    fitFunc->SetParName(1, "a_{1}");
    fitFunc->SetParName(2, "a_{2}");
    fitFunc->SetParName(3, "a_{3}");
    fitFunc->SetParName(4, "#sigma_{0}");

    // Initial guesses: n_knee, a1, a2, a3, sigma_0
    fitFunc->SetParameters(1.1, 4, 0.95, -0.3, 35);
//    fitFunc->SetParameters(560, 4, 0.95, -0.3, 35);
    fitFunc->SetParLimits(0, 0.5, 1.5);
//    fitFunc->SetParLimits(1, 4, 5.0);
//    fitFunc->SetParLimits(2, 0.9, 1.1);
//    fitFunc->SetParLimits(3, -0.25, -0.35);
//    fitFunc->SetParLimits(4, 30, 40);
    

    h_uncorr->Fit(fitFunc, "R");
    
    double chi2 = fitFunc->GetChisquare();
    double ndf = fitFunc->GetNDF();
    std::cout << "Chi^2 / NDF: " << chi2 << " / " << ndf << " = "<<  chi2 / ndf << std::endl;

    TFile* fOut = new TFile("knee_fitting.root","recreate");
    
    fOut->Write();
    fOut->cd();
    lOut->Write();
    delete fOut;
    delete lOut;
    
    cout << "Fitting complete" << endl;
}

