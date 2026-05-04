#include <TFile.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>

double RightKneeFunction(double* x, double* par) {
    double n = x[0];
    double cs2   = par[0];  // speed of sound squared
    double knee  = par[1];  // location of transition
    double sigma = par[2];  // smearing width

    double norm     = sigma * TMath::Sqrt(2.0 / TMath::Pi());
    double exponent = TMath::Exp(-TMath::Power(n - knee, 2) / (2.0 * sigma * sigma));
    double erfc     = TMath::Erfc((n - knee) / (TMath::Sqrt(2.0) * sigma));

    
//    if (erfc < 1e-10) return 1e9;  // avoid division by near-zero
    double f = n - (norm * exponent / erfc);
//    if (f <= 1e-5) return 1e9;     // unphysical

//    if (erfc < 1e-10 || f <= 1e-5) {
//        TF1::RejectPoint();
//        return 0;
//    }
    
    return TMath::Exp(cs2 * TMath::Log(n / f));
}

void sample_fit() {
//    TFile* f = new TFile("CMS_5p02TeV_Measurement_arXiv2401.06896.root", "READ"); //CMS
    TFile* f = new TFile("/Users/calebbroodo/sos_analysis/prepped_files/200COL_Run10_Run11/200COL_Run11_Run10_refMult2/graphs_output_2025-Apr-28_12-41PM.root", "READ"); // STAR 200 GeV refMult2 Run10 + Run11
//    TFile* f = new TFile("Trajectum_output.root", "READ");  //OLD Trajectum
    if (!f || f->IsZombie()) {
        std::cerr << "❌ Could not open file!" << std::endl;
        return;
    }
//    TGraphErrors* g = (TGraphErrors*)f->Get("Mean_pT_vs_Centrality"); //Trajectum OLD 200 GeV
//    TGraphErrors* g = (TGraphErrors*)f->Get("g_dNdeta"); //CMS 5 TeV
    std::cout << "📈 Retrieving 200 GeV refMult2 for comparison" << std::endl;
    TDirectory *gdir1 = (TDirectory*)f->Get("pT study/Normalized (to ref class) quantities/pT-Integrated Nch/fit_range_0.00-0.40_GeV");
    TKey *key1 = gdir1->GetKey("0 #leq p_{T} #leq 10 (GeV/c)");
    TGraphErrors *g = (TGraphErrors*)key1->ReadObject<TGraphErrors>(); //STAR 200 GeV refMult2 Run10 + Run11
    if (!g) {
        std::cerr << "❌ Could not retrieve TGraphErrors!" << std::endl;
        return;
    }

    TF1* fitFunc = new TF1("RightKnee", RightKneeFunction, 1.05, 1.25, 3);
    fitFunc->SetParameter(0, 0.25);      // initial cs²
    fitFunc->SetParameter(1, 1.1);
    fitFunc->SetParameter(2, 0.015);
    fitFunc->SetParLimits(0, 0.2, 0.35);
    fitFunc->SetParLimits(1, 1.0, 1.2);
    fitFunc->SetParLimits(2, 0.005, 0.02);
//    fitFunc->FixParameter(1, 1.11);      // fixed knee
//    fitFunc->FixParameter(2, 0.0592);   // fixed sigma
//    fitFunc->FixParameter(2, 0.015);

    fitFunc->SetLineColor(kRed+1);
    fitFunc->SetLineWidth(3);
    fitFunc->SetNpx(10000);

    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);  // default is 5000
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-6);           // default is 1e-2
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");      // ensure Minuit2
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);               // 0 = fast, 2 = thorough
    
    g->SetMarkerStyle(25); // Open squares
    g->SetMarkerColor(kBlack);
    g->Draw("AP");
    g->Fit(fitFunc, "R");
    
    // Output fit results
    double cs2 = fitFunc->GetParameter(0);
    double cs2_err = fitFunc->GetParError(0);
    std::cout << "✅ Fit complete:" << std::endl;
    std::cout << "c_s^2 = " << cs2 << " ± " << cs2_err << std::endl;

    g->SetTitle("Sample fit;#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%};#LTp_{T}#GT/#LTp_{T}#GT_{0-5%}");
}

//These commands yield cs2 ~ 0.20 (but with poor Chi2):
//TF1* fitFunc = new TF1("RightKnee", RightKneeFunction, 1.05, 1.25, 3);
//fitFunc->SetParameter(0, 0.25);      // initial cs²
//fitFunc->SetParameter(1, 1.1);
//fitFunc->SetParameter(2, 0.015);
//fitFunc->SetParLimits(0, 0.2, 0.35);
//fitFunc->SetParLimits(1, 1.0, 1.2);
//fitFunc->SetParLimits(2, 0.005, 0.02);
//****************************************
//Minimizer is Minuit2 / Migrad
//Chi2                      =      3735.86
//NDf                       =            9
//Edm                       =    7.762e-05
//NCalls                    =          340
//p0                        =          0.2   +/-   6.66614e-06       (limited)
//p1                        =      1.04087   +/-   1.58425e-05       (limited)
//p2                        =         0.02   +/-   1.8717e-06        (limited)
//****************************************


//These commands yield the same cs2 ~ 0.24 as CMS:
//fitFunc->SetParameter(0, 0.25);      // initial cs²
//fitFunc->SetParameter(1, 1.1);
//fitFunc->SetParameter(2, 0.015);
//fitFunc->SetParLimits(0, 0.24, 0.35);
//fitFunc->SetParLimits(1, 1.0, 1.2);
//fitFunc->SetParLimits(2, 0.005, 0.02);
