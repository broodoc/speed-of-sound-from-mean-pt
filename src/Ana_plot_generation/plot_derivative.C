#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TDirectory.h>
#include <iostream>
#include <vector>

void plot_derivative() {
    const char* filepath = "/Users/calebbroodo/Desktop/sos_analysis/prepped_files/200COL_Run11/200COL_Run11_rC_w_corrections_acc_from_data_pT_from_Pip_correction_Exclusive/graphs_output_2025-Jul-23_09-31PM.root";
    
    TFile* file = TFile::Open(filepath);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file.\n";
        return;
    }

    // Navigate into the "general" directory
    TDirectory* generalDir = (TDirectory*)file->Get("general");
    if (!generalDir) {
        std::cerr << "Error: general directory not found in file.\n";
        return;
    }

    // Get the histogram from the general directory
    TH1D* h = nullptr;
    generalDir->GetObject("hRefMultCorrGlauber", h);
    if (!h) {
        std::cerr << "Error: hRefMultCorrGlauber not found inside general directory.\n";
        return;
    }

    const int nBins = h->GetNbinsX();
    std::vector<double> x_vals, dydx_vals;

    for (int i = 2; i < nBins; ++i) { // central difference: skip first and last bin
        double x_prev = h->GetBinCenter(i - 1);
        double y_prev = h->GetBinContent(i - 1);

        double x_next = h->GetBinCenter(i + 1);
        double y_next = h->GetBinContent(i + 1);

        double dx = x_next - x_prev;
        if (dx == 0) continue; // skip if spacing is invalid

        double derivative = (y_next - y_prev) / dx;

        double x_mid = h->GetBinCenter(i); // position of the derivative
        x_vals.push_back(x_mid);
        dydx_vals.push_back(derivative);
    }

    // Plot with TGraph
    TGraph* gDeriv = new TGraph(x_vals.size(), &x_vals[0], &dydx_vals[0]);
    gDeriv->SetTitle("d(EventCount)/d(N_{ch});N_{ch};d(EventCount)/d(N_{ch})");
    gDeriv->SetLineColor(kRed + 1);
    gDeriv->SetLineWidth(2);
    gDeriv->SetMarkerStyle(20);
    gDeriv->SetMarkerSize(0.8);

    TCanvas* c1 = new TCanvas("c1", "First Derivative Plot", 900, 600);
    gDeriv->Draw("APL");

    c1->SaveAs("hRefMultCorrGlauber_derivative.png");
}
