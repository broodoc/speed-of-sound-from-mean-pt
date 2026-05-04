#include <TFile.h>
#include <TH1D.h>
#include <TGraphErrors.h>

#include <TFile.h>
#include <TH1D.h>
#include <TGraphErrors.h>

void write_CMS_meanpT_dNdeta() {
    const int nBins = 19;
    double binCenters[nBins] = {1.038, 1.05, 1.06, 1.075, 1.085, 1.095, 1.105, 1.115, 1.126, 1.136, 1.145, 1.155, 1.162, 1.172, 1.18, 1.191, 1.199, 1.208, 1.22};
    double binError[nBins]   = {0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001};  // x error
    double values[nBins]     = {0.9987, 0.9984, 0.9987, 0.999, 0.9996, 1.0005, 1.0015, 1.0028, 1.004, 1.0055, 1.007, 1.0084, 1.0105, 1.012, 1.014, 1.015, 1.017, 1.0185, 1.0197};
    double errors[nBins]     = {0.00025, 0.00025, 0.00025, 0.00025, 0.00025, 0.00025, 0.00025, 0.00025, 0.00025, 0.00025, 0.00025, 0.00025, 0.00025, 0.00025, 0.00025, 0.00025, 0.0005, 0.001, 0.002}; // y error

    // Approximate fixed bin width for TH1D just for visual purposes
    double binWidth = 3.0;
    double binEdges[nBins + 1];
    for (int i = 0; i <= nBins; ++i) {
        if (i == 0)
            binEdges[i] = binCenters[0] - binWidth / 2.0;
        else
            binEdges[i] = binCenters[i - 1] + binWidth / 2.0;
    }

    // 1. Histogram
    TH1D* h = new TH1D("h_dNdeta", "CMS Pb+Pb 5.02 TeV #LTp_{T}#GT vs dN/d#eta;dN/d#eta / dN/d#eta_{0-5%};#LTp_{T}#GT / #LTp_{T}#GT_{0-5%}(GeV/c)", nBins, binEdges);
    for (int i = 0; i < nBins; ++i) {
        h->SetBinContent(i + 1, values[i]);
        h->SetBinError(i + 1, errors[i]);
    }

    // 2. Graph with X errors
    TGraphErrors* g = new TGraphErrors(nBins, binCenters, values, binError, errors);
    g->SetName("g_dNdeta");
    g->SetTitle("CMS Pb+Pb 5.02 TeV #LTp_{T}#GT vs dN/d#eta;dN/d#eta / dN/d#eta_{0-5%};#LTp_{T}#GT / #LTp_{T}#GT_{0-5%}(GeV/c)");
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kBlue+1);

    // 3. Save both
    TFile* fout = new TFile("CMS_5p02TeV_Measurement_arXiv2401.06896.root", "RECREATE");
    h->Write();
    g->Write();
    fout->Close();
    
    std::cout << "Histogram written to CMS_5p02TeV_Measurement_arXiv2401.06896.root" << std::endl;
}


//void write_STAR_meanpT_dNdeta() {
//    // Bin midpoints and values from the CSV
//    const int nBins = 7;
//    double binCenters[nBins] = {30, 100, 197, 285, 420, 558, 691};
//    double binError[nBins] = {0, 12., 20., 25., 45., 65., 82.};
//    double values[nBins]     = {0.45, 0.485, 0.501, 0.514, 0.516, 0.518, 0.518};
//    double errors[nBins]     = {0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02};
//
//    // Create histogram with bin edges (1 bin per value, centered around the binCenters)
//    double binWidth = 3.0;  // Approximated width based on spacing
//    double binEdges[nBins + 1];
//    for (int i = 0; i < nBins + 1; ++i) {
//        if (i == 0) {
//            binEdges[i] = binCenters[0] - binWidth / 2.0;
//        } else {
//            binEdges[i] = binCenters[i - 1] + binWidth / 2.0;
//        }
//    }
//
//    TH1D* h = new TH1D("h_dNdeta", "STAR 200 GeV mean pT vs dN/d#eta;dN/d#eta;<pT> (GeV/c)", nBins, binEdges);
//
//    for (int i = 0; i < nBins; ++i) {
//        h->SetBinContent(i + 1, values[i]);
//        h->SetBinError(i + 1, errors[i]);
//    }
//
//    // Save to ROOT file
//    TFile* fout = new TFile("STAR_200GeV_Measurement_doi_10.1088slash1742-6596slash5slash1slash003.root", "RECREATE");
//    h->Write();
//    fout->Close();
//
//    std::cout << "Histogram written to STAR_200GeV_Measurement_doi_10.1088slash1742-6596slash5slash1slash003.root" << std::endl;
//}

