#include <TFile.h>
#include <TH1D.h>
#include <TGraphErrors.h>

#include <TFile.h>
#include <TH1D.h>
#include <TGraphErrors.h>

void write_STAR_meanpT_dNdeta() {
    const int nBins = 7;
    double binCenters[nBins] = {30, 100, 197, 285, 420, 558, 691};
    double binError[nBins]   = {0, 12., 20., 25., 45., 65., 82.};  // x error
    double values[nBins]     = {0.45, 0.485, 0.501, 0.514, 0.516, 0.518, 0.522};
    double errors[nBins]     = {0.02, 0.02, 0.02, 0.02, 0.02, 0.025, 0.045}; // y error

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
    TH1D* h = new TH1D("h_dNdeta", "STAR 200 GeV #LTp_{T}#GT vs dN/d#eta;dN/d#eta;#LTp_{T}#GT (GeV/c)", nBins, binEdges);
    for (int i = 0; i < nBins; ++i) {
        h->SetBinContent(i + 1, values[i]);
        h->SetBinError(i + 1, errors[i]);
    }

    // 2. Graph with X errors
    TGraphErrors* g = new TGraphErrors(nBins, binCenters, values, binError, errors);
    g->SetName("g_dNdeta");
    g->SetTitle("STAR 200 GeV #LTp_{T}#GT vs dN/d#eta;dN/d#eta;#LTp_{T}#GT (GeV/c)");
    g->SetMarkerStyle(20);
    g->SetMarkerColor(kBlue+1);

    // 3. Save both
    TFile* fout = new TFile("STAR_200GeV_Measurement_doi_10.1088slash1742-6596slash5slash1slash003.root", "RECREATE");
    h->Write();
    g->Write();
    fout->Close();
    
    std::cout << "Histogram written to STAR_200GeV_Measurement_doi_10.1088slash1742-6596slash5slash1slash003.root" << std::endl;
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

