void write_STAR_200GeV_dNdeta() {
    // Bin midpoints and values from the CSV
    const int nBins = 11;
    double binCenters[nBins] = {1.5, 4.5, 8.0, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5};
    double values[nBins]     = {700, 629, 548, 455, 376, 312, 252, 202, 164, 130, 95};
    double errors[nBins]     = {27,  24,  21,  18,  15,  12,  10,   8,   6,   5,   4};

    // Create histogram with bin edges (1 bin per value, centered around the binCenters)
    double binWidth = 3.0;  // Approximated width based on spacing
    double binEdges[nBins + 1];
    for (int i = 0; i < nBins + 1; ++i) {
        if (i == 0) {
            binEdges[i] = binCenters[0] - binWidth / 2.0;
        } else {
            binEdges[i] = binCenters[i - 1] + binWidth / 2.0;
        }
    }

    TH1D* h = new TH1D("h_dNdeta", "STAR 200 GeV dN/d#eta;Centrality Bin (%);dN/d#eta", nBins, binEdges);

    for (int i = 0; i < nBins; ++i) {
        h->SetBinContent(i + 1, values[i]);
        h->SetBinError(i + 1, errors[i]);
    }

    // Save to ROOT file
    TFile* fout = new TFile("STAR_200GeV_Measurement_arXiv_nucl_ex_0201005v2.root", "RECREATE");
    h->Write();
    fout->Close();

    std::cout << "Histogram written to STAR_200GeV_Measurement_arXiv_nucl_ex_0201005v2.root" << std::endl;
}

