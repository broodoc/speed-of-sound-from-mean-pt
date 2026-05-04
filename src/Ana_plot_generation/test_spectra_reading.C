void test_spectra_reading() {
    // Open ROOT file
    TFile* f = TFile::Open("graphs_output_2025-Jul-22_02-41PM.root");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open file!" << std::endl;
        return;
    }

    // Navigate to the relevant directory
    bool cd_success = f->cd("Best Chi2Ndf on spectra/spectra study/fit_range_0.00-1.25_GeV");
    if (!cd_success) {
        std::cerr << "Error: Could not cd into desired directory!" << std::endl;
        return;
    }

    // Fetch the histogram
    TH1D* h = (TH1D*)gDirectory->Get("(ref) 0-5% [LevyTsallis]");
    if (!h) {
        std::cerr << "Error: Histogram not found!" << std::endl;
        return;
    }

    // Debug: print entries
    std::cout << "Successfully loaded histogram: (ref) 0-5% [LevyTsallis]" << std::endl;
    std::cout << "Entries = " << h->GetEntries() << std::endl;

    // Optional: visualize
    TCanvas* c = new TCanvas("c", "Test Canvas", 800, 800);
    h->Draw();
    c->SaveAs("test_ref_0_5_spectra.png");  // Save output if running in batch
}

