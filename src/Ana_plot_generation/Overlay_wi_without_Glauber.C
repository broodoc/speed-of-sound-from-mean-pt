void Overlay_wi_without_Glauber() {
    // File paths
    TString path_without_Glauber = "/Users/calebbroodo/Desktop/sos_analysis/prepped_files/200COL_Run11/200COL_Run11_rC_w_corrections_acc_from_data_wi_Pip_pT_correction_pTlow0p20_Exclusive/graphs_output_2025-Jul-10_09-03AM.root";
    TString path_wi_Glauber = "/Users/calebbroodo/Desktop/sos_analysis/prepped_files/200COL_Run11/200COL_Run11_rC_w_corrections_acc_from_data_wi_Pip_pT_correction_pTlow0p20_Exclusive/graphs_output_2025-Jul-10_09-14AM.root";

    // Open ROOT files
    TFile* file_excl = TFile::Open(path_without_Glauber);
    TFile* file_incl = TFile::Open(path_wi_Glauber);

    if (!file_excl || !file_incl) {
        std::cerr << "Error opening one of the ROOT files!" << std::endl;
        return;
    }

    // Navigate to the directory
    TString subdirPath = "pT study/Normalized (to ref class) quantities/pT-Integrated Nch/fit_range_0.00-1.10_GeV";
    TDirectory* dir_excl = (TDirectory*)file_excl->Get(subdirPath);
    TDirectory* dir_incl = (TDirectory*)file_incl->Get(subdirPath);

    if (!dir_excl || !dir_incl) {
        std::cerr << "Error finding subdirectory!" << std::endl;
        return;
    }

    TGraphErrors* g_excl = nullptr;
    TGraphErrors* g_incl = nullptr;

//    TString targetName = "0 #leq p_{T} #leq 10 (GeV/c)[LevyTsallis]";
    TString targetName = "0 #leq p_{T} #leq 10 (GeV/c)[Hagedorn]";

    // Loop through all keys in exclusive directory
    TIter next_excl(dir_excl->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next_excl())) {
        if (TString(key->GetName()) == targetName) {
            g_excl = (TGraphErrors*)key->ReadObj();
            break;
        }
    }

    // Loop through all keys in inclusive directory
    TIter next_incl(dir_incl->GetListOfKeys());
    while ((key = (TKey*)next_incl())) {
        if (TString(key->GetName()) == targetName) {
            g_incl = (TGraphErrors*)key->ReadObj();
            break;
        }
    }
    
    if (!g_excl || !g_incl) {
        std::cout << "Exclusive file directory:" << std::endl;
        dir_excl->ls();
        std::cout << "Inclusive file directory:" << std::endl;
        dir_incl->ls();
        std::cerr << "Error retrieving one of the graphs!" << std::endl;
        return;
    }

    // Style the graphs
    g_excl->SetMarkerStyle(20); // filled circle
    g_excl->SetMarkerColor(kBlue);
    g_excl->SetLineColor(kBlue);
    g_excl->SetTitle(";#LTdN/d#eta#GT / #LTdN/d#eta#GT_{0-5%};#LTp_{T}#GT / #LTp_{T}#GT_{0-5%}");

    g_incl->SetMarkerStyle(25); // open square
    g_incl->SetMarkerColor(kRed);
    g_incl->SetLineColor(kRed);

    // Draw
    TCanvas* c1 = new TCanvas("c1", "Overlay Tsallis Fits", 800, 600);
    g_excl->Draw("AP");
    g_incl->Draw("P SAME");

    // Add legend
    TLegend* leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->AddEntry(g_excl, "without Glauber fitting", "pl");
    leg->AddEntry(g_incl, "with Glauber fitting", "pl");
    leg->Draw();

    c1->Update();
}

