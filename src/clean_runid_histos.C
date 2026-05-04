void clean_runid_histos(const char* filename) {
    TFile* f_in = TFile::Open(filename, "READ");
    if (!f_in || f_in->IsZombie()) {
        std::cerr << "Could not open file: " << filename << std::endl;
        return;
    }

    TString outname = TString("cleaned_") + filename;
    TFile* f_out = new TFile(outname, "RECREATE");

    TIter nextDir(f_in->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)nextDir())) {
        TObject* obj = key->ReadObj();
        if (obj->InheritsFrom("TDirectory")) {
            TDirectory* dir_in = (TDirectory*)obj;
            TDirectory* dir_out = f_out->mkdir(dir_in->GetName());
            dir_out->cd();

            TIter nextHist(dir_in->GetListOfKeys());
            TKey* hkey;
            while ((hkey = (TKey*)nextHist())) {
                TObject* hist = hkey->ReadObj();
                TString name = hist->GetName();
                if (!name.Contains("hFXTMult_RunID") &&
                    !name.EqualTo("Centrality__pTlow0p20_hnBTofMatch_RunID") &&
                    !name.EqualTo("dEta1p0__pTlow0p20_hnBTofMatch_RunID")) {
                    hist->Write();
                } else {
                    std::cout << "Skipping: " << name << std::endl;
                }
            }
            f_out->cd();
        } else {
            TString name = obj->GetName();
            if (!name.Contains("hFXTMult_RunID") &&
                !name.EqualTo("Centrality__pTlow0p20_hnBTofMatch_RunID") &&
                !name.EqualTo("dEta1p0__pTlow0p20_hnBTofMatch_RunID")) {
                obj->Write();
            } else {
                std::cout << "Skipping: " << name << std::endl;
            }
        }
    }

    f_out->Close();
    f_in->Close();
}
