void view_all_graphs() {
    TFile *f = TFile::Open("graphs_output.root");
    if (!f || f->IsZombie()) return;

    TIter next(f->GetListOfKeys());
    TKey *key;

    while ((key = (TKey*)next())) {
        TObject *obj = key->ReadObj();
        if (obj->InheritsFrom("TGraph")) {
            TCanvas *c = new TCanvas(Form("c_%s", obj->GetName()), obj->GetTitle(), 1400, 1400);
            obj->Draw("APL");
        }
    }
}

