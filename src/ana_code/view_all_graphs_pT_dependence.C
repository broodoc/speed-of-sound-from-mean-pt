#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TString.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TLatex.h>
#include <iostream>
#include <TBox.h>
#include <TPad.h>
#include <TLine.h>
#include <TPave.h>

void SetPaperStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // Font: 42 = Helvetica-like (clean). Use 132 for Times-like.
  gStyle->SetTextFont(42);
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");

  gStyle->SetTitleSize(0.045, "XYZ");
  gStyle->SetLabelSize(0.040, "XYZ");

  gStyle->SetTitleOffset(1.25, "X");
  gStyle->SetTitleOffset(1.60, "Y");

  gStyle->SetEndErrorSize(4);      // nicer error bar caps
  gStyle->SetLineWidth(4);
  gStyle->SetFrameLineWidth(2);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}

static TString SanitizePtLabel(TString s) {
    s.ReplaceAll("p_{T}", "pT");
    s.ReplaceAll("#leq", "leq");
    s.ReplaceAll("(GeV/c)", "GeV");
    s.ReplaceAll("GeV/c", "GeV");
    s.ReplaceAll("#",""); s.ReplaceAll("{",""); s.ReplaceAll("}","");
    s.ReplaceAll("(",""); s.ReplaceAll(")",""); s.ReplaceAll("/","_");
    s.ReplaceAll(" ", "_");
    // trim common fit-function suffix, optional
    s.ReplaceAll("_cLevyTsallis", "");
    return s;
}


void CopyGraphStyle(TGraph* from, TGraph* to) {
    to->SetMarkerColor(from->GetMarkerColor());
    to->SetMarkerStyle(from->GetMarkerStyle());
    to->SetMarkerSize(from->GetMarkerSize());
    to->SetFillColor(from->GetFillColor());
    to->SetFillStyle(from->GetFillStyle());

    if (from->GetXaxis()) {
        to->GetXaxis()->SetTitle(from->GetXaxis()->GetTitle());
        to->GetXaxis()->SetTitleSize(from->GetXaxis()->GetTitleSize());
        to->GetXaxis()->SetTitleOffset(from->GetXaxis()->GetTitleOffset());
        to->GetXaxis()->SetLabelSize(from->GetXaxis()->GetLabelSize());
    }
    if (from->GetYaxis()) {
        to->GetYaxis()->SetTitle(from->GetYaxis()->GetTitle());
        to->GetYaxis()->SetTitleSize(from->GetYaxis()->GetTitleSize());
        to->GetYaxis()->SetTitleOffset(from->GetYaxis()->GetTitleOffset());
        to->GetYaxis()->SetLabelSize(from->GetYaxis()->GetLabelSize());
    }
}

std::pair<double, double> parseCentrality(const TString& name) {
    TString s = name;
    s.ReplaceAll("(", "");
    s.ReplaceAll(")", "");
    s.ReplaceAll("ref", "0.000-5.000"); // (ref) 0-5% becomes 0-5

    double low = -1, high = -1;
    if (s.Contains("-") && s.Contains("%")) {
        Ssiz_t start = s.Index("0.");
        if (start < 0) start = s.Index("0-");
        TString range = s(start, s.Length());
        range.ReplaceAll("%", "");

        TObjArray* tokens = range.Tokenize("-");
        if (tokens->GetEntries() == 2) {
            low  = ((TObjString*)tokens->At(0))->GetString().Atof();
            high = ((TObjString*)tokens->At(1))->GetString().Atof();
        }
        delete tokens;
    }
    return {low, high};
}

void drawCentralitySlider(TCanvas* c, double low, double high) {
    double sliderXmin = 0.55;
    double sliderXmax = 0.85;
    double sliderY = 0.8;
    double sliderH = 0.02;

    // Create overlay pad using normalized coordinates
    c->cd();
    TPad* overlay = new TPad("overlay", "", 0, 0, 1, 1);
    overlay->SetFillStyle(0);
    overlay->SetFrameFillStyle(0);
    overlay->SetBit(kCanDelete); // Clean up on canvas deletion
    overlay->Draw();
    overlay->cd();

    // Draw base line (slider track)
    TLine* line = new TLine(sliderXmin, sliderY, sliderXmax, sliderY);
    line->SetNDC(true);
    line->SetLineWidth(8);
    line->Draw();

    // Draw centrality slider box (no SetNDC needed in overlay)
    double binXmin = pow(sliderXmin + ((100 - high) / 100.0) * (sliderXmax - sliderXmin) , 4) * 1.5;
    if (binXmin > sliderXmax - 0.04) binXmin = sliderXmax - 0.04;
    if (binXmin < sliderXmin) binXmin = sliderXmin;
    double binXmax = binXmin + 0.04;
//    double binXmax = sliderXmin + ((100 - high) / 100.0) * (sliderXmax - sliderXmin) + 0.02;
//    double binXmin = sliderXmin + (high / 100.0) * (sliderXmax - sliderXmin);
//    double binXmax = sliderXmin + (low / 100.0) * (sliderXmax - sliderXmin);
    

//    std::cout << "Drawing TBox with low = " << low << ", high = " << high << ", binXmin = " << binXmin << ", binXmax = " << binXmax << std::endl;
    TPave* box = new TPave(binXmin, sliderY - sliderH / 2, binXmax, sliderY + sliderH / 2, 0, "brNDC");
    box->SetFillColor(kBlack);
    box->SetLineColor(kYellow);
    box->SetCornerRadius(0.5);
    box->Draw();
//    TBox* box = new TBox(binXmin, sliderY - sliderH / 2, binXmax, sliderY + sliderH / 2);
//    box->SetFillColor(kBlack);
//    box->SetLineColor(kBlack);
////    box->SetCornerRadius(0.35);
//    box->Draw();

    // Labels
    TLatex* left = new TLatex(sliderXmin, sliderY + 0.04, "more");
    TLatex* left_down = new TLatex(sliderXmin, sliderY + 0.02, "peripheral");
    TLatex* right = new TLatex(sliderXmax - 0.05, sliderY + 0.04, "more");
    TLatex* right_down = new TLatex(sliderXmax - 0.05, sliderY + 0.02, "central");
    left->SetNDC(); right->SetNDC(); left_down->SetNDC(); right_down->SetNDC();
    left->SetTextSize(0.03);
    right->SetTextSize(0.03);
    left_down->SetTextSize(0.03);
    right_down->SetTextSize(0.03);
    left->Draw(); right->Draw();
    left_down->Draw(); right_down->Draw();

    // Ultra-central label
//    if (low < 0.05) {
//        TLatex* ultra = new TLatex(sliderXmin, sliderY + 0.04, "Ultra-central class");
//        ultra->SetNDC();
//        ultra->SetTextSize(0.035);
//        ultra->SetTextColor(kRed + 1);
//        ultra->Draw();
//    }

    c->cd(); // Back to base pad
}

void format_ptstudy_norm(TCanvas* c, TObject* obj, const TString& name, const TString& outputDirPath) {
    c->SetTickx(); c->SetTicky();
    gStyle->SetOptStat(0);
    
    SetPaperStyle();
    c->SetLeftMargin(0.14);
    c->SetRightMargin(0.04);
    c->SetTopMargin(0.06);
    c->SetBottomMargin(0.12);
    c->SetTitle("");
    
    const int nFrames = 15;

    double x_mins[nFrames];
    double y_mins[nFrames];

    double x_min = 0.0;
    double x_max = 0.7;
    double y_min = 0.8;
    double y_max = 0.96;

    double dx = (x_max - x_min) / (nFrames - 1);
    double dy = (y_max - y_min) / (nFrames - 1);

    for (int i = 0; i < nFrames; ++i) {
        x_mins[i] = x_min + i * dx;
        y_mins[i] = y_min + i * dy;
    }

    
    if (obj->InheritsFrom("TGraph")) {
        TGraph* g = (TGraph*)obj;
        TGraphErrors* graph = new TGraphErrors();
        CopyGraphStyle(g, graph);
        int j = 0;

        for (int i = 0; i < g->GetN(); ++i) {
            double x, y;
            g->GetPoint(i, x, y);

            double ex = g->GetErrorX(i);
            double ey = g->GetErrorY(i);

            // Keep only nonzero points
            if (x == 0 && y == 0 && ex == 0 && ey == 0) continue;

            graph->SetPoint(j, x, y);
            graph->SetPointError(j, ex, ey);
            ++j;
        }
        
        c->cd();
        graph->Draw("AP");
        
        
        graph->SetTitle("");  // remove internal title
        graph->GetYaxis()->SetTitleOffset(1.8);  // more room
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(2.5);
        int deepGreen = TColor::GetColor(0, 120, 60);
        graph->SetMarkerColor(deepGreen);
        graph->SetLineColor(deepGreen);
//        graph->SetMarkerColor(kBlack);
//        graph->SetLineColor(kBlack);     // error bars
        graph->SetLineWidth(5);
//
//        graph->SetMarkerStyle(20);
//        graph->SetMarkerSize(1.0);
//        graph->SetMarkerColor(kGreen + 1);
        
        for (int i = 0; i < nFrames; ++i) {
            if (i == (nFrames - 1)){
                g->SetTitle("");  // remove internal title
                g->GetYaxis()->SetTitleOffset(1.8);  // more room
                g->SetMarkerStyle(20);
                g->SetMarkerSize(1.0);
                g->SetMarkerColor(kGreen + 1);
                g->GetXaxis()->SetRangeUser(x_mins[i], 1.25);
                g->GetYaxis()->SetRangeUser(y_mins[i], 1.06);
                c->cd();
                g->Draw("P E");
                c->Update();
                c->SaveAs(Form("%s/ptstudy_normalized_zoom_frame_%02d.png", outputDirPath.Data(), i+1));
                continue;
            }
            
            graph->GetXaxis()->SetRangeUser(x_mins[i], 1.25);
            graph->GetYaxis()->SetRangeUser(y_mins[i], 1.06);
            
            c->cd();
            graph->Draw("P E");
            c->Update();
            c->SaveAs(Form("%s/ptstudy_normalized_zoom_frame_%02d.png", outputDirPath.Data(), i+1));
        }
        
    } else if (obj->InheritsFrom("TH1")) {
        TH1* hist = (TH1*)obj;
        hist->SetTitle("");  // remove internal title
        hist->GetYaxis()->SetTitleOffset(1.8);
        hist->SetMarkerStyle(20);
        hist->SetMarkerSize(1.5);

        for (int i = 0; i < nFrames; ++i) {
            hist->GetXaxis()->SetRangeUser(x_mins[i], 1.25);
            hist->GetYaxis()->SetRangeUser(y_mins[i], 1.06);

            c->cd();
            hist->Draw("P E");
            c->SaveAs(Form("%s/ptstudy_normalized_zoom_frame_%02d.png", outputDirPath.Data(), i+1));
        }
    }
    

}

void format_ptstudy_absolute(TCanvas* c, TObject* obj, const TString& name, const TString& outputDirPath) {
    c->SetTickx(); c->SetTicky();
    gStyle->SetOptStat(0);
    c->SetTitle("");
    c->SetRightMargin(0.05);
    c->SetLeftMargin(0.15);
    
    const int nFrames = 20;

    double x_mins[nFrames] = {
        0.000, 34.211, 68.421, 102.632, 136.842,
        171.053, 205.263, 239.474, 273.684, 307.895,
        342.105, 376.316, 410.526, 444.737, 478.947,
        513.158, 547.368, 581.579, 615.789, 650.000
    };

    double y_mins[nFrames] = {
        0.4400, 0.4463, 0.4526, 0.4589, 0.4653,
        0.4716, 0.4779, 0.4842, 0.4905, 0.4968,
        0.5032, 0.5095, 0.5158, 0.5221, 0.5284,
        0.5347, 0.5411, 0.5474, 0.5537, 0.5600
    };
    
    if (obj->InheritsFrom("TGraph")) {
        TGraph* g = (TGraph*)obj;
        TGraphErrors* graph = new TGraphErrors();
        CopyGraphStyle(g, graph);
        int j = 0;

        for (int i = 0; i < g->GetN(); ++i) {
            double x, y;
            g->GetPoint(i, x, y);

            double ex = g->GetErrorX(i);
            double ey = g->GetErrorY(i);

            // Keep only nonzero points
            if (x == 0 && y == 0 && ex == 0 && ey == 0) continue;

            graph->SetPoint(j, x, y);
            graph->SetPointError(j, ex, ey);
            ++j;
        }
        
        c->cd();
        graph->Draw("AP");
        
        
        graph->SetTitle("");  // remove internal title
        graph->GetYaxis()->SetTitleOffset(1.8);  // more room
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(1.5);
        graph->SetMarkerColor(kGreen + 2);
        
        for (int i = 0; i < nFrames; ++i) {
            if (i == (nFrames - 1)){
                g->SetTitle("");  // remove internal title
                g->GetYaxis()->SetTitleOffset(1.8);  // more room
                g->SetMarkerStyle(20);
                g->SetMarkerSize(1.0);
                g->SetMarkerColor(kGreen + 2);
                g->GetXaxis()->SetRangeUser(x_mins[i], 800);
                g->GetYaxis()->SetRangeUser(y_mins[i], 0.58);
                // --- Draw transparent red band for T_eff region ---
                TBox* shadeBox = new TBox(x_mins[i], 0.569, 800.0, 0.576);
                shadeBox->SetFillColorAlpha(kRed, 0.2);  // light transparent red
                shadeBox->SetLineColor(0);  // no border
                shadeBox->Draw("same");

                // --- Draw annotation label ---
                TLatex* tex = new TLatex();
                tex->SetNDC(false);
                tex->SetTextSize(0.03);
                tex->SetTextFont(42);
                tex->DrawLatex(x_mins[i] + 10, 0.5775, "T_{eff} #in (189.7, 192.0) MeV");
                c->cd();
                g->Draw("P E");
                c->Update();
                c->SaveAs(Form("%s/ptstudy_absolute_zoom_frame_%02d.png", outputDirPath.Data(), i+1));
                continue;
            }
            
            graph->GetXaxis()->SetRangeUser(x_mins[i], 800);
            graph->GetYaxis()->SetRangeUser(y_mins[i], 0.58);
            
            c->cd();
            graph->Draw("P E");
            c->Update();
            c->SaveAs(Form("%s/ptstudy_absolute_zoom_frame_%02d.png", outputDirPath.Data(), i+1));
        }
        
    } else if (obj->InheritsFrom("TH1")) {
        TH1* hist = (TH1*)obj;
        hist->SetTitle("");  // remove internal title
        hist->GetYaxis()->SetTitleOffset(1.8);
        hist->SetMarkerStyle(20);
        hist->SetMarkerSize(1.5);

        for (int i = 0; i < nFrames; ++i) {
            hist->GetXaxis()->SetRangeUser(x_mins[i], 1.25);
            hist->GetYaxis()->SetRangeUser(y_mins[i], 1.06);

            c->cd();
            hist->Draw("P E");
            c->SaveAs(Form("%s/ptstudy_absolute_zoom_frame_%02d.png", outputDirPath.Data(), i+1));
        }
    }
    

}

void format_spectra(TCanvas* c, TObject* obj, const TString& name) {
    c->SetLogy(); c->SetGridy(); c->SetTickx(); c->SetTicky();
    gStyle->SetOptStat(0);

    if (obj->InheritsFrom("TGraph")) {
        ((TGraph*)obj)->GetYaxis()->SetTitleOffset(1.2);
        ((TGraph*)obj)->GetYaxis()->SetRangeUser(0.001, 10000);
        ((TGraph*)obj)->GetXaxis()->SetRangeUser(0, 5);
        ((TGraph*)obj)->SetMarkerStyle(20);
        ((TGraph*)obj)->SetMarkerSize(3);
    } else if (obj->InheritsFrom("TH1")) {
        ((TH1*)obj)->GetYaxis()->SetTitleOffset(1.2);
        ((TH1*)obj)->GetYaxis()->SetRangeUser(0.001, 10000);
        ((TH1*)obj)->GetXaxis()->SetRangeUser(0, 5);
        ((TH1*)obj)->SetMarkerStyle(20);
        ((TH1*)obj)->SetMarkerSize(1);
    }

    auto [low, high] = parseCentrality(name);
    if (low >= 0 && high > 0) {
        drawCentralitySlider(c, low, high);
    }

//    TPaveText* pave = new TPaveText(0.75, 4.75, 0.95, 4.95, "NDC");
//    pave->SetFillColor(0);
//    pave->SetBorderSize(0);
//    pave->AddText(name.Data());
//    pave->Draw();
}

void format_centrality_study_dNdeta(TCanvas* c, TObject* obj) {
    c->SetLogx();
    c->SetGridy();
    c->SetTickx();
    c->SetTicky();

    if (obj->InheritsFrom("TGraph")) {
        ((TGraph*)obj)->GetYaxis()->SetRangeUser(0, 800);
        ((TGraph*)obj)->GetYaxis()->SetTitleOffset(1.2);
        ((TGraph*)obj)->GetXaxis()->SetRangeUser(0.0001, 100);
    } else if (obj->InheritsFrom("TH1")) {
        ((TH1*)obj)->GetYaxis()->SetRangeUser(0, 800);
        ((TH1*)obj)->GetXaxis()->SetRangeUser(0.0001, 100);
        ((TH1*)obj)->GetYaxis()->SetTitleOffset(1.2);
    }
}

//void processDirectory(TDirectory* dir, const TString& outputDirPath, const TString& pathPrefix = "") {
void processDirectory(TDirectory* dir,
                      const TString& outputDirPath,
                      const TString& pathPrefix = "",
                      TString quantityClass = "",
                      TString currentFitRange = "",
                      TString currentPtRange = ""){
    gErrorIgnoreLevel = kPrint;
    TIter next(dir->GetListOfKeys());
    TKey* key;

    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();

        // Recurse into subdirectories
        TDirectory* subDir = dynamic_cast<TDirectory*>(obj);
        if (subDir) {
            TString newPrefix = pathPrefix + "_" + subDir->GetName();

            // Update the quantity class if we’re entering that node
            TString newQuantityClass = quantityClass;
            if (TString(subDir->GetName()).Contains("Absolute quantities"))
                newQuantityClass = "Absolute_quantities";
            else if (TString(subDir->GetName()).Contains("Normalized"))
                newQuantityClass = "Normalized_to_ref_class_quantities";

            // Track the current fit-range directory name verbatim
            TString newFit = currentFitRange;
            if (TString(subDir->GetName()).BeginsWith("fit_range_"))
                newFit = subDir->GetName();  // e.g. "fit_range_0.00-1.1_GeV"
            
            // pT range (anything that looks like "... p_{T} ...")
            TString newPt = currentPtRange;
            TString subName = subDir->GetName();
            if (subName.Contains("p_{T}")) {
                // sanitize to a filesystem-friendly name
                TString safePt = subName;
                safePt.ReplaceAll("p_{T}", "pT");
                safePt.ReplaceAll("(GeV/c)", "GeV");
                safePt.ReplaceAll("#", "");
                safePt.ReplaceAll("{", "");
                safePt.ReplaceAll("}", "");
                safePt.ReplaceAll("(", "");
                safePt.ReplaceAll(")", "");
                safePt.ReplaceAll("/", "_");
                safePt.ReplaceAll(" ", "_");
                newPt = safePt;  // e.g. 0_leq_pT_leq_10_GeV
            }
            
            processDirectory(subDir, outputDirPath, newPrefix, newQuantityClass, newFit, newPt);
            continue;
        }
        
        // Recurse into subdirectories
//        TDirectory* subDir = dynamic_cast<TDirectory*>(obj);
//        if (subDir) {
//            TString newPrefix = pathPrefix + "_" + subDir->GetName();
//            processDirectory(subDir, outputDirPath, newPrefix);
//            continue;
//        }

        // Now go back to the original directory
        dir->cd();  // Ensures we're in the correct scope
//        const char* graphName = key->GetName();
        TString graphName = key->GetName();
        TObject* fetched = dir->Get(graphName);

        if (!obj) {
            std::cout << "Could not read object from key: " << key->GetName() << std::endl;
            continue;
        }
        
//        if (!fetched) {
//            std::cout << "Could not fetch " << graphName << "!" << std::endl;
//            continue;
//        }

        TString rawName = pathPrefix + "__" + graphName;

        rawName.ReplaceAll("#", "");
        rawName.ReplaceAll("<", "lt");
        rawName.ReplaceAll(">", "gt");
        rawName.ReplaceAll("(", "");
        rawName.ReplaceAll(")", "");
        rawName.ReplaceAll("[", "");
        rawName.ReplaceAll("]", "");
        rawName.ReplaceAll("/", "_");
        rawName.ReplaceAll(" ", "_");
        
        // Derive a pT folder from the leaf name if we didn't descend into a pT subdir
        TString leafPt = currentPtRange;
        if (leafPt.IsNull()) {
            // Prefer the key name; it’s short (e.g. "0 #leq p_{T} #leq 10 (GeV/c)")
            if (graphName.Contains("p_{T}")) {
                leafPt = SanitizePtLabel(graphName);
            } else if (rawName.Contains("p_{T}")) {
                // fallback: sometimes the generated name carries the pT clause
                leafPt = SanitizePtLabel(rawName);
            }
        }

        TString localOutDir = outputDirPath;
        if (quantityClass.Length())   localOutDir += "/" + quantityClass;      // e.g. Absolute_quantities
        if (currentFitRange.Length()) localOutDir += "/" + currentFitRange;    // e.g. fit_range_0.00-1.1_GeV
        if (leafPt.Length())          localOutDir += "/" + leafPt;
//        if (currentPtRange.Length())  localOutDir += "/" + currentPtRange;
        gSystem->Exec(Form("mkdir -p \"%s\"", localOutDir.Data()));
        
        TCanvas* c;
        if (fetched) {
            c = new TCanvas(Form("c_%s", rawName.Data()), fetched->GetTitle(), 1400, 1400);
            TCanvas* originalCanvas = dynamic_cast<TCanvas*>(fetched);
            if (originalCanvas) {
                originalCanvas->DrawClonePad();
            } else if (fetched->InheritsFrom("TGraph")) {
                fetched->Draw("APL");
            } else if (fetched->InheritsFrom("TH1")) {
                fetched->Draw();
            }
            if (rawName.Contains("spectra_study")) {
//                continue;
                format_spectra(c, fetched, graphName);
            } else if (graphName.Contains("dNdEta")) {
                format_centrality_study_dNdeta(c, fetched);
            } else if (rawName.Contains("pT_study_Normalized_to_ref_class_quantities")) {
                continue;
                format_ptstudy_norm(c, fetched, graphName, localOutDir);
            }
                
        }
        else if (obj){
            c = new TCanvas(Form("c_%s", rawName.Data()), obj->GetTitle(), 1400, 1400);

            TCanvas* originalCanvas = dynamic_cast<TCanvas*>(obj);
            
//            if (originalCanvas) {
//                        originalCanvas->DrawClonePad();
//                    } else if (obj->InheritsFrom("TGraph")) {
//                        obj->Draw("APL");
//                    } else if (obj->InheritsFrom("TH1")) {
//                        obj->Draw();
//                    }
            
            if (rawName.Contains("spectra_study")) {
                format_spectra(c, obj, graphName);
            } else if (graphName.Contains("dNdEta")) {
                format_centrality_study_dNdeta(c, obj);
            } else if (rawName.Contains("pT_study_Normalized_to_ref_class_quantities")) {
                format_ptstudy_norm(c, obj, graphName, localOutDir);
            } else if (rawName.Contains("Absolute")) {
                format_ptstudy_absolute(c, obj, graphName, localOutDir);
            }

        }
        TString savePath = Form("%s/%s.png", localOutDir.Data(), rawName.Data());
        c->SaveAs(savePath);

//        TString savePath = Form("%s/%s.png", outputDirPath.Data(), rawName.Data());
//        c->SaveAs(savePath);
        delete c;
    }
}


void view_all_graphs_pT_dependence(const char* rootFilePath, const char* outputDirPath) {
    gSystem->Exec(Form("mkdir -p \"%s\"", outputDirPath));

    TFile* file = TFile::Open(rootFilePath);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: could not open file " << rootFilePath << std::endl;
        return;
    }

    TDirectory* spectraDir = (TDirectory*)file->Get("Spectra fit ranges results");
    if (!spectraDir) {
        std::cerr << "Error: 'Spectra directory on spectra' directory not found.\n";
        return;
    }

    processDirectory(spectraDir, outputDirPath, "spectrafit_rangesresults");
    file->Close();
}
