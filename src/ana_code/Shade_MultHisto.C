#include "Parameterization.h"
#include "Shade_MultHisto.h"
#include <fstream>

void Shade_MultHisto(string dir,string prepdir,string refMult, bool debug = false){
    string system = "AuAu_200";
    gErrorIgnoreLevel = kFatal;
    Parameterization per{};
    per.SetSystem(system);
    per.SetMultEstimator("TPC");
    
    std::string filename;
    
    for (const auto& entry : fs::directory_iterator(prepdir)) {
        if (entry.is_regular_file()) {
            std::string prepname = entry.path().filename().string();
            if (prepname.find("bin_cuts.txt") != std::string::npos) {
                filename = entry.path().string();
                break; // stop at first match
            }
        }
    }

    if (filename.empty()) {
        std::cerr << "❌ No *bin_cuts.txt file found in: " << prepdir << std::endl;
        return;
    }

    std::cout << "✅ Found bin cuts file: " << filename << std::endl;
    
        std::ifstream infile(filename);
        
        if (!infile.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            std::exit(EXIT_FAILURE);
        }
        
        std::string line;
        std::vector<int> BinCuts;
        
        if (std::getline(infile, line)) {
            std::stringstream ss(line);
            std::string item;
            
            while (std::getline(ss, item, ',')) {
                // Remove any leading or trailing whitespace from the item
                item.erase(0, item.find_first_not_of(" \t\n\r"));
                item.erase(item.find_last_not_of(" \t\n\r") + 1);
                
                if (!item.empty()) {
                    // Convert the item to an integer and add to the vector
                    BinCuts.push_back(std::stoi(item));
                }
            }
        }
        infile.close();
        
        // Print the numbers to verify
        int i = 0;
        cout << "Bin cuts: ";
        for (const auto& num : BinCuts) {
            std::cout << num << ", ";
            i++;
        }
        cout << endl;
        
        
        std::cout << std::endl;
        per.SetFineBinCutsTPC(BinCuts.data());
        if (debug) per.PrintFineBinCuts();
    
    int nperc_bins{per.GetnFinePercentileBinning()};
    
    if (debug) {
        cout << "\t Low Bin percentiles: ";
        for (int i = 0; i < nperc_bins+1; i++){
            cout << per.GetFineAmplitudeCutLow(i)*100 << "%(" << per.GetFineAmplitudeBinCutLow(i) << "), ";
        }
        cout << '\n' << "\t High Bin percentiles: ";
        for (int j = 0; j < nperc_bins; j++){
            cout << per.GetFineAmplitudeCutHigh(j)*100 << "%(" << per.GetFineAmplitudeBinCutHigh(j) << "), ";
        }
    cout << '\n' << '\n';
    }
    
        std::string fullFilename = dir + "temp";
//        std::ifstream infile(fullFilename);
        cout << "\nfullFilename: " << fullFilename << endl;
        TFile *file = TFile::Open(fullFilename.c_str(), "READ");
        
//        TH2D* hPt_Nch =  (TH2D*)file->Get(pTMult.c_str());
        TH1D* Nch_dist = (TH1D*)file->Get(refMult.c_str());
        
        if (!Nch_dist) {
            std::cerr << "Error: Missing histograms in file " << fullFilename << std::endl;
//            if (!hPt_Nch) std::cerr << "  → Missing: " << pTMult << std::endl;
            if (!Nch_dist) std::cerr << "  → Missing: " << refMult << std::endl;
        }
    
//    Form("#%d: %.5f-%.5f%%",
//                                          m_index,
//                                          per.GetFineAmplitudeCutLow(k)*100,
//                                          per.GetFineAmplitudeCutHigh(0)*100)
        
        double integral = Nch_dist->Integral();
        if (integral > 0) {Nch_dist->Scale(1.0 / integral);}
        else {std::cerr << "Warning: Histogram has zero integral!" << std::endl;}
        
//        if (refMult.find("refMult2") != std::string::npos) {
//            Nch_dist->GetYaxis()->SetTitle("P(refMult2) (a.u.)");
//            Nch_dist->GetXaxis()->SetTitle("refMult2");
//        } else {
//            Nch_dist->GetYaxis()->SetTitle("P(refMult) (a.u.)");
//            Nch_dist->GetXaxis()->SetTitle("refMult");
//        }
    
//        Nch_dist->GetYaxis()->SetTitle("P(refMult2) (a.u.)");
//        Nch_dist->GetXaxis()->SetTitle("refMult2");
    
    Nch_dist->GetYaxis()->SetTitle("P(refMult) (a.u.)");
    Nch_dist->GetXaxis()->SetTitle("refMult");
      
        Nch_dist->SetTitle("");
        
        gROOT->SetBatch(kTRUE);
        gStyle->SetOptStat(0);

    
        TCanvas* c1 = new TCanvas("c1", "pdf shading", 800, 800);
        c1->SetLogy();
        c1->SetTicks(1, 1);
        Nch_dist->SetMinimum(1e-7);
        Nch_dist->SetMaximum(1e-1);
        Nch_dist->GetYaxis()->SetTitleSize(0.035);
        Nch_dist->GetYaxis()->SetTitleOffset(1.25);

        Nch_dist->Draw("HIST");
        
        std::cout << "Adding shaded event classes" << std::endl;
        
        const int nCustomColors = sizeof(customColors) / sizeof(int);
    
        TLegend* leg = new TLegend(0.12, 0.7, 0.89, 0.88);
        leg->SetNColumns(4);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.015);
    
    
        for (int k = 0; k < nperc_bins; k++) {
            int low_cut = per.GetFineAmplitudeBinCutLow(k);
            if (low_cut < 0) {low_cut = 1;}
            int high_cut = per.GetFineAmplitudeBinCutHigh(k);
            if (k == 0) { high_cut = Nch_dist->GetNbinsX(); }

            std::cout << "  Bin " << k << " [" << per.GetFineAmplitudeCutLow(k)*100 << "% - " << per.GetFineAmplitudeCutHigh(k)*100 << "%]: " << low_cut << " to " << high_cut << std::endl;
            
//            if (low_cut <= 0 || high_cut <= 0 || high_cut >= Nch_dist->GetNbinsX() || low_cut > high_cut) continue;

            double low_edge = Nch_dist->GetBinLowEdge(low_cut);
            double high_edge = Nch_dist->GetBinLowEdge(high_cut + 1);  // Get upper edge of high_cut bin

            double maxY = Nch_dist->GetMaximum();
            
            TH1D* shaded = (TH1D*)Nch_dist->Clone(Form("shaded_bin_%d", k));
            for (int bin = 1; bin <= shaded->GetNbinsX(); ++bin) {
                if (bin < low_cut || bin > high_cut) {
                    shaded->SetBinContent(bin, 0);
                }
            }
            
//            shaded->SetFillColorAlpha(k % 2 == 0 ? kBlue : kRed, 0.3);
            int colorIndex = customColors[k % nCustomColors];
            shaded->SetFillColorAlpha(colorIndex, 0.3);
//            const int colorIndex = palette[k % nColors];
//            shaded->SetFillColorAlpha(colorIndex, 0.3);
            shaded->SetLineColor(0);
            shaded->Draw("SAME HIST F");
            
            double low_pct = per.GetFineAmplitudeCutLow(k) * 100.0;
            double high_pct = per.GetFineAmplitudeCutHigh(k) * 100.0;
            if (low_pct > high_pct) {
                std::cout << "Degenerate bin " << low_pct << " - " << high_pct << ". Ignoring" << std::endl;
                continue;
            }

            int precision = 0;
            
            if (low_pct < 0.001) precision = 4;
            else if (low_pct < 0.01) precision = 3;
            else if (low_pct < 0.1) precision = 2;
            else if (low_pct < 1.0) precision = 1;

            TString format;
            if (low_pct < 0.0000001) {format = Form("%%.%df - %%.%df%%%%", 0, precision);}
            format = Form("%%.%df - %%.%df%%%%", precision, precision);
            TString entryLabel = Form(format, low_pct, high_pct);
            leg->AddEntry(shaded, entryLabel, "f");
        }
    leg->Draw();
//    std::string outFilePath = prepdir + "/pdf_shaded_" + refMult + ".pdf";
    std::string outFilePath = prepdir + "/pdf_shaded.pdf";
    c1->Write();
    c1->SaveAs(outFilePath.c_str());
    
    TCanvas* cZoom = new TCanvas("cZoom", "pdf shaded zoomed", 800, 800);
    cZoom->SetLogy();
    cZoom->SetTicks(1, 1);

    TH1D* zoomed = (TH1D*)Nch_dist->Clone("zoomed_dist");
    zoomed->SetMinimum(1e-7);
    zoomed->SetMaximum(1e-1);
    zoomed->GetYaxis()->SetTitleSize(0.035);
    zoomed->GetYaxis()->SetTitleOffset(1.25);

    // Determine x-range for <10% centrality (high multiplicity)
//    int bin10pct = per.GetFineAmplitudeBinCutHigh(nperc_bins-1); // 35 bins → bin 34 is 90–100%
//    double xMin = zoomed->GetBinLowEdge(bin10pct);
    double xMax = zoomed->GetBinLowEdge(zoomed->GetNbinsX()) + zoomed->GetBinWidth(zoomed->GetNbinsX());
    
//    zoomed->GetXaxis()->SetRangeUser(400, 700);
    zoomed->GetXaxis()->SetRangeUser(100, 600);
    zoomed->Draw("HIST");

    TLegend* legZoom = new TLegend(0.12, 0.7, 0.89, 0.88);
    legZoom->SetNColumns(4);
    legZoom->SetBorderSize(0);
    legZoom->SetFillStyle(0);
    legZoom->SetTextSize(0.015);

    for (int k = 0; k < nperc_bins; k++) {
        double low_pct = per.GetFineAmplitudeCutLow(k)*100.0;
        double high_pct = per.GetFineAmplitudeCutHigh(k)*100.0;

        if (high_pct > 11.0) continue;  // skip anything less central than 11%
        if (low_pct > high_pct) continue;

        int low_cut = per.GetFineAmplitudeBinCutLow(k);
        if (low_cut < 0) low_cut = 1;
        int high_cut = per.GetFineAmplitudeBinCutHigh(k);
        if (k == 0) high_cut = Nch_dist->GetNbinsX();

        TH1D* shadedZoom = (TH1D*)Nch_dist->Clone(Form("shaded_zoom_bin_%d", k));
        for (int bin = 1; bin <= shadedZoom->GetNbinsX(); ++bin) {
            if (bin < low_cut || bin > high_cut) {
                shadedZoom->SetBinContent(bin, 0);
            }
        }

        int colorIndex = customColors[k % nCustomColors];
        shadedZoom->SetFillColorAlpha(colorIndex, 0.3);
        shadedZoom->SetLineColor(0);
        shadedZoom->Draw("SAME HIST F");
        
        int precision = 0;

        if (low_pct < 0.001) precision = 4;
        else if (low_pct < 0.05) precision = 4;
        else if (low_pct < 0.01) precision = 3;
        else if (low_pct < 0.06) precision = 3;
        else if (low_pct < 0.1) precision = 2;
        else if (low_pct < 1.0) precision = 1;
        
        precision = 1; //override
        
        TString format;
        if (low_pct < 0.000000001) {format = Form("%%.%df - %%.%df%%%%", 0, precision);}
        else {format = Form("%%.%df - %%.%df%%%%", precision, precision);}
        TString entryLabel = Form(format, low_pct, high_pct);
        legZoom->AddEntry(shadedZoom, entryLabel, "f");
    }

    legZoom->Draw();

    std::string outZoomPath = prepdir + "/pdf_shaded_zoomed_in.pdf";
    cZoom->Write();
    cZoom->SaveAs(outZoomPath.c_str());
    

// Construct CSV file path
std::string csvPath = prepdir + "/centrality_breakdown.csv";
std::ofstream csvOut(csvPath);
if (!csvOut.is_open()) {
    std::cerr << "❌ Could not open CSV file for writing: " << csvPath << std::endl;
    return;
}

// Write header
csvOut << "centrality regime (%), event selection,num events,error\n";

// Loop over all bins
for (int k = 0; k < nperc_bins; ++k) {
    int low_cut = per.GetFineAmplitudeBinCutLow(k);
    if (low_cut < 0) low_cut = 1;
    int high_cut = per.GetFineAmplitudeBinCutHigh(k);
    if (k == 0) high_cut = Nch_dist->GetNbinsX();

    double low_pct = per.GetFineAmplitudeCutLow(k) * 100.0;
    double high_pct = per.GetFineAmplitudeCutHigh(k) * 100.0;

    double content = 0;
    double error2 = 0;
    TString nchRange;

    for (int bin = low_cut; bin <= high_cut; ++bin) {
        content += Nch_dist->GetBinContent(bin);
        error2 += std::pow(Nch_dist->GetBinError(bin), 2);
    }

    double err = std::sqrt(error2);

    // Nch range text
    double Nch_min = Nch_dist->GetBinLowEdge(low_cut);
    double Nch_max = Nch_dist->GetBinLowEdge(high_cut + 1);
    nchRange = Form("%.1f–%.1f", Nch_min, Nch_max);

    // Format precision
    int precision = 0;
    if (low_pct < 0.001) precision = 4;
    else if (low_pct < 0.01) precision = 3;
    else if (low_pct < 0.1) precision = 2;
    else if (low_pct < 1.0) precision = 1;

    // Write to file
    csvOut << std::fixed << std::setprecision(precision)
           << low_pct << "-" << high_pct << "%,"
           << nchRange << ","
           << content << ","
           << err << "\n";
}

csvOut.close();
std::cout << "✅ CSV saved to: " << csvPath << std::endl;

}

