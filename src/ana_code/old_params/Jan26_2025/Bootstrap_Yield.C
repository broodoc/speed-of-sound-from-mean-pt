#include "Parameterization.h"
#include "Bootstrap_Yield.h"

void Bootstrap_Yield(string dir, string prepdir, string refMult, string pTMult, string InEx, string DebugVerbose){
    
    bool debug = false;
    if (DebugVerbose == "Y"){debug = true;}
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

    if (debug) std::cout << "✅ Found bin cuts file: " << filename << std::endl;
    
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
        if (debug) cout << "Bin cuts: ";
        for (const auto& num : BinCuts) {
            if (debug) std::cout << num << ", ";
            i++;
        }
        cout << endl;
        
        
        if (debug) std::cout << std::endl;
        per.SetFineBinCutsTPC(BinCuts.data());
        if (debug) per.PrintFineBinCuts();
    
    
    std::vector<std::string> filenames;
    int B = 10;
    for (int b=1; b<=B; ++b) filenames.push_back(Form("bootstrap/boot_%d.root", b));

//    std::vector<std::string> filenames = {
//            "merged_1.root",
//            "merged_2.root",
//            "merged_3.root",
//            "merged_4.root",
//            "merged_5.root",
//            "merged_6.root",
//            "merged_7.root",
//            "merged_8.root",
//            "merged_9.root",
//            "merged_10.root"
//        };
    
    TList* lOut = new TList();
    lOut->SetOwner(true);
    int nperc_bins{per.GetnFinePercentileBinning()};
    std::map<int, std::map<double, std::vector<double>>> centralityBinContentMap; // Centrality -> (pT bin center -> Yields)
    std::map<int, std::vector<double>> centralityEventCount; // k -> [Nevt for each subset]
    int m_index{0}; //subset index (1/10th of total dataset)
    // Loop for all #'s
    
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
    
    for (const auto& mfilename : filenames) {
        if (debug) cout << "Processing subset file: " << mfilename.c_str() << "..." << endl;
        m_index++;
        std::string fullFilename = dir + mfilename;
        std::ifstream infile(fullFilename);
        if (debug) cout << "\nfullFilename: " << fullFilename << endl;
        TFile *file = TFile::Open(fullFilename.c_str(), "READ");
        
        TH2D* hPt_Nch =  (TH2D*)file->Get(pTMult.c_str());
        TH1D* Nch_dist = (TH1D*)file->Get(refMult.c_str());
        
        if (!hPt_Nch || !Nch_dist) {
            std::cerr << "Error: Missing histograms in file " << mfilename << std::endl;
            if (!hPt_Nch) std::cerr << "  → Missing: " << pTMult << std::endl;
            if (!Nch_dist) std::cerr << "  → Missing: " << refMult << std::endl;
            continue;
        }
        
        hPt_Nch->SetName(Form("hPt_Nch_%s", mfilename.c_str()));
        Nch_dist->SetName(Form("Nch_dist_%s", mfilename.c_str()));
        
        lOut->Add(hPt_Nch);
        lOut->Add(Nch_dist);
        
        if (debug) std::cout << "Adding histogram: ";
        for (int k = 0; k < nperc_bins; k++) {
            
            TH1D* hpt = nullptr;
            
            if (InEx == "I"){ //Inclusive analysis
                hpt = hPt_Nch->ProjectionY(Form("#%d: %.5f-%.5f%%",
                                                      m_index,
                                                      per.GetFineAmplitudeCutLow(k)*100,
                                                      per.GetFineAmplitudeCutHigh(0)*100),
                                                 per.GetFineAmplitudeBinCutLow(k),
                                                 per.GetFineAmplitudeBinCutHigh(0));
            }
            else{             //Exclusive analysis
                hpt = hPt_Nch->ProjectionY(Form("#%d: %.5f-%.5f%%",
                                                      m_index,
                                                      per.GetFineAmplitudeCutLow(k)*100,
                                                      per.GetFineAmplitudeCutHigh(k)*100),
                                                 per.GetFineAmplitudeBinCutLow(k),
                                                 per.GetFineAmplitudeBinCutHigh(k));
            }
            
            hpt->Sumw2();
//            hpt->SetYTitle("1/N_{evt} 1/(2#pi p_{T}) d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-2}"); /// Rutik/Rene's version (differential yield)
            hpt->SetYTitle("1/N_{evt} d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-1}"); /// Omar's version (differential yield)
            hpt->SetXTitle("p_{T} (GeV/c)");
            
            double Nevt = 0.0;
            if (InEx == "I") {
              Nevt = Nch_dist->Integral(per.GetFineAmplitudeBinCutLow(k),
                                        per.GetFineAmplitudeBinCutHigh(0));
            } else {
              Nevt = Nch_dist->Integral(per.GetFineAmplitudeBinCutLow(k),
                                        per.GetFineAmplitudeBinCutHigh(k));
            }
            centralityEventCount[k].push_back(Nevt); // k -> [Nevt for each subset]

            
            if (InEx == "I"){ //Inclusive analysis
            hpt->Scale(1./Nch_dist->Integral(per.GetFineAmplitudeBinCutLow(k),per.GetFineAmplitudeBinCutHigh(0))); //Inclusive
            }
            else{             //Exclusive analysis
            hpt->Scale(1./Nch_dist->Integral(per.GetFineAmplitudeBinCutLow(k),per.GetFineAmplitudeBinCutHigh(k)));
            }
            per.NormalizeByBinWidth(hpt);
            lOut->Add(hpt);
            if (debug) std::cout << hpt->GetName() << ",";
            
            //Storing values
            for (int i = 1; i <= hpt->GetNbinsX(); ++i) {
                double binCenter = hpt->GetBinCenter(i);
                double binContent = hpt->GetBinContent(i);
                centralityBinContentMap[k][binCenter].push_back(binContent);
            }
            
        }
        int k{nperc_bins};
        //Reference bin (0-5%)
        int ref_bin{12}; //9 -> ref_bin = 8, 17 -> ref_bin = 16
//        int ref_bin{17}; //9 -> ref_bin = 8, 17 -> ref_bin = 16
        const int low_ref{per.GetFineAmplitudeBinCutLow(ref_bin)};
        //        const int low_ref{( per.GetFineAmplitudeBinCutHigh(ref_bin) + per.GetFineAmplitudeBinCutLow(ref_bin) ) / 2}; /// Obtaining reference bins ---> There are 49 elements (30 corresponds to 5% bin):: Here is where the data needs to be pulled from .txt
        const int high_ref{per.GetFineAmplitudeBinCutHigh(0)}; /// Obtaining reference bins:: Here is where the data needs to be pulled from .txt
        
        TH1D* hpt_ref = hPt_Nch->ProjectionY(Form("#%d: 0-5%%", m_index), low_ref,high_ref);
        hpt_ref->Sumw2();
        hpt_ref->SetTitle("0-5%");
//        hpt_ref->SetYTitle("1/N_{evt} 1/(2#pi p_{T}) d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-2}"); /// Rutik/Rene's version (differential yield)
        hpt_ref->SetYTitle("1/N_{evt} d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-1}"); /// Omar's version (differential yield)
        hpt_ref->SetXTitle("p_{T} (GeV/c)");
        
        double Nevt = 0.0;
          Nevt = Nch_dist->Integral(low_ref, high_ref);
        centralityEventCount[k].push_back(Nevt); // k -> [Nevt for each subset]

        hpt_ref->Scale(1./Nch_dist->Integral(low_ref,high_ref));
        per.NormalizeByBinWidth(hpt_ref);
        lOut->Add(hpt_ref);
        if (debug) std::cout << hpt_ref->GetName() << ",";
        
        //Storing values: REF is addressed to the last centrality (k) bin
        for (int i = 1; i <= hpt_ref->GetNbinsX(); ++i) {
            double binCenter = hpt_ref->GetBinCenter(i);
            double binContent = hpt_ref->GetBinContent(i);
            centralityBinContentMap[k][binCenter].push_back(binContent);
        }
        if (debug) std::cout << std::endl;
    }
    
    std::string yieldErrDir = prepdir + "/yield_errors";
    std::filesystem::create_directories(yieldErrDir);
    
    cout << "Computing errors..." << endl;
    for (int k = 0; k < nperc_bins + 1; k++) {
        
            std::string txtFilePath = yieldErrDir + "/" + Form("yield_error_%d.txt", k);
            std::ofstream outFile(txtFilePath);
        
//            std::ofstream outFile(Form("yield_errors/yield_error_%d.txt", k));
            if (!outFile.is_open()) {
                std::cerr << "Error: Unable to open file for writing for centrality " << k << std::endl;
                continue;
            }
        if (k == nperc_bins){outFile << Form("#%d: 0-5%%\n", k);}
        else {
            if (InEx == "I"){ //Inclusive analysis
                outFile << Form("#%d: %.3f-%.3f%%\n", k, per.GetFineAmplitudeCutLow(0) * 100, per.GetFineAmplitudeCutHigh(k) * 100);
            }
            else{             //Exclusive analysis
                outFile << Form("#%d: %.3f-%.3f%%\n", k, per.GetFineAmplitudeCutLow(k) * 100, per.GetFineAmplitudeCutHigh(k) * 100);
            }
            
        }
        
            
            outFile << "pT (GeV/c) | std dev | mean\n";

            TGraphErrors* graph = new TGraphErrors();
            graph->SetName(Form("bootstrapped_yield_%d", k));
            graph->GetXaxis()->SetTitle("p_{T} (GeV/c)");
//            graph->GetYaxis()->SetTitle("1/N_{evt} 1/(2#pi p_{T}) d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-2}");
            graph->GetYaxis()->SetTitle("1/N_{evt} d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-1}");
            
        const auto& w = centralityEventCount[k]; // size = Nsub
        const int Nsub = (int)w.size();
        
            int pointIndex = 0;
            for (const auto& [binCenter, binContents] : centralityBinContentMap[k]) {
                if ((int)binContents.size() != Nsub) continue;
                
                //event count weights
                double W = 0.0, Wy = 0.0, W2 = 0.0;
                for (int i=0;i<Nsub;i++){
                  const double wi = w[i];
                  W  += wi;
                  Wy += wi * binContents[i];
                  W2 += wi * wi;
                }
                if (W <= 0) continue;
                
                double mean = Wy/W;
                
//                for (double content : binContents) {
//                    mean += content;
//                }
//                mean /= binContents.size();

                double num = 0.0;
                for (int i=0;i<Nsub;i++){
                  const double wi = w[i];
                  const double di = binContents[i] - mean;
                  num += wi * di * di;
                }

                const double denom = W - (W2 / W);
                double s = 0.0;
                if (denom > 0) s = std::sqrt(num / denom);
                
//                const double Neff = (W2 > 0) ? (W*W / W2) : 0.0;
                const double Neff = binContents.size();

                const double err_spread = s;
                const double err_mean   = (Neff > 0) ? (s / std::sqrt(Neff)) : 0.0;
                
                const double stdDev = err_spread; // or err_mean
//                const double stdDev = err_mean;
                
//                double variance = 0.0;
//                for (double content : binContents) {
//                    variance += std::pow(content - mean, 2);
//                }
//                variance /= (binContents.size() - 1); // N-1 for sample std deviation
//                double stdDev = std::sqrt(variance);
//                stdDev /= std::sqrt(binContents.size());

                outFile << binCenter << " | " << stdDev << " | " << mean << "\n";

                graph->SetPoint(pointIndex, binCenter, mean);
                graph->SetPointError(pointIndex, 0.0, stdDev); // No error on x
                pointIndex++;
            }

            outFile.close();
            lOut->Add(graph);
        }
    
    std::string outFilePath = yieldErrDir + "/differential_yields.root";
    TFile* fOut = new TFile(outFilePath.c_str(), "recreate");
    
    if (debug){
        std::cout
        << "Number of objects in lOut: "
        << lOut->GetSize()
        << " (" << nperc_bins << " centrality classes X "
        << "10 subsets)"
        << std::endl;
    }

    if (debug){
        TIter nextDebug(lOut);
        TObject* debugObj;
        while ((debugObj = nextDebug())) {
            std::cout << "Object: " << debugObj->GetName()
            << ", Type: " << debugObj->ClassName()
            << ", IsZombie: " << debugObj->IsZombie() << std::endl;
        }
    }


    fOut->Write();
    fOut->cd();
    lOut->Write();
    delete fOut;
    delete lOut;
    
    gROOT->SetBatch(kFALSE);
    
    std::cout << "\nYield errors written to yield_errors/yield_errors.txt" << std::endl;
    std::cout << "Subset charged differential yields written to yield_errors/differential_yeilds.root" << std::endl;
}
