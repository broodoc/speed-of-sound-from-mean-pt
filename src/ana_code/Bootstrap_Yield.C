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
    
    
    TList* lOut = new TList();
    lOut->SetOwner(true);
    int nperc_bins{per.GetnFinePercentileBinning()};
    std::map<int, std::map<double, std::vector<double>>> centralityBinContentMap; // Centrality -> (pT bin center -> Yields)
    std::map<int, std::map<double, std::vector<double>>> centralityBinContentMapY;
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

            
            int cLow;
            int cHigh;
            
            for (int i = 1; i <= hpt->GetNbinsX(); ++i) {
                double binCenter = hpt->GetBinCenter(i);
                double binContent = hpt->GetBinContent(i);
                centralityBinContentMapY[k][binCenter].push_back(binContent);
            }
            if (InEx == "I"){ //Inclusive analysis
            hpt->Scale(1./Nch_dist->Integral(per.GetFineAmplitudeBinCutLow(k),per.GetFineAmplitudeBinCutHigh(0))); //Inclusive
                cLow = k;
                cHigh = 0;
            }
            else{             //Exclusive analysis
            hpt->Scale(1./Nch_dist->Integral(per.GetFineAmplitudeBinCutLow(k),per.GetFineAmplitudeBinCutHigh(k)));
                cLow = k;
                cHigh = k;
            }
            per.NormalizeByBinWidth(hpt);
            lOut->Add(hpt);
            if (debug) std::cout << hpt->GetName() << ",";


            int xLow  = per.GetFineAmplitudeBinCutLow(k);
            int xHigh = (InEx=="I") ? per.GetFineAmplitudeBinCutHigh(0)
                                    : per.GetFineAmplitudeBinCutHigh(k);
            
            double N0 = Nch_dist->Integral(xLow, xHigh); //total events for the centrality bin
            if (N0 <= 0) continue;
            
            TRandom3 rng(12345);
            // Loop over pT bins
            for (int j = 1; j <= hPt_Nch->GetNbinsY(); ++j) {

                double Y_CBWC = 0.0;
                double denom = 0.0, numer = 0.0;
                // Loop over microbins c inside K
                for (int c = xLow; c <= xHigh; ++c) {
                    
                    TH1D* hpt_c = hPt_Nch->ProjectionY(
                        Form("hpt_c_k%d_c%d_b%d", k, c, m_index),
                        c,
                        c
                    );
                    
                    double Nc = Nch_dist->GetBinContent(c);
                    const double w = (double)rng.Poisson(Nc);
                    

                    if (Nc <= 0) continue;

                    hpt_c->Scale(1.0 / Nc);
                    per.NormalizeByBinWidth(hpt_c);

                    double Yc = hpt_c->GetBinContent(j);
                    numer += w * Yc;
                    denom += w;
//                    Y_CBWC += (Nc / N0) * Yc;
//                    Y_CBWC += (Nc / N0) * Yc;

                    delete hpt_c;
                }
                
                Y_CBWC += ( numer / denom );

                double binCenter = hPt_Nch->GetYaxis()->GetBinCenter(j);
                centralityBinContentMap[k][binCenter].push_back(Y_CBWC);
            }


            //Storing values
//            for (int i = 1; i <= hpt->GetNbinsX(); ++i) {
//                double binCenter = hpt->GetBinCenter(i);
//                double binContent = hpt->GetBinContent(i);
//                centralityBinContentMapY[k][binCenter].push_back(binContent);
//            }
            
        }
        int k{nperc_bins};
        //Reference bin (0-5%)
        int ref_bin{27}; //9 -> ref_bin = 8, 17 -> ref_bin = 16
        const int low_ref{per.GetFineAmplitudeBinCutLow(ref_bin)};
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
            centralityBinContentMapY[k][binCenter].push_back(binContent);
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
        if (debug){
            outFile << "pT (GeV/c) | "
                    << "err unweighted mean | "
                    << "mean unweighted | "
                    << "sqrt(W2/W^2) | "
                    << "mean weighted | "
                    << "err weighted_mean | "
                    << "W = sum(w_i) | "
                    << "W2 = sum(w_i^2) | "
                    << "Effective Ntrks | "
                    << " Y | "
                    << " mean Nevts | "
                    << " RMS Nevts \n";
        }

        
        

            TGraphErrors* graph = new TGraphErrors();
            graph->SetName(Form("bootstrapped_yield_%d", k));
            graph->GetXaxis()->SetTitle("p_{T} (GeV/c)");
//            graph->GetYaxis()->SetTitle("1/N_{evt} 1/(2#pi p_{T}) d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-2}");
            graph->GetYaxis()->SetTitle("1/N_{evt} d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-1}");
            
        const auto& w = centralityEventCount[k];
        const int Nsub = (int)w.size();
        
            int pointIndex = 0;
            for (const auto& [binCenter, binContents] : centralityBinContentMap[k]) {
                if ((int)binContents.size() != Nsub) continue;
                
                auto itY = centralityBinContentMapY[k].find(binCenter);
                if (itY == centralityBinContentMapY[k].end()) continue;

                const std::vector<double>& binContentsY = itY->second;
            
                double NtrkEff = 0.0;
                double Ntrk_W = 0.0;
                double Ntrk_W2 = 0.0;
                for (int i=0;i<Nsub;i++){
                  const double wi = w[i];
                  Ntrk_W += binContentsY[i];
                  Ntrk_W2 += binContentsY[i] * binContentsY[i] ;
                }
                
                NtrkEff = (Ntrk_W*Ntrk_W) / Ntrk_W2;
                Ntrk_W = Ntrk_W / Nsub;
                //https://en.wikipedia.org/wiki/Weighted_arithmetic_mean
//                event count weights
                double W = 0.0, Wy = 0.0, W2 = 0.0; double mean = 0;
                for (int i=0;i<Nsub;i++){
                  const double wi = w[i];
                  W  += wi;
                  Wy += wi * binContents[i];
                  W2 += wi * wi;
                  mean += binContents[i];
                }
                if (W <= 0) continue;
                
                double Wm = W / Nsub;
                double sig_W = 0;
                for (int i=0;i<Nsub;i++){
                  const double wi = w[i];
                  const double di = wi - Wm;
                    sig_W += di*di;
                }
                
                sig_W = sig_W / (Nsub - 1);
                sig_W = std::sqrt(sig_W);
                mean = mean / Nsub;
                double mean_w = Wy/W; //event weighted mean

                double num = 0.0;
                for (int i=0;i<Nsub;i++){
//                  const double wi = w[i];
                  const double di = binContents[i] - mean;
//                  num += wi * di * di;
                  num += di * di;
                }

                double V_WEIGHT = W2 / (W * W);
                double S_WEIGHT = std::sqrt(V_WEIGHT);
//                const double denom = W - (W2 / W);
                const double denom = Nsub - 1;
                double s = 0.0;
                if (denom > 0) s = std::sqrt(num / denom);
                
                const double Neff = binContents.size();
                const double err_spread = s;
//                const double stdDev = err_spread * S_WEIGHT;
                const double stdDev = err_spread;
//                const double stdDev = 0.5;
                
                if (debug) std::cout << "S_WEIGHT: " << S_WEIGHT << std::endl;
                if (debug){
                    outFile << binCenter << " | " << stdDev << " | " << mean << " | " << S_WEIGHT << " | " << mean_w << " | " << s*S_WEIGHT
                            << " | " << W << " | " << W2 << " | " << NtrkEff << " | " << Ntrk_W << " | " << W / Nsub << " | " << sig_W << "\n";
                }
                else outFile << binCenter << " | " << stdDev << " | " << mean << "\n";
                
//                double mean = 0;
//                for (int i=0;i<Nsub;i++) mean += binContents[i];
//                mean /= Nsub;
//
//                double var = 0;
//                for (int i=0;i<Nsub;i++){
//                  double d = binContents[i] - mean;
//                  var += d*d;
//                }
//                var /= (Nsub - 1);
//                double stdDev = std::sqrt(var);        // spread across replicas
////                double err_mean = s/std::sqrt(Nsub); // optional

//                graph->SetPoint(pointIndex, binCenter, mean);
                
                
                graph->SetPoint(pointIndex, binCenter, mean_w);
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
