#include "Run_Analysis.h"
#include "Parameterization.h"

void Run_Analysis(string dir,string prepdir,string refMult, string pTMult, string InEx, string DebugVerbose, string DisplayFitRanges, string MomentAna){
    
    string system = "AuAu_200";
    Parameterization per{};
    per.SetSystem(system);
    per.SetMultEstimator("TPC");
    per.SetCentralitySetting(InEx); //"Inclusive" or "Exclusive" centrality setting
    bool debug = false; bool momentAna = false;
    if(DebugVerbose == "Y"){per.SetDebugSetting(true); debug = true;}
    else{per.SetDebugSetting(false);}
    if(DisplayFitRanges == "Y"){per.SetDisplayFitRangesSetting(true);}
    else{per.SetDisplayFitRangesSetting(false);}
    if(MomentAna == "Y") momentAna = true;
    
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
        
        std::string ErrorFilename;
    
    for (const auto& entry : fs::directory_iterator(prepdir)) {
        if (entry.is_regular_file()) {
            std::string prepname = entry.path().filename().string();
            if (prepname.find("Error.txt") != std::string::npos) {
                ErrorFilename = entry.path().string();
                break; // stop at first match
            }
        }
    }
        
        std::ifstream infileError(ErrorFilename);
        std::vector<double> ErrorCuts;
        
        if (std::getline(infileError, line)) {
            std::stringstream ss(line);
            std::string item;
            
            while (std::getline(ss, item, ',')) {
                item.erase(0, item.find_first_not_of(" \t\n\r"));
                item.erase(item.find_last_not_of(" \t\n\r") + 1);
                
                if (!item.empty()) {
                    ErrorCuts.push_back(std::stod(item));
                }
            }
        }
        infileError.close();
        
        // Print the numbers to verify
        
        i = 0;
        cout << "Relative errors: ";
        for (const auto& num : ErrorCuts) {
            std::cout << num << ", ";
            i++;
        }
        cout << endl;
        
        std::cout << std::endl;
        per.SetFineBinCutsTPC(BinCuts.data());
        per.SetRelErrorCutsTPC(ErrorCuts.data()); //Set Relative Error here
        per.PrintFineBinCuts();
        bool isCoarseBinning = true;
        
    string dirfilename = dir + "temp";
    TFile *fIn = TFile::Open(dirfilename.c_str(), "READ");
    
    if (fIn->IsZombie()) {
        std::cout << "Error opening file" << std::endl; std::exit(EXIT_FAILURE);
    }
        
    TH1D* Nch_dist = (TH1D*)fIn->Get(refMult.c_str());
    if (!Nch_dist){
        cout << "\n\tError pulling Nch histogram from " << fIn->GetName() << '\n';
        cout << "\nModified Filename: " << refMult.c_str() << '\n';
        std::exit(EXIT_FAILURE);
    }
        
    TH2D* hPt_Nch =  (TH2D*)fIn->Get(pTMult.c_str());
    if (!hPt_Nch){
        cout << "\nError pulling hPt_Nch histogram from " << fIn->GetName() << '\n';
        cout << "Modified Filename: " << pTMult.c_str() << '\n';
        std::exit(EXIT_FAILURE);
    }
        cout << '\n' << "Running analysis on " << dirfilename << "... " << '\n';
        per.Ratio_compute_pt_Nch(hPt_Nch, Nch_dist, per); ///Fine
        std::cout << '\n' << "-=-==-=-==-=-==-=-==-=-==-=-==-=-==-=-==-=-==-=-==-=-==-=-==-"<< std::endl;
    if (momentAna) {
        TH2D* hPt_Mult_hist =  (TH2D*)fIn->Get("<pT>_Nch");
        TH2D* hNch_Mult_hist =  (TH2D*)fIn->Get("Nch_refMultCorr");
        cout << '\n' << "Running moment analysis on " << dirfilename << "... " << '\n';
        per.MomentAna(hPt_Nch, hPt_Mult_hist, hNch_Mult_hist, Nch_dist, per); 
    }
        gErrorIgnoreLevel = kFatal; //Ignoring useless and nonFatal writing-based errors
        std::cout << '\n' << "END OF Run_Analysis()"<< std::endl;
}
