#ifndef PARAMETERIZATION_H
#define PARAMETERIZATION_H

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

using std::cout;
using std::endl;
using std::string;

#include "TH1.h"
#include "TH2.h"
#include "TGraphErrors.h"
#include "TRandom3.h"
#include <TLegend.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include "Minuit2/Minuit2Minimizer.h"

//-----------Analysis Headers--------------
#include <chrono>
#include <thread>
#include <unistd.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <TROOT.h>
#include "TFile.h"
#include "TStyle.h"
#include "TString.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include <TH1.h>
#include <TH2.h>
#include "TF1.h"
#include "TLatex.h"
#include <TCanvas.h>
#include <TStyle.h>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <random>
#include "TRandom.h"
#include <TError.h>
#include <TGraph.h>
#include <TMath.h>
#include <TTree.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include "TPaveText.h"
//-----------Analysis Headers--------------

class Parameterization {
 public:
    std::vector<TFitResultPtr> fitResults;
    int index{0};
  enum class Function {
    LevyTsallis,
    Tsallis,
    Hagedorn,
    TBW,
    FokkerPlanck,
    PowerLaw,
    ExpDecay,
    RightKneeFunction,
    RightKneeFunctionApproximation
  };

  static std::string FunctionToString(Function f);
  static Function StringToFunction(const std::string& name);
    
  enum class PercFirstBin { kFirst, kLast };
  Parameterization();
  Parameterization(Function f);
  virtual ~Parameterization(){};

  double operator()(double *x, double *);
  void SetTrigger(string trigger) { fTrigger = trigger; }
  void SetBinCutFile(string bincutfile) { fbincutfile = bincutfile; }
  void SetMultEstimator(string estimator) { fMultEstimator = estimator; }
  void SetCentralitySetting(string setting) { fCentralitySetting = setting; } //Inclusive/Exclusive setting
  void SetDebugSetting(bool setdebug){ fDebugSetting = setdebug; }
  void SetDisplayFitRangesSetting(bool setDisplayFitRangesSetting){ fDisplayFitRangesSetting = setDisplayFitRangesSetting; }
  void SetSystem(string system) { fSystem = system; }
  void SetMC(bool isMC, bool isMCtrue) {
    fIsMC = isMC;
    fIsMCTrue = isMCtrue;
  }
  void SetFineBinCutsTPC(int* BinCuts) {
        for ( int i = 0; i < ( 49 ); ++i ) {
            FineBinCutsTPC[ i ] = BinCuts[ i ];
        }
    }
  void SetRelErrorCutsTPC(double* RelError) {
        for ( int i = 0; i < ( 49 ); ++i ) {
            RelErrorCutsTPC[ i ] = RelError[ i ];
        }
    }
  void FindPercentiles(TH1 *h, PercFirstBin nn, bool d, bool i);
  void Ratio_compute_pt_Nch(TH2D* Pt_Nch_hist, TH1D* Nch_dist, Parameterization params);
  void WriteBestChi2(Parameterization::Function BestChi2fitfunctype, std::pair<double, double> BestChi2fitrange, std::pair<double, double> pTrange, TFile* fout, TH2D* Pt_Nch_hist, TH1D* Nch_dist, Parameterization params);
  void NormalizeByBinWidth(TH1 *h);
  void UnnormalizeByBinWidth(TH1 *h);
  void PrintFineBinCuts();
    bool hasPointsInRange(TGraphErrors* graph, double xmin, double xmax);
    void saveToCSV(int nperc_bins, std::vector<int> high_bin_cuts, std::vector<int> low_bin_cuts, std::vector<double> low_perc_bins, std::vector<double> high_perc_bins, TH1* Nch_dist, Parameterization params);
  string GetTrigger() const { return fTrigger; }
  string GetMultEstimator() const { return fMultEstimator; }
  string GetCentralitySetting() const {return fCentralitySetting;}
  bool GetDebugSetting(){ return fDebugSetting; }
  bool GetDisplayFitRangesSetting(){ return fDisplayFitRangesSetting; }
  string GetBinCutFile() const { return fbincutfile; }
  string GetSystem() const { return fSystem; }
//  double GetdNdEta(const TH1 *h = nullptr, bool usecut = true,
//                   const double &cut = 10.0);
    TF1* GetFit(TH1* h, double low_cut = 0, double high_cut = 1., Parameterization::Function fitFuncType = Parameterization::Function::LevyTsallis, const TString& name = "");
//  TF1* GetFit(TH1 *h, double low_cut = 0, double high_cut = 1.0, Parameterization::Function fitFuncType = Parameterization::Function::LevyTsallis);
  int GetNParams(Parameterization::Function fitFuncType = Parameterization::Function::LevyTsallis);
  void SetDefaultFitParameters(TF1* fit, Parameterization::Function fitFuncType = Parameterization::Function::LevyTsallis);
//  int GetNParams(Function f);
//  void SetDefaultFitParameters(TF1* fit, Function f);
  double GetErrorMeanpT(const TH1 *h, TF1 *fit, bool usecut,
                                          const double &cut, bool debug = false);
  double GetErrordNdEta(const TH1 *h = nullptr, TF1 *fit = nullptr, bool usecut = true, const double &cut = 10.0, bool debug = false);
    //Is get error dNdEta the same as get error for GetMean? ( They're dNdEta are both evaluated over the same range )
//  double GetMean(const TH1 *h = nullptr, bool usecut = true, const double &cut = 10.0);
  void SetErrorYield(TH1D*& h, int centrality_index);
  double GetMean(TH1 *h = nullptr, TF1 *fit = nullptr, bool usecut = true, const double &high_cut = 10.0, const double &low_cut = 0.0);
  double GetdNdEta(TH1 *h = nullptr, TF1 *fit = nullptr, bool usecut = true, const double &high_cut = 10.0, const double &low_cut = 0.0);
  double GetNch(TH1 *h, const int low_ref, const int high_ref);
  double GetRange(const int k, Parameterization params, bool display);
  double GetDensity(const int k, Parameterization params);
  double GetCentrality(const int k, Parameterization params, bool display);
  static Double_t RightKneeFunction(Double_t *x, Double_t *par);
  static Double_t PowerLaw(Double_t *x, Double_t *par);
//    TGraphErrors* GetSoundFit(TGraphErrors* g, double low_cut, double high_cut);
    TGraphErrors* GetSoundFit(TGraphErrors* g, double low_cut, double high_cut, TList* outputList, double pTlow, double pThigh);
    bool IsMonotonicWithinError(double prev_x, double prev_y, double curr_x, double curr_y, double prev_x_err, double prev_y_err, double curr_x_err, double curr_y_err);
    double ComputeFitScore(double chi2ndf);
    TF1* GetFitcs2(TGraphErrors* g);
  double GetRelativeError(int bin) const;
  TH1* GetSpectraRatio(TH1 *h, TH1 *href, TF1 *fit, bool usecut,
                       const double &high_cut, const double &low_cut);
  TH1D* GetKneeParams(TH1D *h, int ref, bool norm, bool fit, bool ratio);
    TH1* GetDataFitRatio(TH1 *h, TF1* fit, std::string title);
    TH1* GetSpectraDifference(TH1 *h, TH1* ref, std::string title);
    
    void MomentAna(TH2D* Pt_Nch_hist, TH2D* Pt_Mult_hist, TH2D* Nch_Mult_hist, TH1D* Nch_dist, Parameterization params);
    double GetCumulant(TH1D* dist, int order);
    double GetCumulantError(TH1D* dist, int order);
        
  static double P_n(const double* x, const double* par);
  TH1D* GetRedCurve();
  TH1D* GetBlackCurve();
  TH1D* GetBlueCurve();
  TH1D* GetOrangeCurve();
  double f(double n);
  float GetCoarseAmplitudeCutLow(int bin) const;
  float GetCoarseAmplitudeCutHigh(int bin) const;
  float GetFineAmplitudeCutLow(int bin) const;
  float GetFineAmplitudeCutHigh(int bin) const;
  int GetCoarseAmplitudeBinCutLow(int bin) const;
  int GetCoarseAmplitudeBinCutHigh(int bin) const;
  int GetCoarseAmplitudeBinCutLowMC(int bin) const;
  int GetCoarseAmplitudeBinCutHighMC(int bin) const;
  int GetFineAmplitudeBinCutLow(int bin) const;
  int GetFineAmplitudeBinCutHigh(int bin) const;
  const char *GetEstimatorLatex() const;
  const char *GetCoarseAmplitudeBinLatex(int bin) const {
    return CoarsePercentileBinningLatex[bin];
  }
  const char *GetFineAmplitudeBinLatex(int bin) const;
  int GetnFinePercentileBinning() const;
  int GetnCoarsePercentileBinning() const;

  int GetColor() const;
  int GetMarker() const;
  string GetNameInputFile() const;
  static const std::vector<std::pair<double, double>>& GetPtRanges();
  static const std::vector<std::pair<double, double>>& GetfitRanges();
  static const std::vector<std::pair<double, double>>& GetsosRanges();
  static const std::vector<int>& Getbadpoints();

 private:
  Function type_f;
  string fTrigger;
  string fbincutfile;
  string fMultEstimator;
  string fCentralitySetting;
  bool fDebugSetting;
  bool fDisplayFitRangesSetting;
  string fSystem;
  bool fIsMC;
  bool fIsMCTrue;
  static constexpr double mass_pion{0.139570};
  static const int nCoarsePercentileBinning{9};
  static const int nFinePercentileBinning{30};
//  static const int nCentralFinePercentileBinning{29}; //28 //20 //7
  static const int nCentralFinePercentileBinning{29}; //28 //20 //7
    
    //-=-=-=-=-=-=-=-=-pT ranges (spectra selection)-=-=-=-=-=-=-=-=-
    inline static const std::vector<std::pair<double, double>> pT_ranges = {
        {0.0, 10.0},
//        {0.05, 10.0},
//        {0.1, 10.0},
//        {0.15, 10.0},
        {0.2, 10.0},
//        {0.0, 1.5},
//        {0.0, 1.7},
//        {0.0, 1.9},
//        {0.0, 2.0},
//        {0.0, 2.2},
//        {0.0, 2.4},
//        {0.0, 2.8},
//        {0.0, 3.0},
//        {0.0, 3.2},
//        {0.3, 10.0},
//        {0.4, 10.0},
//        {0.5, 10.0},
//        {0.6, 10.0},
//        {0.7, 10.0},
//        {0.8, 10.0},
//        {0.9, 10.0},
//        {1., 10.0},
//        {1.1, 10.0},
//        {1.2, 10.0},
//        {1.3, 10.0},
//        {1.4, 10.0},
//        {1.5, 10.0},
//        {2.0, 10.0},
//        {3.0, 10.0},
//        {4.0, 10.0},
//        {5.0, 10.0},
//        {6.0, 10.0},
//        {7.0, 10.0},
//        {8.0, 10.0},
//        {0.0, 0.2},
//        {0.0, 1.6},
//        {0.3, 10.0},
//        {0.4, 10.0},
        
//        
//        {0.0, 1.8},
//        {0.2, 1.8},
//        {0.4, 1.8},
//        
//        {0.0, 1.9},
//        {0.2, 1.9},
//        {0.4, 1.9},
        
//        {0.0, 0.4},
//        {0.0, 0.5},
//        {0.0, 0.6},
//        {0.0, 0.7},
//        {0.0, 0.8},
//        {0.0, 0.9},
//        {0.0, 1.},
        };
    //-=-=-=-=-=-=-=-=-pT ranges-=-=-=-=-=-=-=-=-
    
    //-=-=-=-=-=-=-=-=-fit ranges (for spectra)-=-=-=-=-=-=-=-=-
    inline static const std::vector<std::pair<double, double>> fit_ranges = {
        {0.0, 0.3},
        {0.0, 0.35},
        {0.0, 0.4},
        {0.0, 0.45},
        {0.0, 0.5},
        {0.0, 0.55},
        {0.0, 0.6},
        {0.0, 0.65},
        {0.0, 0.7},
        {0.0, 0.75},
        {0.0, 0.8},
        {0.0, 0.85},
        {0.0, 0.9},
        {0.0, 0.95},
        {0.0, 1.0},
//        {0.0, 1.05},
//        {0.0, 1.1},
//        {0.0, 1.15},
//        {0.0, 1.2},
//        {0.0, 1.25},
//        {0.0, 1.3},
//        {0.2, 1.},
//        {0.3, 1.},
//        {0.2, 1.2},
//        {0.0, 1.35},
//        {0.0, 1.4},
//        {0.0, 1.45},
//        {0.0, 1.5},
//        {0.0, 1.55},
//        {0.0, 1.6},
//        {0.0, 1.65},
//        {0.0, 1.7},
//        {0.0, 1.75},
//        {0.0, 1.8},
//        {0.0, 1.9},
//        {0.0, 2.0},
//        {0.0, 2.1},
//        {0.0, 2.2},
//        {0.0, 2.3},
//        {0.0, 2.4},
//        {0.0, 2.5},
//        {0.0, 3.0},
        };
    //-=-=-=-=-=-=-=-=-fit ranges (for spectra)-=-=-=-=-=-=-=-=-
    
    //-=-=-=-=-=-=-=-=-fit ranges (for SOS calculation)-=-=-=-=-=-=-=-=-
    inline static const std::vector<std::pair<double, double>> sos_ranges = {
//        {0.75, 1.25},
//        {0.8, 1.25},
//        {0.85, 1.25},
//        {0.90, 1.25},
//        {0.95, 1.25},
//        {0.975, 1.25},
//        {1.0, 1.3},
//        {1.02, 1.3},
//        {1.05, 1.3},
//        {1.075, 1.25},
//        {1.10, 1.25},
//        {1.125, 1.25},
//        {1.15, 1.25},
//        {1.175, 1.25},
        {1.02, 1.25},
        {1.03, 1.25},
        {1.04, 1.25},
        {1.05, 1.25},
        {1.06, 1.25},
        {1.065, 1.25},
        {1.07, 1.25},
        {1.08, 1.25},
        {1.09, 1.25},
//        {1.095, 1.25},
        {1.10, 1.25},
//        {1.10, 1.2},
//        {1.10, 1.25},
//        {1.10, 1.3},
//        {1.15, 1.3},
//        {1.105, 1.25},
        {1.11, 1.25},
        {1.12, 1.25},
        {1.13, 1.25},
//        {1.14, 1.3},
//        {1.15, 1.3},
//        {1.16, 1.3},
//        {1.17, 1.3},
//        {1.18, 1.3},
//        {1.19, 1.3},
//        {1.2, 1.3},
//        {1.21, 1.3},
//        {1.22, 1.3},
        
//        {1.11, 1.25},
//        {1.12, 1.25},
//        {1.13, 1.25},
//        {1.14, 1.25},
//        {1.15, 1.25},
//        {1.16, 1.25},
//        {1.17, 1.25},
//        {1.18, 1.25},
//        {1.19, 1.25},
//        {1.2, 1.25},
//        {1.21, 1.25},
//        {1.22, 1.25},
        };
    //-=-=-=-=-=-=-=-=-fit ranges (for SOS calculation)-=-=-=-=-=-=-=-=-
    
    //-=-=-=-=-=-=-=-=-data point omission (for bad points)-=-=-=-=-=-=-=-=-

//        inline static const std::vector<int> bad_points = { -1, 24, 25, 26, 27, 28, 29, 30, 31, };
        inline static const std::vector<int> bad_points = { -1, 24, 25, 26, 27, 28, 29, 30, 31, };
//        inline static const std::vector<int> bad_points = { -1, 0,1,2, 24, 25, 26, 27, 28, 29, 30, 31, };

    //-=-=-=-=-=-=-=-=-data point omission (for bad points)-=-=-=-=-=-=-=-=-
    
    //-=-=-=-=-=-=-=-=-Points to alter pT fit region-=-=-=-=-=-=-=-=-
//    inline static const std::vector<int> ater_ponts = { -1, 6, 0, };
    //-=-=-=-=-=-=-=-=-Points to alter pT fit region-=-=-=-=-=-=-=-=-
 
//    static const int nCentralFinePercentileBinning{37};
  //! The coarse binning is used to find the amplitude cut corresponding to
  //! the 5% events with the highest V0M amplitude.
  static constexpr float CoarsePercentileBinning[nCoarsePercentileBinning + 1] =
      {0.0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0};
  static constexpr float FinePercentileBinning[nFinePercentileBinning + 1] = {
      0.0,   0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007,
      0.008, 0.009, 0.01,  0.015, 0.02,  0.025, 0.03,  0.035,
      0.04,  0.045, 0.05,  0.055, 0.06,  0.065, 0.07,  0.075,
      0.08,  0.085, 0.09,  0.095, 0.1,   0.2,   1.0};
    /*
  // static constexpr float FinePercentileBinning[nFinePercentileBinning + 1] =
  // {
  //     0.0,   0.00025, 0.0005, 0.00075, 0.001, 0.002, 0.003, 0.004, 0.005,
  //     0.006, 0.007,   0.008,  0.009,   0.01,  0.015, 0.02,  0.025, 0.03,
  //     0.035, 0.04,    0.045,  0.05,    0.055, 0.06,  0.065, 0.07,  0.075,
  //     0.08,  0.085,   0.09,   0.095,   0.1,   0.2,   1.0};
     
     //-=-=-=-=-=-=-=-=-FINE PERCENTILE BINNING-=-=-=-=-=-=-=-=-
     
  */ // Other Fine Percentile Binning
  static constexpr float CentralFinePercentileBinning[nCentralFinePercentileBinning + 1] = {
      0.0, 0.001, 0.002, 0.003, 0.005, 0.006, 0.007, 0.008, 0.01,0.02, 0.03,0.04, 0.05, 0.055, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0  };
//    0.0, 0.0002, 0.0003, 0.0005, 0.001, 0.002, 0.003, 0.005, 0.01,0.02, 0.03,0.04, 0.05, 0.055, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0  };
//    0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03,0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0  };
    //0.0, 0.05, 0.1, 0.15 , 0.2, 0.25, 0.3, 1.0 // 0-5 is 0
    //0.0, 0.001, 0.002, 0.003, 0.005, 0.01,0.02, 0.03,0.04, 0.05, 0.055, 0.06, 0.07, 0.08, 0.09, 0.1, 0.3, 0.5, 0.7, 0.8,   1.0  // => nCentralFinePercentileBinning = 20, 0-5 is 9
    //0.0, 0.0001, 0.0002, 0.0003, 0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0055, 0.006, 0.007, 0.008, 0.009, 0.01, 0.03, 0.05, 0.07, 0.08, 0.09, 0.1, 0.3, 0.5, 0.7, 0.8,   1.0 // => nCentralFinePercentileBinning = 26, 0-5 is 17
    //0.0, 0.005, 0.01,0.03,0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.3, 0.5, 0.7, 0.8, 1.0}; // => nCentralFinePercentileBinning = 14, 0-5 is 4
    // Note, change accordingly: int ref_bin{ here };
    
//    {
//            0.0,   0.00005, 0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.00075,
//            0.001, 0.002,   0.003,  0.004,  0.005,  0.006,  0.007,  0.008,
//            0.009, 0.01,    0.015,  0.02,   0.025,  0.03,   0.035,  0.04,
//            0.045, 0.05,    0.055,  0.06,   0.065,  0.07,   0.075,  0.08,
//            0.085, 0.09,    0.095,  0.1,    0.2,    1.0}; // 38 Elements

  // ------------------------ Au+Au 200 GeV ------------------------

  static constexpr int CoarseBinCuts[nCoarsePercentileBinning + 1] = {
      650, 362, 325, 292, 232, 178, 132, 93, 38, 1}; //Edited - 10 elements
  static constexpr int CoarseBinCutsTPC[nCoarsePercentileBinning + 1] = {
      650, 362, 325, 292, 232, 178, 132, 93, 38, 1}; //Edited - 10 elements
  static constexpr int CoarseBinCutsEtaGapTPC[nCoarsePercentileBinning + 1] = {
      1250, 363, 311, 264, 193, 140, 100, 69, 28, 1};
  static constexpr int CoarseBinCutsTracklets10[nCoarsePercentileBinning + 1] =
      {1200, 632, 534, 449, 324, 233, 165, 113, 47, 1};
  static constexpr int CoarseBinCutsTracklets14[nCoarsePercentileBinning + 1] =
      {1200, 850, 705, 591, 426, 306, 216, 148, 61, 1};
  static constexpr int CoarseBinCutsEtaGapSPD[nCoarsePercentileBinning + 1] = {
      1200, 432, 338, 275, 195, 139, 97, 66, 27, 1};
    

  // MC cuts for efficiencies: Pb-Pb 5.02 TeV ------------------------

  static constexpr int CoarseBinCutsTPCMCTrue[nCoarsePercentileBinning + 1] = {
      2500, 769, 658, 557, 408, 298, 214, 149, 63, 1};
  static constexpr int CoarseBinCutsTPCMCMeas[nCoarsePercentileBinning + 1] = {
      2500, 492, 422, 358, 264, 194, 139, 97, 41, 1};
  static constexpr int
      CoarseBinCutsTrackletsEtaGapMCTrue[nCoarsePercentileBinning + 1] = {
          1200, 633, 544, 461, 338, 247, 177, 123, 52, 1};
  static constexpr int
      CoarseBinCutsTrackletsEtaGapMCMeas[nCoarsePercentileBinning + 1] = {
          1200, 422, 329, 268, 191, 137, 97, 67, 28, 1};
  static constexpr int
      CoarseBinCutsTracksEtaGapTPCMCTrue[nCoarsePercentileBinning + 1] = {
          2500, 585, 501, 424, 311, 227, 163, 113, 47, 1};
  static constexpr int
      CoarseBinCutsTracksEtaGapTPCMCMeas[nCoarsePercentileBinning + 1] = {
          2500, 364, 311, 263, 193, 141, 101, 70, 29, 1};
 
  // MC cuts for efficiencies: Pb-Pb 5.02 TeV ------------------------
 // MC Cuts for efficiencies
    
    /// ------------------------------------Au-Au 3 GeV-------------------------------------------------
    int FineBinCutsTPC[nCentralFinePercentileBinning + 1];
    double RelErrorCutsTPC[nCentralFinePercentileBinning + 1];
  static constexpr int FineBinCuts[nFinePercentileBinning + 1] = {
      650, 391, 383, 378, 374, 371, 368, 366, 364, 362, 360, 353,
      347, 342, 337, 333, 329, 325, 321, 317, 313, 309, 305, 301,
      297, 293, 289, 285, 281, 222, 1}; // Edited - 31 Elements (30 bins)
//    static constexpr int FineBinCutsTPC[nCentralFinePercentileBinning + 1] = {
//        1501, 211, 202, 197, 193, 189, 178, 165, 156, 149, 143, 133, 128, 121, 118, 116, 114, 112, 111, 110, 109, 108, 105, 102, 99, 97, 95, 93, 91, 89, 87, 85, 83, 81, 79, 77, 75, 73, 72, 71, 62, 54, 41, 30, 21, 14, 9, 5, 1}; /// Edited - 49 Elements in this set <---- Central Fine percentile binning (TPC)
    /// Values computed from Find Percentiles: 200_GeV_Run11_Analysis_SummaryEta0p5__pTlow0p0__May24.root (Run11)
    /// --------------------------------------------------------------------------------------------------------------------
    
//    {
//        1000, 516, 510, 503, 499, 496, 493, 489, 486, 478, 472, 468, 464, 461, 458, 455, 452, 450, 442, 435, 429, 423, 417, 412, 407, 402, 397, 392, 387, 382, 377, 372, 367, 362, 357, 352, 312, 275, 208, 151, 105, 69, 42, 23, 1}; // Edited - 45 Elements <---- Central Fine percentile binning (TPC)
    // Values computed from Find Percentiles: 200_GeV_Run11_Analysis_Summary_Feb28.root (Run11)
    
//    {1000, 516, 510, 503, 499, 496, 493, 489, 486, 478, 472, 468, 464, 461, 458, 455, 452, 450, 442, 435, 429, 423, 417, 412, 407, 402, 397, 392, 387, 382, 377, 372, 367, 362, 357, 352, 275, 1}; // Edited - 38 Elements <---- Central Fine percentile binning (TPC)
    // Values computed from Find Percentiles: 200_GeV_Run11_Analysis_Summary_Feb28.root (Run11)
//    {
//      650, 416, 411, 405, 401, 398, 396, 392, 389, 382, 377, 373,
//      370, 367, 365, 363, 361, 359, 352, 346, 341, 336, 332, 328,
//      324, 320, 316, 312, 308, 304, 300, 296, 292, 288, 284, 280,
//      221, 1}; // Edited - 38 Elements <---- Central Fine percentile binning (TPC)
    // Values computed from Find Percentiles: 200_GeV_Analysis_Summary1.root (Run10)
    static constexpr int FineBinCutsEtaGapTPC[nCentralFinePercentileBinning + 1] = {};
//      {1250, 416, 411, 406, 403, 400, 398, 394, 391, 384, 379, 375, 372,
//       369,  366, 364, 362, 360, 352, 344, 337, 330, 324, 318, 312, 306,
//       300,  295, 290, 285, 280, 275, 270, 265, 260, 255, 187, 1};
  static constexpr int
    FineBinCutsTracklets10[nCentralFinePercentileBinning + 1] = {};
//          1200, 726, 723, 715, 709, 705, 701, 695, 690, 677, 668, 660, 654,
//          648,  643, 638, 633, 629, 612, 598, 585, 573, 561, 550, 539, 529,
//          519,  509, 500, 491, 482, 473, 465, 457, 449, 441, 318, 1};
  static constexpr int
FineBinCutsTracklets14[nCentralFinePercentileBinning + 1] = {};
//          1200, 989, 979, 968, 961, 955, 950, 941, 934, 916, 903, 892, 883,
//          875,  867, 860, 853, 847, 821, 799, 779, 761, 744, 728, 713, 699,
//          686,  673, 661, 649, 637, 626, 615, 604, 593, 583, 420, 1};
  static constexpr int
FineBinCutsTrackletsEtaGap[nCentralFinePercentileBinning + 1] = {};
//          1200, 520, 514, 507, 502, 499, 496, 491, 487, 475, 467, 460, 454,
//          448,  443, 438, 433, 429, 411, 396, 383, 372, 362, 352, 343, 335,
//          327,  319, 312, 305, 298, 292, 286, 280, 274, 268, 190, 1};
    
  // ------------------------ Xe-Xe 5.44 TeV ------------------------

  static constexpr int CoarseBinCutsXeXe544[nCoarsePercentileBinning + 1] = {
      1720, 883, 761, 649, 482, 355, 257, 181, 79, 1};
  static constexpr int CoarseBinCutsTPCXeXe544[nCoarsePercentileBinning + 1] = {
      1250, 314, 265, 225, 165, 120, 86, 60, 26, 1};
  static constexpr int
      CoarseBinCutsTrackletsEtaGapXeXe544[nCoarsePercentileBinning + 1] = {
          1200, 266, 207, 168, 119, 86, 61, 42, 18, 1};
  static constexpr int FineBinCutsXeXe544[nFinePercentileBinning + 1] = {
      1720, 952, 935, 925, 916, 908, 901, 895, 889, 884, 879,
      859,  842, 826, 811, 796, 782, 769, 756, 743, 731, 719,
      707,  695, 684, 673, 662, 651, 641, 476, 1};
  static constexpr int FineBinCutsTPCXeXe544[nFinePercentileBinning + 1] = {
      1250, 348, 339, 333, 329, 325, 322, 319, 316, 313, 311,
      302,  295, 288, 282, 276, 271, 266, 261, 256, 251, 246,
      241,  237, 233, 229, 225, 221, 217, 159, 1};
  static constexpr int
      FineBinCutsTrackletsEtaGapXeXe544[nFinePercentileBinning + 1] = {
          1200, 306, 296, 290, 285, 280, 276, 272, 269, 266, 263,
          252,  243, 235, 228, 221, 215, 209, 204, 199, 194, 189,
          185,  181, 177, 173, 169, 165, 162, 115, 1};

  // ------------------------ Pb-Pb 2.76 TeV ------------------------
  //
  static constexpr int CoarseBinCutsPbPb276[nCoarsePercentileBinning + 1] = {
      1720, 786, 675, 572, 421, 308, 221, 153, 64, 1};
  static constexpr int CoarseBinCutsTPCPbPb276[nCoarsePercentileBinning + 1] = {
      1250, 368, 305, 256, 185, 133, 94, 64, 26, 1};
  static constexpr int
      CoarseBinCutsTrackletsEtaGapPbPb276[nCoarsePercentileBinning + 1] = {
          1200, 292, 225, 181, 127, 90, 63, 43, 17, 1};
  static constexpr int FineBinCutsPbPb276[nFinePercentileBinning + 1] = {
      1720, 838, 826, 818, 811, 805, 800, 795, 791, 787, 783,
      766,  750, 735, 721, 707, 694, 681, 669, 657, 645, 634,
      623,  612, 602, 592, 582, 572, 563, 414, 1};
  static constexpr int FineBinCutsTPCPbPb276[nFinePercentileBinning + 1] = {
      1250, 418, 406, 398, 392, 386, 381, 377, 373, 369, 366,
      354,  344, 336, 328, 321, 314, 308, 302, 296, 290, 284,
      279,  274, 269, 264, 259, 254, 249, 180, 1};
  static constexpr int
      FineBinCutsTrackletsEtaGapPbPb276[nFinePercentileBinning + 1] = {
          1200, 336, 326, 319, 313, 308, 304, 300, 296, 292, 289,
          276,  265, 256, 248, 241, 234, 228, 222, 216, 211, 206,
          201,  196, 192, 188, 184, 180, 176, 124, 1};
 //Xe+Xe & Pb+Pb 2.76 TeV
    
  static constexpr const char
      *CoarsePercentileBinningLatex[nCoarsePercentileBinning + 1] = {
          "0-1%",   "1-5%",   "5-10%",  "10.20%", "20-30%",
          "30-40%", "40.50%", "50-70%", "70-100%"};
  static constexpr const char *EstimatorLatex[6] = {
      "V0M",
      "#it{N}_{ch}^{0#leq#eta#leq0.8}",
      "#it{N}_{tracklets}^{|#eta|<1",
      "#it{N}_{tracklets}^{|#eta|<1.4}",
      "#it{N}_{tracklets}^{0.7<|#eta|<1.4}",
      "#it{N}_{ch}^{0.5<|#eta|<0.8}"};

  static constexpr const char
      *FinePercentileBinningLatex[nFinePercentileBinning] = {
          "0-0.1%",   "0.1-0.2%", "0.2-0.3%", "0.3-0.4%", "0.4-0.5%",
          "0.5-0.6%", "0.6-0.7%", "0.7-0.8%", "0.8-0.9%", "0.9-1%",
          "1-1.5%",   "1.5-2%",   "2-2.5%",   "2.5-3%",   "3-3.5%",
          "3.5-4%",   "4-4.5%",   "4.5-5%",   "5-5.5%",   "5.5-6%",
          "6-6.5",    "6.5-7%",   "7-7.5%",   "7.5-8%",   "8-8.5%",
          "8.5-9%",   "9-9.5%",   "9.5-10%",  "10-20%",   "20-100%"};

  // static constexpr const char
  //     *FinePercentileBinningLatex[nFinePercentileBinning] = {
  //         "0-0.05%",   "0.2-0.25%", "0.25-0.3%", "0.3-0.35%", "0.35-0.4%",
  //         "0.4-0.45%", "0.45-0.5%", "0.5-0.6%",  "0.6-0.7%",  "0.7-0.8%",
  //         "0.8-0.9%",  "0.9-1%",    "1-1.5%",    "1.5-2%",    "2-2.5%",
  //         "2.5-3%",    "3-3.5%",    "3.5-4%",    "4-4.5%",    "4.5-5%",
  //         "5-5.5%",    "5.5-6%",    "6-6.5",     "6.5-7%",    "7-7.5%",
  //         "7.5-8%",    "8-8.5%",    "8.5-9%",    "9-9.5%",    "9.5-10%",
  //         "10-20%",    "20-100%"};

  static constexpr const char
    *CentralFinePercentileBinningLatex[nCentralFinePercentileBinning] = {};
//          "0-0.005%",   "0.005-0.01%", "0.01-0.02%", "0.02-0.03%", "0.03-0.04%",
//          "0.04-0.05%", "0.05-0.075%", "0.075-0.1%", "0.1-0.2%",   "0.2-0.3%",
//          "0.3-0.4%",   "0.4-0.5%",    "0.5-0.6%",   "0.6-0.7%",   "0.7-0.8%",
//          "0.8-0.9%",   "0.9-1%",      "1-1.5%",     "1.5-2%",     "2-2.5%",
//          "2.5-3%",     "3-3.5%",      "3.5-4%",     "4-4.5%",     "4.5-5%",
//          "5-5.5%",     "5.5-6%",      "6-6.5%",     "6.5-7%",     "7-7.5%",
//          "7.5-8%",     "8-8.5%",      "8.5-9%",     "9-9.5%",     "9.5-10%",
//          "10-20%",     "20-100%"};

  static constexpr const char *NameInputFile[3] = {"AnalysisResults_PbPb276",
                                                   "AnalysisResults_PbPb502",
                                                   "AnalysisResults_XeXe544"};
};

#endif
