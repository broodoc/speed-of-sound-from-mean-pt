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

class Parameterization {
 public:
  enum class Function {
    LevyTsallis,
    Tsallis,
    Hagedorn,
    FokkerPlanck,
    PowerLaw,
    ExpDecay,
    RightKneeFunction
  };
  enum class PercFirstBin { kFirst, kLast };
  Parameterization();
  Parameterization(Function f);
  virtual ~Parameterization(){};

  double operator()(double *x, double *);
  void SetTrigger(string trigger) { fTrigger = trigger; }
  void SetBinCutFile(string bincutfile) { fbincutfile = bincutfile; }
  void SetMultEstimator(string estimator) { fMultEstimator = estimator; }
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
  void FindPercentiles(TH1 *h, PercFirstBin nn, bool d, bool i);
    void Ratio_compute_pt_Nch(TH2D* Pt_Nch_hist, TH1D* Nch_dist, bool isCoarseBinning, Parameterization params);
  void NormalizeByBinWidth(TH1 *h);
  void PrintFineBinCuts();
    void saveToCSV(int nperc_bins, std::vector<int> high_bin_cuts, std::vector<int> low_bin_cuts, std::vector<double> low_perc_bins, std::vector<double> high_perc_bins, TH1* Nch_dist, Parameterization params);

  string GetTrigger() const { return fTrigger; }
  string GetMultEstimator() const { return fMultEstimator; }
  string GetBinCutFile() const { return fbincutfile; }
  string GetSystem() const { return fSystem; }
//  double GetdNdEta(const TH1 *h = nullptr, bool usecut = true,
//                   const double &cut = 10.0);
  TF1* GetFit(TH1 *h);
  double GetErrordNdEta(const TH1 *h = nullptr, bool usecut = true, const double &cut = 10.0);
    //Is get error dNdEta the same as get error for GetMean? ( They're dNdEta are both evaluated over the same range )
//  double GetMean(const TH1 *h = nullptr, bool usecut = true, const double &cut = 10.0);
  
  double GetMean(TH1 *h = nullptr, TF1 *fit = nullptr, bool usecut = true, const double &high_cut = 10.0, const double &low_cut = 0.0);
  double GetdNdEta(TH1 *h = nullptr, TF1 *fit = nullptr, bool usecut = true, const double &high_cut = 10.0, const double &low_cut = 0.0);
    double GetNch(TH1 *h, const int low_ref, const int high_ref);
  static Double_t RightKneeFunction(Double_t *x, Double_t *par);
    TF1* GetFitcs2(TGraphErrors* g);
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

 private:
  Function type_f;
  string fTrigger;
  string fbincutfile;
  string fMultEstimator;
  string fSystem;
  bool fIsMC;
  bool fIsMCTrue;
  static constexpr double mass_pion{0.139570};
  static const int nCoarsePercentileBinning{9};
  static const int nFinePercentileBinning{30};
  static const int nCentralFinePercentileBinning{23};

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
  */ // Other Fine Percentile Binning
  static constexpr float CentralFinePercentileBinning[nCentralFinePercentileBinning + 1] = {
          0.0, 0.0001, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01,0.02, 0.03,0.04, 0.05, 0.055, 0.06, 0.07, 0.08,
          0.09, 0.1, 0.3, 0.5, 0.7, 0.8,   1.0}; //
    
//    0.0, 0.005, 0.01,0.03,0.05, 0.06, 0.07, 0.08,
//    0.09, 0.1, 0.3, 0.5, 0.7, 0.8,   1.0}; // => nCentral = 14, 0-5 is 4
    
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
