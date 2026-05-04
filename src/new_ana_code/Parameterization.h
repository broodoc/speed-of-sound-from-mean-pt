#pragma once

#include <string>
#include <vector>
#include <utility>
#include <memory>

class TH1;
class TH1D;
class TH2D;
class TGraphErrors;
class TList;

class OutputManager;
struct AnalysisConfig;

// This class becomes the "orchestrator" for the SoS pipeline.
// - No ROOT file ownership.
// - No giant physics utilities living here.
// - Minimal state: config + cross-step caches (e.g. fitResults if you still need them).
class Parameterization {
public:
  // Keep your enum here or move to a common header if multiple modules need it.
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

  explicit Parameterization(const AnalysisConfig& cfg, OutputManager& out);

  // Main pipeline entry point (you can keep your original name).
  void Ratio_compute_pt_Nch(TH2D* Pt_Nch_hist, TH1D* Nch_dist);

  // If you still want these ranges owned here, fine.
  // Alternatively: store in cfg and return cfg ranges from OutputManager/config.
//  const std::vector<std::pair<double,double>>& GetPtRanges()  const;
//  const std::vector<std::pair<double,double>>& GetFitRanges() const;
//  const std::vector<std::pair<double,double>>& GetSoSRanges() const;

private:
  // ---- Pipeline step helpers (orchestrator-only; implementations in .cpp) ----
  void ValidateInputs_(TH2D* Pt_Nch_hist, TH1D* Nch_dist) const;

  // This is where your nested loops live; each step delegates to module functions.
  void RunFullScan_(TH2D* Pt_Nch_hist, TH1D* Nch_dist);

  // If you track "best fit" across loops, keep a small struct here.
  struct BestFitTracker {
    bool   hasBest = false;
    double bestChi2Ndf = 1e300;
    double bestCs2 = 0.0;
    double bestCs2Err = 0.0;

    // Optional bookkeeping
    std::pair<double,double> bestPtRange{0,0};
    std::pair<double,double> bestFitRange{0,0};
    Function bestFunction = Function::LevyTsallis;
  };

  // Use a small “graphs bundle” type so you don’t pass 12 pointers around.
  struct GraphBundle {
    // Fill these out with your actual graphs
    TGraphErrors* gMeanPt = nullptr;
    TGraphErrors* gNch    = nullptr;
    TGraphErrors* gSoS    = nullptr;

    // Optional: lists/containers used for writing
    TList* list = nullptr;
  };

  // Create graphs for one (fit_range, pt_range, function) combo.
  GraphBundle MakeGraphs_(const std::pair<double,double>& fitRange,
                          const std::pair<double,double>& ptRange,
                          Function f);

  // Finalization per combo
  void WriteAndCleanup_(GraphBundle& gb);

private:
  const AnalysisConfig& cfg_;
  OutputManager& out_;

  // If you still store fitResults for error propagation, store them here.
  // Ideally FitUtils returns what ErrorPropagation needs (TFitResultPtr/cov matrix)
  // and you keep them local; but if you insist on member storage, keep it contained.
  std::vector<void*> fitResultsOpaque_; // placeholder: replace with TFitResultPtr if you include ROOT headers
  int fitResultIndex_ = -1;

  // Ranges: either owned here or pulled from cfg
  std::vector<std::pair<double,double>> pt_ranges_;
  std::vector<std::pair<double,double>> fit_ranges_;
  std::vector<std::pair<double,double>> sos_ranges_;
};
