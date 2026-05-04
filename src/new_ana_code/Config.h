//Copyright 2026 Caleb Broodo
//Analysis configuration
#pragma once
#include <string>
#include <vector>
#include <set>
#include <utility>
#include <stdexcept>
#include <algorithm>
#include "Models.h"

enum class CentralitySetting { Inclusive, Exclusive };
enum class PercFirstBin { kFirst, kLast };
enum class Observable { MeanPT, dNdEta, Nevts, Nch, Cs2 };

struct Config {
  // Output / bookkeeping
  std::string outputFile = "output.root";
  std::string outputTag  = "default";

  // System / run info
  std::string system = "AuAu_200";
  std::string trigger = "NA";
  std::string multEstimator = "refMult";
  CentralitySetting centralitySetting = CentralitySetting::Inclusive;
  bool isMC = false;
  bool isMCTrue = false;

  std::set<int> goodRunIds;
  std::set<int> badRunIds;

  bool batchMode = true;
  bool debug = false;
  bool displayFitRanges = false;

  std::vector<int> BadPoints;
    
  PercFirstBin percentileDirection = PercFirstBin::kLast;
  std::vector<float> PercentileBinning;

  std::vector<Models::Function> FitFunctions;
  
  // Range scans
  std::vector<std::pair<double,double>> PTRanges;  // pt range for observable calc.
  std::vector<std::pair<double,double>> FitRanges; // pt range for spectra fit
  std::vector<std::pair<double,double>> SOSRanges; // x-axis range for power law fit
    
  // reference multiplicity window used for mean Nch
  int refLow  = 0;
  int refHigh = 999999;
  
  // Observable / spectrum settings
  double EtaAbsMax = 0.5;                // |eta| < ...
  double SpectraLowCut = 0.2;             // low pT cut for fit or data usage
  double SpectraHighCut = 10.0;           // high pT cut for fit
  double FitCutUseDataAbove = 0.2;        // your "fit_cut" threshold where you switch to data bins

  double acceptanceFactor = 1.0;
  bool useAcceptanceCorrection = true;

  bool useExternalYieldErrors = false;
  std::string externalYieldErrorDir = "";

  // Observables to compute/write
  std::vector<Observable> observables;

  // API
  static Config Defaults();
  void Validate() const;
  bool IsRunAllowed(int runId) const;

  static const char* CentralitySettingName(CentralitySetting s);
  static const char* PercFirstBinName(PercFirstBin p);
  static const char* ObservableName(Observable o);
};
