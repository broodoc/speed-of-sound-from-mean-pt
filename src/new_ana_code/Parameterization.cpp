#include "Parameterization.h"

#include "Config.h"
#include "OutputManager.h"

// Modules you said you have:
#include "Models.h"
#include "FitUtils.h"
#include "SpectraObservables.h"
#include "ErrorPropagation.h"
#include "CentralityBinning.h"

#include <stdexcept>
#include <iostream>

// ROOT forward-declared in header; include in .cpp
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TList.h"

Parameterization::Parameterization(const Config& cfg, OutputManager& out)
: cfg_(cfg), out_(out)
{
  // Prefer ranges in cfg; but if you still want defaults here:
  pt_ranges_  = cfg_.ptRanges;   // if you add these to config
  fit_ranges_ = cfg_.fitRanges;
  sos_ranges_ = cfg_.sosRanges;
}

const std::vector<std::pair<double,double>>& Parameterization::GetPtRanges() const  { return pt_ranges_; }
const std::vector<std::pair<double,double>>& Parameterization::GetFitRanges() const { return fit_ranges_; }
const std::vector<std::pair<double,double>>& Parameterization::GetSoSRanges() const { return sos_ranges_; }

void Parameterization::ValidateInputs_(TH2D* Pt_Nch_hist, TH1D* Nch_dist) const {
  if (!Pt_Nch_hist) throw std::runtime_error("Pt_Nch_hist is null");
  if (!Nch_dist)    throw std::runtime_error("Nch_dist is null");
}

void Parameterization::Ratio_compute_pt_Nch(TH2D* Pt_Nch_hist, TH1D* Nch_dist) {
  ValidateInputs_(Pt_Nch_hist, Nch_dist);

  // OutputManager owns “where” things go; Parameterization just asks it to prep.
  out_.BeginFile();          // optional no-op if file already open
  out_.WriteConfig(cfg_);    // store cfg into the ROOT output

  RunFullScan_(Pt_Nch_hist, Nch_dist);

  out_.FinalizeFile();       // optional: write lists, close dirs, flush
}

void Parameterization::RunFullScan_(TH2D* Pt_Nch_hist, TH1D* Nch_dist) {
  BestFitTracker best;

  // Example: global objects that are written once.
  auto globalList = out_.GetOrMakeList("global");

  // You likely have “MakePercentiles” logic; that should be in CentralityBinning now.
  // Example usage:
  // CentralityBinning::FindPercentiles(Nch_dist, ...);

  for (const auto& fit_range : fit_ranges_) {
    for (const auto& pt_range : pt_ranges_) {

      // 1) Build reference spectrum for this pT range (your code)
      //    This part is analysis-specific: you may build TH1D* from Pt_Nch_hist projections.
      //    Keep the *construction* here, but everything math/normalization should be in SpectraObservables.
      TH1D* refSpec = SpectraObservables::MakeRefSpectrum(Pt_Nch_hist, /*...*/, pt_range);

      for (auto f : cfg_.fitFunctions) { // vector<Parameterization::Function>

        // 2) Graph container for this combo
        auto graphs = MakeGraphs_(fit_range, pt_range, f);

        // 3) Fit reference spectrum (FitUtils owns TF1 building + parameter seeds)
        FitUtils::FitOutput refFitOut = FitUtils::FitSpectrum(
          *refSpec,
          f,
          /*low_cut=*/cfg_.spectraLowCut,
          /*high_cut=*/cfg_.spectraHighCut,
          cfg_
        );

        // 4) Loop percentile bins (example)
        const int nperc = CentralityBinning::GetNumBins(cfg_);
        for (int k = 0; k < nperc; ++k) {

          TH1D* spec = SpectraObservables::MakeSpectrumForBin(Pt_Nch_hist, k, /*...*/, pt_range);

          // Optionally apply per-bin errors from external file
          if (cfg_.useExternalYieldErrors) {
            ErrorPropagation::SetErrorYield(*spec, k, cfg_);
          }

          // Fit per-bin spectrum
          FitUtils::FitOutput fitOut = FitUtils::FitSpectrum(
            *spec,
            f,
            cfg_.spectraLowCut,
            cfg_.spectraHighCut,
            cfg_
          );

          // 5) Compute observables from spec+fit (mean pT, dNdEta, Nch, etc.)
          const double meanpt = SpectraObservables::GetMean(*spec, *fitOut.fit, cfg_);
          const double dndeta = SpectraObservables::GetdNdEta(*spec, *fitOut.fit, cfg_);
          const double nch    = SpectraObservables::GetNch(*Nch_dist, cfg_.refLow, cfg_.refHigh);

          // 6) Statistical errors (data-only, fit-only, etc.)
          const double meanptErr = ErrorPropagation::GetErrorMeanpT(*spec, *fitOut, cfg_);
          const double dndetaErr = ErrorPropagation::GetErrordNdEta(*spec, *fitOut, cfg_);

          // 7) Fill graphs
          SpectraObservables::FillGraphs(graphs, k, nch, dndeta, meanpt, dndetaErr, meanptErr, cfg_);

          // cleanup per-bin
          delete spec;
        }

        // 8) Run c_s^2 fits (FitUtils or a dedicated SoS module)
        FitUtils::SoSFitOutput sosOut = FitUtils::FitSoundSpeed(
          *graphs.gMeanPt, *graphs.gNch,
          fit_range.first, fit_range.second,
          cfg_
        );

        // 9) Update best-fit tracker
        if (sosOut.isValid && sosOut.chi2ndf < best.bestChi2Ndf) {
          best.hasBest = true;
          best.bestChi2Ndf = sosOut.chi2ndf;
          best.bestCs2 = sosOut.cs2;
          best.bestCs2Err = sosOut.cs2Err;
          best.bestPtRange = pt_range;
          best.bestFitRange = fit_range;
          best.bestFunction = f;
        }

        // 10) Write combo outputs
        out_.WriteCombo(fit_range, pt_range, f, graphs, sosOut);

        // 11) Cleanup combo
        WriteAndCleanup_(graphs);

      } // func loop

      delete refSpec;
    }
  }

  // Write global summary (best fit, etc.)
  if (best.hasBest) {
    out_.WriteBestFitSummary(best.bestCs2, best.bestCs2Err, best.bestChi2Ndf,
                             best.bestFitRange, best.bestPtRange,
                             static_cast<int>(best.bestFunction));
  }
}

Parameterization::GraphBundle
Parameterization::MakeGraphs_(const std::pair<double,double>& fitRange,
                              const std::pair<double,double>& ptRange,
                              Function f)
{
  GraphBundle gb;

  // OutputManager decides directory structure + owns lists.
  gb.list = out_.GetOrMakeList(out_.ComboKey(fitRange, ptRange, f).c_str());

  gb.gMeanPt = SpectraObservables::MakeGraph("gMeanPt", "<pT> vs Nch");
  gb.gNch    = SpectraObservables::MakeGraph("gNch",    "Nch (ref) vs bin");
  gb.gSoS    = SpectraObservables::MakeGraph("gSoS",    "cs2 extraction");

  // Optionally register graphs into list now
  if (gb.list) {
    gb.list->Add(gb.gMeanPt);
    gb.list->Add(gb.gNch);
    gb.list->Add(gb.gSoS);
  }

  return gb;
}

void Parameterization::WriteAndCleanup_(GraphBundle& gb) {
  // OutputManager can write lists itself; Parameterization just deletes local allocations if needed.
  // If you transfer ownership to OutputManager/list, do NOT delete here.
  // Decide one ownership model and stick to it.

  // Example: if graphs are owned by list, do nothing:
  gb.gMeanPt = nullptr;
  gb.gNch = nullptr;
  gb.gSoS = nullptr;
  gb.list = nullptr;
}
