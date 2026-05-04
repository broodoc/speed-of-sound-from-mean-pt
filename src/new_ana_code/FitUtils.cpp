// Copyright 2026 Caleb Broodo
//GetFit, GetNParams, default parameters/limits, fit logging/statbox if you want)

#include "FitUtils.h"

#include <vector>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

#include "TMath.h"
#include "Math/Functor.h"

// ------------------------------------------------------------
// Internal: keep Models objects alive for TF1 functors.
// This avoids dangling captures when TF1 outlives local scope.
// ------------------------------------------------------------
namespace {
  std::vector<std::shared_ptr<Models>>& ModelKeepAlive() {
    static std::vector<std::shared_ptr<Models>> store;
    return store;
  }

  std::string AutoName(Models::Function f) {
    std::ostringstream os;
    os << "f_" << FitUtils::FunctionToString(f);
    return os.str();
  }
} // anonymous namespace

namespace FitUtils {

int GetNParams(Models::Function f) {
  switch (f) {
    case Models::Function::LevyTsallis:                    return 4;
    case Models::Function::Tsallis:                        return 3; // (your Tsallis used V,q,T)
    case Models::Function::Hagedorn:                       return 4;
    case Models::Function::TBW:                            return 5;
    case Models::Function::FokkerPlanck:                   return 6;
    case Models::Function::PowerLaw:                       return 2;
    case Models::Function::ExpDecay:                       return 3;
    case Models::Function::RightKneeFunction:              return 3; // cs2,knee,sigma in your code
    case Models::Function::RightKneeFunctionApproximation:  return 2; // (xx/par0)^par1 in your code
    default:                                               return 0;
  }
}

const char* FunctionToString(Models::Function f) {
  switch (f) {
    case Models::Function::LevyTsallis:                   return "LevyTsallis";
    case Models::Function::Tsallis:                       return "Tsallis";
    case Models::Function::Hagedorn:                      return "Hagedorn";
    case Models::Function::TBW:                           return "TBW";
    case Models::Function::FokkerPlanck:                  return "FokkerPlanck";
    case Models::Function::PowerLaw:                      return "PowerLaw";
    case Models::Function::ExpDecay:                      return "ExpDecay";
    case Models::Function::RightKneeFunction:             return "RightKneeFunction";
    case Models::Function::RightKneeFunctionApproximation:return "RightKneeFunctionApproximation";
    default:                                              return "UnknownFunction";
  }
}

void SetParNames(TF1* fit, Models::Function f) {
  if (!fit) return;

  switch (f) {
    case Models::Function::LevyTsallis:
      fit->SetParNames("m", "C", "n", "dNdy");
      break;
    case Models::Function::Tsallis:
      fit->SetParNames("V", "q", "T");
      break;
    case Models::Function::Hagedorn:
      fit->SetParNames("m", "pT0", "b", "Const");
      break;
    case Models::Function::TBW:
      fit->SetParNames("A", "m", "T", "q", "beta_s");
      break;
    case Models::Function::FokkerPlanck:
      fit->SetParNames("m_unused", "A", "b", "c", "d", "T");
      break;
    case Models::Function::PowerLaw:
      fit->SetParNames("exp", "scale");
      break;
    case Models::Function::ExpDecay:
      fit->SetParNames("p0", "p1", "p2");
      break;
    case Models::Function::RightKneeFunction:
      fit->SetParNames("cs2", "knee", "sigma");
      break;
    case Models::Function::RightKneeFunctionApproximation:
      fit->SetParNames("p0", "p1");
      break;
    default:
      break;
  }
}

void SetDefaultParameters(TF1* fit, Models::Function f) {
  if (!fit) return;

  switch (f) {
    case Models::Function::LevyTsallis:
      // m, C, n, dN/dy
      fit->SetParameters(0.139, 0.2, 12.0, 6000.0);
      break;

    case Models::Function::Tsallis:
      // V, q, T
      fit->SetParameters(1.0, 1.10, 0.15);
      break;

    case Models::Function::Hagedorn:
      // m, pT0, b, const
      fit->SetParameters(0.139, 0.2, 7.0, 100.0);
      break;

    case Models::Function::TBW:
      // A, m, T, q, beta_s
      fit->SetParameters(1.0, 0.139, 0.114, 1.086, 0.3);
      break;

    case Models::Function::FokkerPlanck:
      // m (unused in your operator for FP), A, b, c, d, T
      fit->SetParameters(0.139, 0.5, 1.0, 2.0, 2.0, 0.15);
      break;

    case Models::Function::PowerLaw:
      // exp, scale
      fit->SetParameters(0.2, 1.0);
      break;

    case Models::Function::ExpDecay:
      fit->SetParameters(1.0, 0.1, 2.0);
      break;

    case Models::Function::RightKneeFunction:
      // cs2, knee, sigma
      fit->SetParameters(0.2, 1.12, 0.06);
      break;

    case Models::Function::RightKneeFunctionApproximation:
      // p0, p1
      fit->SetParameters(1.0, 0.2);
      break;

    default:
      break;
  }
}

void SetDefaultLimits(TF1* fit, Models::Function f) {
  if (!fit) return;

  // These are intentionally conservative. Tune as you learn failure modes.
  switch (f) {
    case Models::Function::LevyTsallis:
      fit->SetParLimits(0, 0.0, 2.0);        // m
      fit->SetParLimits(1, 1e-4, 5.0);       // C
      fit->SetParLimits(2, 2.0, 100.0);      // n
      fit->SetParLimits(3, 1e-12, 1e12);     // dN/dy
      break;

    case Models::Function::Tsallis:
      fit->SetParLimits(0, 1e-12, 1e12);     // V
      fit->SetParLimits(1, 1.0001, 3.0);     // q (>1)
      fit->SetParLimits(2, 1e-4, 2.0);       // T
      break;

    case Models::Function::Hagedorn:
      fit->SetParLimits(0, 0.0, 2.0);        // m
      fit->SetParLimits(1, 1e-4, 10.0);      // pT0
      fit->SetParLimits(2, 0.0, 200.0);      // b
      fit->SetParLimits(3, 1e-12, 1e12);     // const
      break;

    case Models::Function::TBW:
      fit->SetParLimits(0, 1e-12, 1e12);     // A
      fit->SetParLimits(1, 0.0, 2.0);        // m
      fit->SetParLimits(2, 1e-4, 2.0);       // T
      fit->SetParLimits(3, 1.0001, 3.0);     // q
      fit->SetParLimits(4, 0.0, 0.999);      // beta_s
      break;

    case Models::Function::FokkerPlanck:
      // m unused — but keep finite
      fit->SetParLimits(0, 0.0, 2.0);
      fit->SetParLimits(1, 1e-12, 1e12);     // A
      fit->SetParLimits(2, 0.0, 100.0);      // b
      fit->SetParLimits(3, 0.0, 100.0);      // c
      fit->SetParLimits(4, 0.0, 100.0);      // d
      fit->SetParLimits(5, 1e-4, 2.0);       // T
      break;

    case Models::Function::PowerLaw:
      fit->SetParLimits(0, -10.0, 10.0);     // exponent
      fit->SetParLimits(1, 1e-12, 1e12);     // scale
      break;

    case Models::Function::ExpDecay:
      fit->SetParLimits(0, -1e12, 1e12);
      fit->SetParLimits(1, -1e12, 1e12);
      fit->SetParLimits(2, 1e-6, 1e6);       // decay rate in exp(par2 * x)
      break;

    case Models::Function::RightKneeFunction:
      fit->SetParLimits(0, -5.0, 5.0);       // cs2-ish exponent
      fit->SetParLimits(1, 0.5, 2.0);        // knee (typical)
      fit->SetParLimits(2, 1e-4, 1.0);       // sigma
      break;

    case Models::Function::RightKneeFunctionApproximation:
      fit->SetParLimits(0, 1e-12, 1e12);
      fit->SetParLimits(1, -10.0, 10.0);
      break;

    default:
      break;
  }
}

TF1* MakeTF1(const std::string& name,
             Models::Function f,
             double xmin,
             double xmax,
             bool applyDefaults,
             bool applyLimits) {
  const int npar = GetNParams(f);
  if (npar <= 0) return nullptr;

  // Keep the model alive and bind it into a ROOT functor.
  auto model = std::make_shared<Models>(f);
  ModelKeepAlive().push_back(model);

  // ROOT::Math::Functor expects a callable (double const*, double const*) -> double
  ROOT::Math::Functor functor(
    [model](const double* x, const double* p) -> double {
      // Models::operator() uses (double*, double*)
      return (*model)(const_cast<double*>(x), const_cast<double*>(p));
    },
    npar
  );

  TF1* fit = new TF1(name.c_str(), functor, xmin, xmax, npar);
  fit->SetNpx(2000); // can be overwritten by Config later if you want
  SetParNames(fit, f);

  if (applyDefaults) SetDefaultParameters(fit, f);
  if (applyLimits)   SetDefaultLimits(fit, f);

  return fit;
}

TF1* MakeTF1(Models::Function f,
             double xmin,
             double xmax,
             bool applyDefaults,
             bool applyLimits) {
  return MakeTF1(AutoName(f), f, xmin, xmax, applyDefaults, applyLimits);
}

TFitResultPtr FitGraph(TGraphErrors* g,
                       TF1* fit,
                       const std::string& options) {
  if (!g || !fit) return TFitResultPtr();
  return g->Fit(fit, options.c_str());
}

void PrintFitSummary(const TF1* fit, const char* label) {
  if (!fit) return;

  const double chi2 = fit->GetChisquare();
  const int    ndf  = fit->GetNDF();
  const int    npar = fit->GetNpar();

  if (label) std::cout << label << "\n";
  std::cout << "  TF1: " << fit->GetName()
            << "  chi2/ndf = ";
  if (ndf > 0) std::cout << (chi2 / ndf);
  else         std::cout << "NaN";
  std::cout << "  (chi2=" << chi2 << ", ndf=" << ndf << ")"
            << "  npar=" << npar
            << "\n";
}

} // namespace FitUtils
