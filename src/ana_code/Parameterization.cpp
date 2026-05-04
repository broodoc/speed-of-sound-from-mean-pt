// Parameterization - Bulk Speed of sound analysis functions developed by Caleb Broodo
#include "Parameterization.h"
#include "MomentHistogram.h"

Parameterization::Parameterization()
    : type_f(Function::LevyTsallis),
      fTrigger("V0_CT"),
      fMultEstimator("TPC"),
      fSystem("AuAu_3"),
      fIsMC(false),
      fIsMCTrue(false) {}
Parameterization::Parameterization(Parameterization::Function type)
    : type_f(type),
      fTrigger("V0_CT"),
      fMultEstimator("TPC"),
      fSystem("AuAu_3"),
      fIsMC(false),
      fIsMCTrue(false) {}

void Parameterization::FindPercentiles(TH1 *hin,
                                       PercFirstBin startBin, bool debug,
                                       bool isCoarseBinning) {
    string name = hin->GetName();
    if (name.length() >= 4) {
        // Replace the last four characters with "bin_cuts.txt"
        name.replace(name.length() - 4, 4, "bin_cuts");
    }
    std::string output_label = name + ".txt";
    std::ofstream outFile(output_label);
    std::string outputError_label = name + "Error.txt";
    std::ofstream outFileError(outputError_label);
    
    int nBins{50};
  if (isCoarseBinning) {
    nBins = this->GetnCoarsePercentileBinning();
      cout << "Getting coarse percentile binning for " << output_label << "... " << '\n';
  } else {
    nBins = this->GetnFinePercentileBinning();
      cout << "Getting fine percentile binning for " << output_label << "... " << '\n';
  }
  double newBin[nBins + 1];
  double errorBin[nBins + 1];
  for (int bin_y_perc = 0; bin_y_perc < nBins + 1; ++bin_y_perc) {
      newBin[bin_y_perc] = 0; errorBin[bin_y_perc] = 0;
  }
  const int maxBinX{hin->GetNbinsX()};
  const double integral{hin->Integral(1, maxBinX)};

  if (integral <= 1E3) {
    cout << "Less than 1E3 entries, not enough statistics" << endl;
  }
    if (debug){cout << "	Total number of bins = " << maxBinX << '\n'
        << "	Number of percentile bins = " << nBins << '\n'
        << " 	Total number of events (total entries in the histogram) = "
        << integral << '\n';
    }
  int foundbin{-1};
  int firstbin{-1};

  /// Find percentiles from the first bin
  if (startBin == PercFirstBin::kFirst) {
    firstbin = 1;
    for (int bin_y_perc = 0; bin_y_perc < nBins; ++bin_y_perc) {
      if (foundbin > 0) {
        firstbin = foundbin + 1;
      }

      float high_perc{0};
      float low_perc{0};
      if (isCoarseBinning) {
        low_perc = this->GetCoarseAmplitudeCutLow(bin_y_perc);
        high_perc = this->GetCoarseAmplitudeCutHigh(bin_y_perc);
      } else {
        low_perc = this->GetFineAmplitudeCutLow(bin_y_perc);
        high_perc = this->GetFineAmplitudeCutHigh(bin_y_perc);
      }

      double intergral_tmp{0.0};
      double events{0.0};
      for (int bin_y = firstbin; bin_y <= maxBinX; ++bin_y) {
        intergral_tmp += hin->GetBinContent(bin_y) / integral;
        events += hin->GetBinContent(bin_y);
        if (intergral_tmp > abs(high_perc - low_perc)) {
          foundbin = bin_y;
          break;
        }
      }
      // integ_classes += intergral_tmp;
      // cout << "bin_y_perc=" << bin_y_perc << "  integ_classes=" <<
      // integ_classes << endl; cout << "Perc. [" << RhoPercBining[bin_y_perc]
      // << "-" << RhoPercBining[bin_y_perc+1] << "]" << " foundbin=" <<
      // foundbin << endl;
      if (bin_y_perc == nBins - 1) {
        foundbin = maxBinX;
      }
      newBin[bin_y_perc] = firstbin;

      double low_edge{hin->GetXaxis()->GetBinLowEdge(firstbin)};
      double up_edge{hin->GetXaxis()->GetBinUpEdge(foundbin)};

        if (debug) {cout << " Real Percentage: " << 100 * abs(high_perc - low_perc) << "% |"
            << " Amplitude range: " << low_edge << " - " << up_edge << " |"
            << " Measured percentage: " << 100 * intergral_tmp << "% |"
            << " Number of events: " << events << '\n';
            
            cout << " Percentage window: " << 100 * low_perc << " - "
            << 100 * high_perc << " % |" << '\n'
            << "======================="
            << "======================="
            << "======================="
            << "======================="
            << "=======================" << '\n';
        }
    }
    newBin[nBins] = maxBinX;
  }

  //! Find percentiles from the last bin
  if (startBin == PercFirstBin::kLast) {
    firstbin = maxBinX;
    for (int bin_y_perc = 0; bin_y_perc < nBins; ++bin_y_perc) {
      if (foundbin > 0) {
        firstbin = foundbin - 1;
//        firstbin = foundbin;
      }

      double high_perc{0};
      double low_perc{0};
      if (isCoarseBinning) {
        low_perc = this->GetCoarseAmplitudeCutLow(bin_y_perc);
        high_perc = this->GetCoarseAmplitudeCutHigh(bin_y_perc);
        // high_perc = CoarsePercentileBinning[bin_y_perc + 1];
        // low_perc = CoarsePercentileBinning[bin_y_perc];
      } else {
        low_perc = this->GetFineAmplitudeCutLow(bin_y_perc);
        high_perc = this->GetFineAmplitudeCutHigh(bin_y_perc);
      }

      double intergral_tmp{0.0};
      double events{0.0};
      for (int bin_y = firstbin; bin_y >= 1; bin_y--) {
        intergral_tmp += hin->GetBinContent(bin_y) / integral;
        events += hin->GetBinContent(bin_y);
          if (high_perc > 0.99) {foundbin = 1; break;}
        if (intergral_tmp > abs(high_perc - low_perc)) {
          foundbin = bin_y; //one bin element above the percentile class
          break;
        }
      }

      if (bin_y_perc == (nBins - 1) && foundbin <= 0) {
          foundbin = 1;
      }
//      if (bin_y_perc == (nBins - 1)) {
//        foundbin = 1;
//      }

      newBin[bin_y_perc] = firstbin;
      errorBin[bin_y_perc] = abs ( (intergral_tmp - abs(high_perc - low_perc)) / abs(high_perc - low_perc) );

      double up_edge{hin->GetXaxis()->GetBinUpEdge(firstbin)};
      double low_edge{hin->GetXaxis()->GetBinLowEdge(foundbin)};
        
        if (up_edge <= low_edge) {
            std::cerr << "Degenerate bin range for percentile bin " << bin_y_perc
            << " ( " << low_edge << " - " << up_edge << " ) "
            << std::endl;
            continue;
        }

        if (debug) {cout << " Real Percentage: " << 100 * abs(high_perc - low_perc) << "% |"
            << " N_{charged} range: " << low_edge << " - " << up_edge << " |"
            << " Measured percentage: " << 100 * intergral_tmp << "% |"
            << " Number of events: " << events << '\n';
            
            cout << " Percentage window: " << 100 * low_perc << " - "
            << 100 * high_perc << " % | " << " Relative percentage error: " << errorBin[bin_y_perc] << '\n'
            << "======================="
            << "======================="
            << "======================="
            << "======================="
            << "=======================" << '\n';
        }
    }
//    errorBin[nBins] = 0;
    newBin[nBins] = 1;
  }

  /// Except from the first and last bin cut, the lower bin cut has to start one
  /// bin ahead of what is given by these bin cuts.
  cout << "	THESE ARE THE BIN CUTS: ";
    if (outFile.is_open()){
        for (int bin_y_perc = 0; bin_y_perc < nBins + 1; ++bin_y_perc) {
            cout << newBin[bin_y_perc] << ", ";
            outFile << newBin[bin_y_perc] << ", ";
        }
        outFile.close();
        std::cout << '\n' << "Output file processed for " << output_label << std::endl;
    }
    else {
            std::cerr << "Unable to open file" << std::endl;
        }
  cout << endl;
    
    cout << "    THESE ARE THE RELATIVE PERCENTAGE ERRORS: ";
      if (outFileError.is_open()){
          for (int bin_y_perc = 0; bin_y_perc < nBins + 1; ++bin_y_perc) {
              cout << errorBin[bin_y_perc] << ", ";
              outFileError << errorBin[bin_y_perc] << ", ";
          }
          outFileError.close();
          std::cout << '\n' << "Output file processed for " << outputError_label << std::endl;
      }
      else {
              std::cerr << "Unable to open file" << std::endl;
          }
    cout << endl;
  
}

//-----------------------------------------------------

double Parameterization::operator()(double *x, double *par) {
    double xx = x[0];
    switch (type_f) {
            // Type-A Tsallis
            // https://arxiv.org/pdf/1905.08623.pdf
            // m - par[0]
            // C - par[1]
            // n - par[2]
            // dN/dy - par[3]
        case Function::LevyTsallis:
            return par[3] * (par[2] - 1) * (par[2] - 2) /
            (2 * 3.14159 * par[1] * par[2] *
             (par[1] * par[2] + par[0] * (par[2] - 2))) *
            pow(1 + (sqrt(xx * xx + par[0] * par[0]) - par[0]) /
                (par[1] * par[2]),
                -par[2]) *
            xx;
            break;
        case Function::Tsallis:
            // https://arxiv.org/abs/2008.02559
            // V - par[0]
            // q - par[1]
            // T - par[2]
            return (2.0 * pow(xx, 2.0) * par[0] / 39.478418) *
            (pow((1.0 +
                  (par[1] - 1) *
                  (sqrt(pow(0.139570, 2.0) + pow(xx, 2.0)) / par[2])),
                 (-1.0 * par[1] / (par[1] - 1.0))) +
             pow((1.0 +
                  (par[1] - 1) *
                  (sqrt(pow(0.493677, 2.0) + pow(xx, 2.0)) / par[2])),
                 (-1.0 * par[1] / (par[1] - 1.0))) +
             2.0 * pow((1.0 + (par[1] - 1) *
                        (sqrt(pow(0.938272, 2.0) + pow(xx, 2.0)) /
                         par[2])),
                       (-1.0 * par[1] / (par[1] - 1.0))));
            break;
        case Function::Hagedorn:
            // https://arxiv.org/pdf/1007.0719.pdf
            // m - par[0]
            // pT0 - par[1]
            // b - par[2]
            // constant - par[3]
            return par[4] * (xx * xx / sqrt(par[0] * par[0] + xx * xx)) *
                        pow((1 + xx / par[1]), -par[2]) * par[3];
            // https://drupal.star.bnl.gov/STAR/system/files/meanpT_Coll_March2024_0.pdf slide 10
            // m - par[0]
            // beta_T - par[1]
            // n - par[2]
            // T - par[3]
            // A - par[4]
//            return par[4] * xx * pow(1 + ( 1/sqrt(1 - par[1] * par[1]) )
//                                         * ( sqrt( xx * xx + par[0]*par[0] ) - ( xx * par[1] ) )
//                                        / ( par[2] * par[3] ), -par[2] );
            break;
        case Function::BW: {
            // Blast-wave (Boltzmann-Gibbs) parameterization
            // https://arxiv.org/pdf/nucl-th/0312024
            // https://arxiv.org/pdf/2409.19758v2
            // m      - par[0]
            // beta_T - par[1]   (interpreted as surface velocity beta_s)
            // n      - par[2]   (profile exponent in beta(r)=beta_s*(r/R)^n)
            // T      - par[3]
            // A      - par[4]
                const double xx = x[0];           // pT
                const double m  = par[0];
                const double beta_s = par[1];
                const double n  = par[2];
                const double T  = par[3];
                const double A  = par[4];

                if (xx < 0.0) return 0.0;
                if (T <= 0.0) return 0.0;
                if (beta_s < 0.0) return 0.0;
                if (beta_s >= 1.0) return 0.0;   // avoid atanh(>=1)
                if (n < 0.0) return 0.0;

                const double mT = std::sqrt(xx*xx + m*m);

                // Numerical integral over xi = r/R in [0,1]
                // Using Simpson's rule
                const int N = 200;                 // must be even
                const double a = 0.0, b = 1.0;
                const double h = (b - a) / N;

                auto integrand = [&](double xi) -> double {
                    // beta(xi) = beta_s * xi^n
                    double beta = beta_s * std::pow(xi, n);
                    if (beta >= 0.999999999) beta = 0.999999999;  // safety clamp

                    const double rho = 0.5 * std::log((1.0 + beta) / (1.0 - beta)); // atanh(beta)

                    const double argI0 = (xx * std::sinh(rho)) / T;
                    const double argK1 = (mT * std::cosh(rho)) / T;

                    // ROOT provides these:
                    const double I0 = TMath::BesselI0(argI0);
                    const double K1 = TMath::BesselK1(argK1);

                    return xi * I0 * K1;
                };

                double S = integrand(a) + integrand(b);
                for (int i = 1; i < N; ++i) {
                    const double xi = a + i*h;
                    S += (i % 2 ? 4.0 : 2.0) * integrand(xi);
                }
                const double integral = (h / 3.0) * S;

                // dN/mTdmT form: A * pT * mT * integral
                // dN/pT form: A * pT * integral
                return A * xx * mT * integral;
        }
        case Function::TBW: {
            // Tsallis Blast-Wave model (https://arxiv.org/pdf/0812.1609)
            // par[0] = A (normalization)
            // par[1] = m (particle mass, GeV) -> TOTAL_CHARGED_MASS
            // par[2] = T (temperature, GeV)
            // par[3] = q (nonextensivity parameter, q > 1)
            // par[4] = beta_s (surface flow velocity)

            const double A      = par[0];
            const double m      = par[1];
            const double T      = par[2];
            const double q      = par[3];
            const double beta_s = par[4];
            const double pT     = x[0];
            if (T <= 0.0 || q <= 1.0 || beta_s < 0.0 || beta_s >= 1.0) return 0.0;

            // ----- TBW integrand -----
//            const int n_r   = 24;
//            const int n_phi = 36;
//            const int n_y   = 24;
//            const int n_r   = 10;
//            const int n_phi = 10;
//            const int n_y   = 10;
            const int n_r   = 4;
            const int n_phi = 5;
            const int n_y   = 5;
            const double R  = 1.0;
            const double dr   = R / n_r;
            const double dphi = 2.0 * TMath::Pi() / n_phi;
            const double Ywin = 0.10; // rapidity half-window: integrates y ∈ [-Ywin, +Ywin]
            const double dy   = (2.0 * Ywin) / n_y;

            const double mT = std::sqrt(pT * pT + m * m);
            double sum = 0.0;

            for (int ir = 0; ir < n_r; ++ir) {
                const double r  = (ir + 0.5) * dr;
                const double br = beta_s * (r / R);
                if (br >= 0.999999) continue;
                const double rho = std::atanh(br);

                for (int ip = 0; ip < n_phi; ++ip) {
                    const double phi = -TMath::Pi() + (ip + 0.5) * dphi;
                    const double cosphi = std::cos(phi);

                    for (int iy = 0; iy < n_y; ++iy) {
                        const double y   = -Ywin + (iy + 0.5) * dy;
                        const double cy  = std::cosh(y);

                        const double expo = (mT * cy * std::cosh(rho)
                                           - pT * std::sinh(rho) * cosphi) / T;
                        const double base = 1.0 + (q - 1.0) * expo;
                        if (base <= 0.0) continue;

                        const double weight = std::pow(base, -1.0 / (q - 1.0));
                        sum += r * cy * weight;
                    }
                }
            }

            const double shape = mT * sum * dr * dphi * dy;
//            std::cout << "A * shape: " << A * shape << std::endl;
            return A * shape * pT;
            break;
        }
        case Function::FokkerPlanck: {
//                    double m = par[0], A = par[1], b = par[2], c = par[3], d = par[4], T = par[5];
                    double A = par[1], b = par[2], c = par[3], d = par[4], T = par[5];
                    double arg1 = xx + b;
                    double arg2 = -d * xx / T;
                    if (arg1 <= 0 || T <= 0 || !std::isfinite(arg1) || !std::isfinite(arg2))
                        return 1e-12;
                    double result = A * std::pow(arg1, -c) * std::exp(arg2);
                    return (std::isfinite(result) && result >= 0) ? result : 1e-12;
                }
        case Function::PowerLaw:
            // n - par[0]
//            return (par[1] / pow(xx, par[0]));
//            return par[1] * pow(xx, par[0]);
            return par[1] * pow(xx / par[1], par[0]);
            break;
        case Function::ExpDecay:
            // n - par[0] //Boltzmann
            return xx * (par[0] + par[1] * xx * xx) / exp(par[2] * xx);
            break;
        case Function::RightKneeFunction:{
            double n = x[0];
            double cs2 = par[0];
            double knee = par[1];
            double sigma = par[2];
//            double sigma = 0.0592; //par2
//            double knee = 1.1168; //par1

            double norm = sigma * TMath::Sqrt(2 / TMath::Pi());
            double exponent = TMath::Exp(-TMath::Power(n - knee, 2) / (2 * TMath::Power(sigma, 2)));
            double erfc = TMath::Erfc((n - knee) / (TMath::Sqrt(2) * sigma));
            if (erfc < 1e-10) return 1e9;
            
            double f = n - (norm * exponent / erfc);
            
            
            if (f <= 1e-5) return 1e9;
            
            
            double value = exp(cs2 * log(n / f));
//            double value = pow( ( n / f ), cs2);
            return value;
            break;
        }
        case Function::RightKneeFunctionApproximation:
            return TMath::Power( (xx/par[0]), par[1]);
            break;
    }
}
//double value = pow( ( xx / Parameterization::f(xx) ), par[0]);
//-----------------------------------------------------

Double_t Parameterization::RightKneeFunction(Double_t *x, Double_t *par){
    Double_t xx = x[0];
    return TMath::Power( (xx / par[0]), par[1]);
}

//-----------------------------------------------------

Double_t Parameterization::PowerLaw(Double_t *x, Double_t *par){
    Double_t xx = x[0];
    return par[1] * pow(xx, par[0]);
}

//-----------------------------------------------------

TF1* Parameterization::GetFit(TH1* h, double low_cut, double high_cut, Parameterization::Function fitFuncType, const TString& name, int index) {
//    std::cout << "Inside GetFit()" << std::endl;
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(111);
//    bool debug{false};
    freopen("/dev/null", "w", stderr);  // Silence ROOT warnings
    gErrorIgnoreLevel = kWarning;
    h->GetXaxis()->SetRangeUser(low_cut, high_cut);
    
    Parameterization model(fitFuncType);
    int nParams = GetNParams(fitFuncType);

    TF1* fit = new TF1(name, model, low_cut, high_cut, nParams);
    fit->SetNpx(50000);
    fit->SetLineColor(kRed);
    fit->SetLineWidth(6.0);
    string fitname = "";
    double TOTAL_CHARGED_MASS = 0.7*0.13957039 + 0.14*0.493677 + 0.16*0.938272081; //Estimated total charged mass (0.7*m_pion + 0.14*m_kaon + 0.16*M_proton)

    
    switch (fitFuncType) {
        
        case Function::LevyTsallis:
//            fit->SetParameters(5e-3, 0.158, 7.9, 5500);
            fit->SetParameters(0.139, 0.158, 7.9, 5500);
            fit->FixParameter(0, 0.139);
//            fit->FixParameter(0, 0.01);
//            fit->FixParameter(0, 0.);
//            fit->FixParameter(0, 0.12);
//            fit->FixParameter(0, 4e-07);
//            fit->FixParameter(0, 0.5e-01);
//            fit->FixParameter(0, 0.9e-01);
//            fit->SetParLimits(0, 0.08, 1.139);
//            fit->FixParameter(0, 0.);
//            fit->FixParameter(0, -0.1);
//            fit->FixParameter(0, 5e-6);
//            fit->FixParameter(0, 0.139);
//            fit->SetParLimits(3, 2000, 10e4);
//            fit->SetParLimits(3, 2000, 9000);
//            fit->SetParLimits(3, 2000, 9500);

            fit->SetParNames("m", "C", "n", "dN/dy");
            fitname = "LevyTsallis fit";
            break;
        case Function::BW:
            fit->SetParameters(
                0.1,   // m  (GeV)
                0.7,     // beta_s  (start ~0.6–0.8)
                1.0,     // n  (profile exponent)
                0.1,     // T  (GeV)
                1000.0   // A  (normalization guess)
            );

            fit->SetParNames("m", "beta_s", "n", "T", "A");

            fit->FixParameter(0, 0.139);
//            fit->FixParameter(0, TOTAL_CHARGED_MASS);
//            fit->FixParameter(0, 0.01);
            fit->SetParLimits(1, 0.0, 0.99);  // beta_s
            fit->SetParLimits(2, 0.0, 5.0);   // n
            fit->SetParLimits(3, 0.01, 0.3);  // T
            fit->SetParLimits(4, 0.0, 1e8);   // A

            fitname = "Blast-Wave fit";
            break;
        case Function::Tsallis:
            fit->SetParameters(10000, 1.1, 0.1); // V, q, T
            
            fit->SetParNames("V", "q", "T", "");
            fitname = "Tsallis fit";
            break;
        case Function::Hagedorn:
            fit->SetParameters(0.139, 0.11, 7.5, 1000, 3000);
//            fit->FixParameter(0, 0.139);                       // IMPORTANT (see next section)
            fit->FixParameter(0, TOTAL_CHARGED_MASS);                       // IMPORTANT (see next section)
//            fit->FixParameter(0, 1e-01);                       // IMPORTANT (see next section)
//            fit->SetParLimits(1, 0.0, 0.95);                   // beta must be < 1
//            fit->SetParLimits(2, 2.0, 50.0);                   // n reasonable
//            fit->SetParLimits(3, 0.02, 0.5);                   // T in GeV
            fit->SetParLimits(4, 0.0, 1e7);                    // A scale
            
            
            fit->SetParNames("m", "beta_T", "n", "T", "norm.");
//            fit->SetParNames("m", "pT0", "b", "const.");
            fitname = "Hagedorn fit";
            break;
        case Function::FokkerPlanck:
            fit->SetParameters(0.139, 0.1, 0.5, 1.5, 2.0, 0.1);
            
            fit->SetParNames("m", "A", "b", "c", "d", "T");
            fitname = "FokkerPlanck fit";
            break;
        case Function::PowerLaw:
            fit->SetParameters(6.0, 1000);
            
            fit->SetParNames("n", "Norm");
            fitname = "PowerLaw fit";
            break;
        case Function::ExpDecay:
            //Good for 0-0.35 w/ around -
            //                        fit->SetParLimits(0, 0, 1.9e4);
            //                        fit->SetParLimits(1, 0, 1e6);
//                                    fit->SetParLimits(2, 7.1, 7.3);
                fit->SetParameters(1.916e4, 8.375e4, 7.1);
//                fit->SetParLimits(0, 500, 2e4);
//                fit->SetParLimits(1, 8.76e4, 9e4);
//            if (index < 0) fit->SetParLimits(1, 8.2e4, 8.685e4);
                
//                fit->SetParLimits(2, 0, 2.74);


            fit->SetParNames("Const", "Slope", "Decay");
            fitname = "ExpDecay fit";
            break;
        case Function::TBW:
            // based on /Users/calebbroodo/aux_sos_analysis/embedding/RCF/200GeV/PID_new/TBW_curves.cpp
            fit->SetParameters(0, 500000.0);  // A is a free parameter
            fit->FixParameter(1, TOTAL_CHARGED_MASS); //m
            fit->FixParameter(2, 0.07); //T
            fit->FixParameter(3, 1.086); //q
            fit->SetParameters(4, 0.0);   //beta_s
            fit->SetParLimits(4, 0.3, 1.); //q
            fit->SetParNames("A", "m", "T", "q", "beta_s");
            fitname = "TBW fit";
            break;
        case Function::RightKneeFunction:
            fit->SetParameters(1.0);
            fit->SetParNames("Slope");
            break;
        case Function::RightKneeFunctionApproximation:
            fit->SetParameters(1.0, 0.5, 1.0);
            fit->SetParNames("Const", "Slope", "Decay");
            break;
    }

    if (!lastSeedParams_.empty() && (int)lastSeedParams_.size() == fit->GetNpar() && index > 1) {
        fit->SetParameters(lastSeedParams_.data());
        
        switch (fitFuncType) {
            case Function::LevyTsallis: {
                const double upperC = lastSeedParams_[1];
                const double uppern = lastSeedParams_[2];
//                fit->SetParLimits(1, 0.125, upperC);
//                fit->FixParameter(0, 1e-01);
//                fit->SetParLimits(1, 0.15, upperC);
//                fit->SetParLimits(2, 7.0, uppern);
//                fit->SetParLimits(2, 5.0, uppern);
//                fit->SetParLimits(2, 4.0, uppern);
                break;
            }
            case Function::Hagedorn:{
//                const double upperm = lastSeedParams_[0];
//                const double lowerbeta = lastSeedParams_[1];
//                const double lowerM = lastSeedParams_[2];
//                fit->FixParameter(0, 0.06);
//                fit->SetParLimits(0, 5e-6, upperm);
//                fit->SetParLimits(1, lowerbeta, 2.);
//                fit->SetParLimits(2, lowerM, 10.);
                break;
            }
            case Function::BW:{
                fit->SetParNames("m", "beta_s", "n", "T", "A");
                const double lowerbeta = lastSeedParams_[1];
                const double uppern = lastSeedParams_[2];
                const double upperT = lastSeedParams_[3];
                const double upperA = lastSeedParams_[4];
                fit->SetParLimits(1, lowerbeta, 0.999);
                fit->SetParLimits(2, 0.0, uppern);
                fit->SetParLimits(3, 0.01, upperT);
                fit->SetParLimits(4, 0.0, upperA);
                break;
            }
            case Function::Tsallis:
                break;
            case Function::FokkerPlanck:
                break;
            case Function::PowerLaw:
                break;
            case Function::ExpDecay:{
//                const double upperConst = lastSeedParams_[0];
                const double lowerSlope = lastSeedParams_[1];
                const double lowerDecay = lastSeedParams_[2];
//                fit->SetParLimits(1, lowerSlope, 1.2e5);
//                fit->SetParLimits(2, lowerDecay, 9.);
//                fit->SetParLimits(2, 7.1, 7.3);
                break;
            }
                
            case Function::TBW:
                break;
            case Function::RightKneeFunction:
                break;
            case Function::RightKneeFunctionApproximation:
                break;
        }
    }

//    TFitResultPtr fitResult = h->Fit(fit, "SMQR");
    TFitResultPtr fitResult = h->Fit(fit, "SQR");
    lastSeedParams_.assign(fit->GetNpar(), 0.0);
    for (int i = 0; i < fit->GetNpar(); ++i)
        lastSeedParams_[i] = fit->GetParameter(i);
    
    fitResults.push_back(fitResult);

    // Create fit log
//    std::cout << "Creating fitlog" << std::endl;
    std::string histName = h->GetName();
    std::string logFileName = "Fitlogs/" + histName + "_fitlog.txt";
    mkdir("Fitlogs", 0777);
    std::ofstream fitLog(logFileName);
    if (!fitLog.is_open()) {
        std::cerr << "Error opening " << logFileName << " for writing." << std::endl;
        return fit;
    }

    fitLog << "****************************************" << std::endl;
    fitLog << "Minimizer is Minuit2 / Migrad" << std::endl;
    fitLog << "Chi2                      = " << fitResult->Chi2() << std::endl;
    fitLog << "NDf                       = " << fitResult->Ndf() << std::endl;
    fitLog << "Chi2 / NDF                = " << fitResult->Chi2() / fitResult->Ndf() << std::endl;
    fitLog << "covariance                = " << fitResult->CovMatrixStatus() << std::endl;
    fitLog << "edm                       = " << fitResult->Edm() << std::endl;

    for (int i = 0; i < fit->GetNpar(); ++i) {
        double lowLimit, highLimit;
        fit->GetParLimits(i, lowLimit, highLimit);
        fitLog << fit->GetParName(i) << " = "
               << fit->GetParameter(i) << " +/- " << fit->GetParError(i);
        if (lowLimit != highLimit)
            fitLog << "  (limits: [" << lowLimit << ", " << highLimit << "])";
        fitLog << std::endl;
    }

    fitLog << "\n\n---------------FIT FUNCTION------------------------\n";
//    double* pars = fit->GetParameters();
//    if (pars != nullptr) {
//        for (double x = 0; x < 2.01; x += 0.01) {
//            fitLog << "f(" << x << ") = " << model(&x, pars) << std::endl;
//        }
//    } else {
//        fitLog << "ERROR: Fit parameters are null — skipping function evaluation.\n";
//    }

    fitLog << "\n\n---------------DATASET------------------------\n";
    for (double x = 0; x < 2.01; x += 0.01) {
        int bin = h->FindBin(x);
        double y = h->GetBinContent(bin);
        fitLog << "Data(" << x << ") = " << y << std::endl;
    }
    fitLog.close();

    // Draw stat box
    double chi2 = fitResult->Chi2();
    int ndf = fitResult->Ndf();
    double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;

    TPaveText* statBox = new TPaveText(0.5, 0.47, 0.9, 0.75, "NDC");
    statBox->SetBorderSize(1);
    statBox->SetFillColor(kWhite);
    statBox->AddText(Form("#chi^{2} / NDF = %.2f / %d = %.3f", chi2, ndf, chi2ndf));
    statBox->AddText(fitname.c_str());
    for (int i = 0; i < fit->GetNpar(); ++i) {
        statBox->AddText(Form("%s = %.4g #pm %.4g",
                              fit->GetParName(i),
                              fit->GetParameter(i),
                              fit->GetParError(i)));
    }
    statBox->SetTextSize(0.03);
    h->GetListOfFunctions()->Add(statBox);

//    std::cout << "Done fitting" << std::endl;
    return fit;
}

int Parameterization::GetNParams(Parameterization::Function f) {
    switch (f) {
        case Parameterization::Function::LevyTsallis:       return 4;
        case Parameterization::Function::Hagedorn:          return 5;
        case Parameterization::Function::BW:                return 5;
        case Parameterization::Function::TBW:               return 5;
        case Parameterization::Function::FokkerPlanck:      return 6;
        case Parameterization::Function::Tsallis:           return 4;
        case Parameterization::Function::PowerLaw:          return 2;
        case Parameterization::Function::ExpDecay:          return 3;
        case Parameterization::Function::RightKneeFunction: return 1;
        case Parameterization::Function::RightKneeFunctionApproximation: return 3;
        default:                                            return 0;
    }
}

void Parameterization::SetDefaultFitParameters(TF1* fit, Parameterization::Function f) {
    switch (f) {
        case Parameterization::Function::LevyTsallis:
            fit->SetParameters(0.1, 0.2, 12, 6000);
            break;
        case Parameterization::Function::Tsallis:
            fit->SetParameters(1.0, 1.1, 0.15);
            break;
        case Parameterization::Function::Hagedorn:
            fit->SetParameters(0.139, 0.2, 7, 100);
            break;
        case Parameterization::Function::TBW:
            fit->SetParameters(1.0, 0.139, 0.114, 1.086, 0.);
            break;
        case Parameterization::Function::FokkerPlanck:
            fit->SetParameters(0.139, 0.5, 1.0, 2.0, 2.0, 0.15);
            break;
        case Parameterization::Function::PowerLaw:
            fit->SetParameters(6, 100);
            break;
        case Parameterization::Function::ExpDecay:
            fit->SetParameters(1.0, 0.1, 2.0);
            break;
        default:
            break;
    }
}



//-----------------------------------------------------

TH1* Parameterization::GetDataFitRatio(TH1 *h, TF1* fit, std::string title){
    TH1* ratio = new TH1F(title.c_str(), title.c_str(),
                              h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
        ratio->GetYaxis()->SetTitle("Data/Fit");
        ratio->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        
        for (int i = 1; i <= h->GetNbinsX(); ++i) {
            double data = h->GetBinContent(i);
            double fitValue = fit->Eval(h->GetBinCenter(i));
            if (fitValue != 0 && data != 0) {
                double ratioValue = data / fitValue;
                ratio->SetBinContent(i, ratioValue);
            } else {
                ratio->SetBinContent(i, 0);
            }
        }
        
        return ratio;
}

//-----------------------------------------------------

void Parameterization::NormalizeByBinWidth(TH1 *h) {
  for (int bin = 0; bin <= h->GetNbinsX(); bin++) {
    double content{h->GetBinContent(bin)};
    double center{h->GetBinCenter(bin)};
    double width{h->GetBinWidth(bin)};
      double acceptance_factor = 1.;
      // acceptance_factor = 1.842; //determined by self-weight scheme (or comparison to PHOBOS dN/deta) <- w/ and w/o 1/pT
      
      acceptance_factor = 2 * M_PI / acceptance_factor;
//      acceptance_factor = 1 / acceptance_factor;
//      content = content / ( width * center * acceptance_factor);
      content = content / ( width * acceptance_factor);
//      content = content / ( width * center * M_PI);
//      content = content / ( width * center * 2 * M_PI);
//      content = content / ( width * center );
//      content = content / ( width * acceptance_factor );
//    h->Sumw2();
    h->SetBinContent(bin, content);
    
  }
}

//-----------------------------------------------------

void Parameterization::UnnormalizeByBinWidth(TH1 *h) {
    for (int bin = 1; bin <= h->GetNbinsX(); bin++) {
        double content{h->GetBinContent(bin)};
        double center{h->GetBinCenter(bin)};
        double width{h->GetBinWidth(bin)};
        h->SetBinContent(bin, content * ( width * center * 2 * M_PI));
    }
}

//-----------------------------------------------------

double Parameterization::GetMean(TH1 *h, TF1 *fit, bool usecut,
                                 const double &high_cut, const double &low_cut) {
    double integral{0}; bool debug = false;
    double fit_cut{0.2}; //200 MeV pT cutoff
//    double fit_cut{0.4}; //200 MeV pT cutoff
//    double fit_cut{0.15}; //150 MeV pT cutoff
//    double fit_cut{0.0};  //ignore fit altogether
//    double fit_cut{10.};
    double mean{0};
    if(debug){
        std::cout << "\nhigh cut: " << high_cut << " | " << "low cut: " << low_cut << " | " << "fit cut: " << fit_cut << "usecut: " << usecut << std::endl;
        std::cout << "pT bin center | h->GetBinContent(bin) | fit->Eval(center) | content" << std::endl;
    }

  for (int bin = 0; bin <= h->GetNbinsX(); ++bin) {
    double center{h->GetBinCenter(bin)}; // pT
    double width{h->GetBinWidth(bin)};   // dpT
    double content{fit->Eval(center)};   // f(pT)
//    double error(h->GetBinError(bin));
//    double content{h->GetBinContent(bin)};
//      content = 0; //Removing fit contributions
      if (usecut && center > fit_cut){
          content = h->GetBinContent(bin);
      }
    if (usecut && center > high_cut) {
        continue;
    }
    if (usecut && center < low_cut) {
        continue;
      }
      if (content < 0) continue;
      if(debug && bin < 20){
      std::cout << center << " | " << h->GetBinContent(bin) << " | " << fit->Eval(center) << " | " << content << std::endl;
      }
//
//    integral += width * center * content;
//    mean += center * center * (width * content); // https://arxiv.org/pdf/0809.2482 page 20 (1/pT dN/dptdeta)
      integral += width * content;
      mean += center * (width * content);  //Omar's approach (when you have dN/dptdeta)
  }
  mean = mean / integral;
    if(debug){ std::cout << "mean: " << mean << std::endl;}

  return mean;
}

//-----------------------------------------------------

double Parameterization::GetdNdEta(TH1 *h, TF1 *fit, bool usecut,
                                   const double &high_cut, const double &low_cut) {
    //Note: low and high cuts are not applied for dNdEta class.
//    double fit_cut{0.2};
//      double fit_cut{0.15};
    double fit_cut{0.0}; //ignore fit altogether
//    double fit_cut{10.};
    double integral{0};
    bool debug = false;
    if (debug){
        std::cout << "\nhigh cut: " << high_cut << " | " << "low cut: " << low_cut << " | " << "fit cut: " << fit_cut << "usecut: " << usecut << std::endl;
        std::cout << "pT bin center | h->GetBinContent(bin) | fit->Eval(center) | content" << std::endl;
    }

    for (int bin = 0; bin <= h->GetNbinsX(); ++bin) {
      double center{h->GetBinCenter(bin)}; // pT
      double width{h->GetBinWidth(bin)};   // dpT
      double content{fit->Eval(center)};   // f(pT)
        if (usecut && center > fit_cut){
            content = h->GetBinContent(bin);
        }
        if (content < 0) continue;
//      if (usecut && center > high_cut) {
//          continue;
//      }
//      if (usecut && center < low_cut) {
//          continue;
//        }
//      integral += center * content; //cross-checked with https://arxiv.org/pdf/nucl-ex/0201005 page 5
//      integral += width * center * content; //cross-checked with https://arxiv.org/pdf/0809.2482 page 20
        integral += width * content;
        if(debug && bin < 20){
        std::cout << center << " | " << h->GetBinContent(bin) << " | " << fit->Eval(center) << " | " << content << std::endl;
        }
    }

    double acceptance_factor = 1.;
//    acceptance_factor = 10.606;
//    acceptance_factor = 4.375; //default
//    acceptance_factor = 0.2855;
    acceptance_factor = 2 * M_PI / acceptance_factor;
    integral /= acceptance_factor;
    if (debug) cout << "<dN/deta>_{pT Integrated} = " << integral << endl;
    return integral;
}

//-----------------------------------------------------

double Parameterization::GetNch(TH1 *h, const int low_ref, const int high_ref) {
    double integral{0};
    double mean{0}; bool debug = false;
    int bin_low = h->GetXaxis()->FindBin(low_ref);
    int bin_high = h->GetXaxis()->FindBin(high_ref);
    if (bin_high < bin_low){return kFatal;}
    for (int bin = bin_low; bin <= bin_high; ++bin) {
      double center{h->GetBinCenter(bin)};
      double width{h->GetBinWidth(bin)};
  //    double content{fit->Eval(center)};
      double content{h->GetBinContent(bin)};
//      integral += center * width * content;
//      mean += center * center * (width * content);
        integral += width * content;
      mean += center * (width * content);
    }
    mean /= integral;
    if (debug) std::cout << "low_ref = " << low_ref << ", high ref = " << high_ref << " <dN/deta>_{StRefMultCorr}: " << mean << std::endl;
    return mean;
}

//-----------------------------------------------------

TH1* Parameterization::GetSpectraDifference(TH1 *h, TH1* ref, std::string title){
    TH1* difference = new TH1F(title.c_str(), title.c_str(),
                              h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
        difference->GetYaxis()->SetTitle("Relative Difference (Note: Ideally <<1)");
        difference->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        
        for (int i = 1; i <= h->GetNbinsX(); ++i) {
            double hcontent = h->GetBinContent(i);
            double refcontent = ref->GetBinContent(i);
            if (hcontent != 0 && refcontent != 0) {
                double diff = (hcontent - refcontent)/(refcontent); //Relative difference
                difference->SetBinContent(i, diff);
            } else {
                difference->SetBinContent(i, 0);
            }
        }
        return difference;
}

//-----------------------------------------------------

double Parameterization::GetRange(const int k, Parameterization params, bool display) {
    double range{0};
    range = abs(static_cast<double>(params.GetFineAmplitudeBinCutHigh(k))
          - static_cast<double>(params.GetFineAmplitudeBinCutLow(k)));
    if (display){
        cout << "\nGetRange: k = " << k+1 << " | ";
        double centrality{params.GetCentrality(k, params, true)};
        cout << "(" << centrality << "%) | ";
        cout << "high Ref Mult bin: " << params.GetFineAmplitudeBinCutHigh(k) << " | ";
        cout << "low Ref Mult bin: " << params.GetFineAmplitudeBinCutLow(k);
        cout << " | Reference multiplicity range: " << range << endl;
    }
    return range;
}

//-----------------------------------------------------

double Parameterization::GetDensity(const int k, Parameterization params) {
    //Getter for Ncharged bin density (for centrality resolution eval)
    double density{0};
//    density = static_cast<double>(params.GetFineAmplitudeBinCutHigh(k))
//          - static_cast<double>(params.GetFineAmplitudeBinCutLow(k));
//    cout << "Reference multiplicity bin density: " << density << endl;
    return density;
}

//-----------------------------------------------------

double Parameterization::GetCentrality(const int k, Parameterization params, bool display) {
    //Getter for Ncharged bin density (for centrality resolution eval)
    double centrality{0};
    double high = params.GetFineAmplitudeCutHigh(k)*100;
    double low = params.GetFineAmplitudeCutLow(k)*100;
    centrality = low + abs(( high - low ) / 2 );
    if (display) {
        cout << "centrality: " << low << " - " << high << "%";
    }
    return centrality;
}

//-----------------------------------------------------

double Parameterization::GetErrorMeanpT(const TH1 *h, TF1 *fit, bool usecut,
                                        const double &cut, bool debug) {
    // ^Error progation-> https://en.wikipedia.org/wiki/Propagation_of_uncertainty
    //                    https://arxiv.org/pdf/1609.04150#page24
    double uncertainty{0}; double fit_cut{0.0};
    
    //fit error contribution
    TFitResultPtr fitResult = fitResults[index];
    if (!fitResult) {
        std::cerr << "\n⚠️ [GetErrorMeanpT] No fit result found (fitResults[index] is nullptr)" << std::endl;
    }
    const double* params = fitResult->GetParams(); // Parameters of the fit
    const TMatrixDSym& covMatrix = fitResult->GetCovarianceMatrix(); // Covariance matrix

    if (debug){
        std::cout << "\n*🔧 Processing error on the ⟨pT⟩ in [0," << cut << "] (GeV/c)*" << std::endl;
        std::cout << "Fit Parameters:\n";
        for (int i = 0; i < fitResult->NPar(); ++i) {
            std::cout << "  Param[" << i << "]: " << params[i] << std::endl;
        }
        std::cout << "Covariance Matrix:\n";
        for (int i = 0; i < fitResult->NPar(); ++i) {
            for (int j = 0; j < fitResult->NPar(); ++j) {
                std::cout << std::setw(10) << covMatrix(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }
    
    // Numerical integrals and their errors
    double fit_num{0}, fit_num_err{0};
    double fit_den{0}, fit_den_err{0};
    
    if (debug) {
        std::cout << "[DEBUG-FIT] error = fit->IntegralError(pt, pt + width, params, covMatrix.GetMatrixArray()) = " << std::endl;
    }
    for (int bin = 0; bin <= h->GetNbinsX(); ++bin) {
        double pt = h->GetBinCenter(bin);
        double width = h->GetBinWidth(bin);
        if (pt > fit_cut) continue;

        double value = fit->Eval(pt);  // Fit value at pT
        double error = fit->IntegralError(pt, pt + width, params, covMatrix.GetMatrixArray());
        if (debug && bin < 3){
                 std::cout << "[DEBUG-FIT] [" << bin << "] error = " << error << std::endl;
        }
        fit_num += width * pt * pt * value; // numerator: pT^2 * del(pT) * sum Y(pT)
        fit_den += width * pt * value;      // denominator: pT * del(pT) * sum Y(pT)
        
        fit_den_err += pow(width * pt * error, 2.0);       // pT * del(pT) * error(pT)
        fit_num_err += pow(width * pt * pt * error, 2.0);  // pT^2 * del(pT) * error(pT)
    }
    double fit_err;
    if (debug) {
        std::cout << "[DEBUG-FIT] fit_num = " << fit_num << ", fit_den = " << fit_den << std::endl;
        std::cout << "[DEBUG-FIT] fit_num_err = " << sqrt(fit_num_err) << ", fit_den_err = " << sqrt(fit_den_err) << std::endl;
    }
    if (fit_num <= 1e-20 || fit_den <= 1e-20) {
        if (debug){
            std::cout << "⚠️  [GetErrorMeanpT] Invalid numerator/denominator for fit error calc: fit_num = " << fit_num << ", fit_den = " << fit_den << std::endl;
        }
        fit_err = 0; // bad fit error calculation so ignores it
    } else {
            fit_err = ( fit_num / fit_den ) * sqrt(
            pow(fit_num_err / fit_num, 2.0 ) + pow(fit_den_err / fit_den, 2.0 )
        );
    }
    
    //dataset error contribution
    double data_num_err{0}; double data_den_err{0};
    double data_num{0}; double data_den{0};

    if (debug) {
        std::cout << "[DEBUG-DATA] error = h->GetBinError(bin)  = " << std::endl;
    }
    
    int ctr{0};
  for (int bin = 0; bin <= h->GetNbinsX(); ++bin) {
    double center{h->GetBinCenter(bin)};
    double width{h->GetBinWidth(bin)};
    double error{h->GetBinError(bin)};
    double content{h->GetBinContent(bin)};
      if (usecut && center < fit_cut) { 
//          if (debug && ctr < 1) { std::cout << "[DEBUG-DATA] Excluding at lower cut " << fit_cut << "... (bin error = " << error << ")" << std::endl; }
          continue;
      }
      // high_pT cutoff (for the data error contribution):
      if (usecut && center > cut) { 
          if (debug && ctr < 1) { std::cout << "[DEBUG-DATA] Excluding at upper cut " << cut << "... (bin error = " << error << ")" << std::endl; }
          continue;
      }
      ctr++;
      if (debug && ctr < 3){
               std::cout << "[DEBUG-DATA] [" << bin << "] error = " << error << std::endl;
      }
//      data_num += ( pow(center, 2.0) * width * content ); // numerator: pT^2 * del(pT) * sum Y(pT)
//      data_den += center * width * content;               // denominator: pT * del(pT) * sum Y(pT)
//      
//      data_num_err += pow( ( pow(center, 2.0) * width * error ), 2.0 ); // pT^2 * del(pT) * error(pT)
//      data_den_err += pow( ( center * width * error ), 2.0 );           // pT * del(pT) * error(pT)
      
      //Omar's approach (when you have dN/dptdeta)
      data_num += ( center * width * content ); // numerator: pT^2 * del(pT) * sum Y(pT)
      data_den += width * content;               // denominator: pT * del(pT) * sum Y(pT)
      
      data_num_err += pow( ( center * width * error ), 2.0 ); // pT^2 * del(pT) * error(pT)
      data_den_err += pow( ( width * error ), 2.0 );           // pT * del(pT) * error(pT)
      
  }
    //Final calculation
    
    double data_err = ( data_num / data_den ) * sqrt ( pow( data_num_err / data_num , 2.0) + pow ( data_den_err / data_den, 2.0 ));
//    if (std::isnan(data_err)) { data_err = 0; }
    fit_err = 0.;
        
    uncertainty = pow(fit_err, 2.0) + pow(data_err,2.0);
    
    if (debug){
        cout << " fit error | data error | uncertainty" << endl;
        cout << fit_err << " | " << data_err << " | " << sqrt(uncertainty) << endl;
        cout << "Error (<pT>): " << sqrt(uncertainty) << endl;
        std::cout << "*✅ Finished processing error on the ⟨pT⟩*" << endl;
    }
    
    
  return sqrt(uncertainty);
}

//-----------------------------------------------------

double Parameterization::GetErrordNdEta(const TH1 *h, TF1 *fit, bool usecut,
                                        const double &cut, bool debug) {
    double uncertainty{0}; double fit_cut{0.2}; debug = false;
  
    
    //fit error contribution
    TFitResultPtr fitResult = fitResults[index];
    if (!fitResult) {
        std::cerr << "⚠️ [GetErrordNdEta] No fit result found (fitResults[index] is nullptr)" << std::endl;
    }
    const double* params = fitResult->GetParams(); // Parameters of the fit
    const TMatrixDSym& covMatrix = fitResult->GetCovarianceMatrix(); // Covariance matrix

    if (debug){
        std::cout << "*Processing error on the ⟨dN/dη⟩*" << endl;
        std::cout << "Fit Parameters:\n";
        for (int i = 0; i < fitResult->NPar(); ++i) {
            std::cout << "  Param[" << i << "]: " << params[i] << std::endl;
        }
        std::cout << "Covariance Matrix:\n";
        for (int i = 0; i < fitResult->NPar(); ++i) {
            for (int j = 0; j < fitResult->NPar(); ++j) {
                std::cout << std::setw(10) << covMatrix(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }
    
    double fit_err{0}; //Fit error contribution
    
    if (debug) {
        std::cout << "[DEBUG] error = fit->IntegralError(pt, pt + width, params, covMatrix.GetMatrixArray()) = " << std::endl;
    }
    for (int bin = 0; bin <= h->GetNbinsX(); ++bin) {
        double pt = h->GetBinCenter(bin);
        double width = h->GetBinWidth(bin);
        if (pt > fit_cut) continue;

        double error = fit->IntegralError(pt, pt + width, params, covMatrix.GetMatrixArray());
        if (debug && bin < 10){
                 std::cout << "[" << bin << "] error = " << error << std::endl;
        }
        
        fit_err += pow(width * pt * error, 2.0);       // (pT * del(pT) * error(pT))^2
        //Omar's approach (when you have dN/dptdeta)
        fit_err += pow(width * error, 2.0);       // (pT * del(pT) * error(pT))^2
    }
    
    //Data error contribution
  for (int bin = 0; bin <= h->GetNbinsX(); ++bin) {
    double center{h->GetBinCenter(bin)};
    double width{h->GetBinWidth(bin)};
    double error{h->GetBinError(bin)};
      if (usecut && center < cut) continue;
      
//    uncertainty += pow(center,2.0) * pow(width, 2.0) * pow(error, 2.0); // (pT * del(pT) * error(pT))^2
      
    uncertainty += pow(width, 2.0) * pow(error, 2.0); // (del(pT) * error(pT))^2
  }
    fit_err = 0.; //No statistical uncertainty should come from fit, only systematic (which is determined during cs2 variations)
    uncertainty = uncertainty + pow(fit_err, 2.0);
    
    if (debug) {
        cout << "Error (dN/dEta): " << sqrt(uncertainty) << endl;
        std::cout << "*⚙️ Finished processing error on the ⟨dN/dη⟩*" << endl;
    }
        
        
  return sqrt(uncertainty);
}

//-----------------------------------------------------

void Parameterization::SetErrorYield(TH1D*& h, int centrality_index){
    bool debug = false;
    std::string directory =
    "/Users/calebbroodo/Desktop/sos_analysis/src/ana_code/yield_errors/yield_error_";
    std::string filename = directory + std::to_string(centrality_index) + ".txt";
    
    std::ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
            std::this_thread::sleep_for(std::chrono::milliseconds(10000));
        }
    
    std::string line;
        // Skip the first two lines
        std::getline(infile, line); // Skip the centrality header
        std::getline(infile, line); // Skip the column header
    
    TH1D* h_out = (TH1D*)h->Clone();
    h_out->SetDirectory(0);
    
    
    if (debug){
        std::cout << "\n Initial bin errors in pT spectra: " << std::endl;
        for (int bin = 0; bin < 40; bin++) {
            std::cout << " [" << bin << "] error: " << h->GetBinError(bin) << endl;
        }
    }
    
    while (std::getline(infile, line)) {
            std::stringstream ss(line);
            double binCenter, stdDev, mean;
            char separator;

            if (ss >> binCenter >> separator >> stdDev >> separator >> mean) {
                // Find the bin corresponding to the bin center
                int bin = h_out->FindBin(binCenter);
                h_out->SetBinError(bin, stdDev);
//                h_out->SetBinContent(bin, mean); // only set the error on the mean
            } else {
                std::cerr << "Error: Malformed line in file: " << line << std::endl;
            }
        }
    
    infile.close();
    delete h;
    h = h_out;
    

    if (debug){
        std::cout << "\n Setting error in pT spectra: " << std::endl;
        for (int bin = 1; bin < 10; bin++) {
            std::cout << " [" << bin << "] error: " << h->GetBinError(bin) << endl;
        }
    }
    
}

//-----------------------------------------------------

double Parameterization::GetRelativeError(int bin) const {
    int nBins{this->GetnFinePercentileBinning()};
    double valreturn{-1};
    if (bin < 0 || bin > nBins) {
        cout << "Error: Bin out of bounds." << '\n';
        valreturn = -1;
    }
    
    if (bin != nBins - 1) {
        if (this->GetSystem() == "AuAu_200") {
            if (this->GetMultEstimator() == "TPC"){
                valreturn = RelErrorCutsTPC[bin];
            }
            else {
                cout << "MultEstimator not defined!" << '\n';
            }
        }
        if (this->GetSystem() == "AuAu_3") {
            if (this->GetMultEstimator() == "TPC") {
                valreturn = RelErrorCutsTPC[bin];
            }
            else {
                cout << "MultEstimator not defined!" << '\n';
            }
        }
    }
    return valreturn;
}

//-----------------------------------------------------

TH1* Parameterization::GetSpectraRatio(TH1 *h, TH1 *href, TF1 *fit, bool usecut,
                                       const double &high_cut, const double &low_cut) {
    if (!h || !href) {
        std::cerr << "Error: One or both histograms are null." << std::endl;
        return nullptr;
    }

    // custom binning for better illustration of high pT effects
    std::vector<double> custom_bins;
    double bin_high_edge = h->GetXaxis()->GetXmax();

    int rebin_factor = 10;
    for (int bin = h->GetXaxis()->FindBin(0.1); bin <= h->GetXaxis()->FindBin(1); bin += rebin_factor) {
        custom_bins.push_back(h->GetBinLowEdge(bin));
    }
    rebin_factor = 20;
    for (int bin = h->GetXaxis()->FindBin(1); bin <= h->GetXaxis()->FindBin(2); bin += rebin_factor) {
        custom_bins.push_back(h->GetBinLowEdge(bin));
    }
    rebin_factor = 200;
    for (int bin = h->GetXaxis()->FindBin(2); bin <= h->GetXaxis()->FindBin(6); bin += rebin_factor) {
        custom_bins.push_back(h->GetBinLowEdge(bin));
    }
    
    rebin_factor = 400;
    for (int bin = h->GetXaxis()->FindBin(6); bin <= h->GetNbinsX(); bin += rebin_factor) {
        custom_bins.push_back(h->GetBinLowEdge(bin));
    }
    custom_bins.push_back(bin_high_edge); // Include the upper edge of the last bin

    // Create a new histogram with the custom binning
    int n_bins = custom_bins.size() - 1;
    TH1* h_rebinned = new TH1D("h_rebinned", "Rebinned Histogram", n_bins, custom_bins.data());
    TH1* href_rebinned = new TH1D("href_rebinned", "Rebinned Reference Histogram", n_bins, custom_bins.data());
    
    cout << "\n href bins: " << href->GetNbinsX();
    cout << "\n href bins: " << h->GetNbinsX();
    cout << "\n Content: ";
    for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
        double content = h->GetBinContent(bin);
        double error = h->GetBinError(bin);
        double x = h->GetBinCenter(bin);
        h_rebinned->Fill(x, content);
        h_rebinned->SetBinError(h_rebinned->FindBin(x), error);
        if (bin % 100 == 0) cout << "[ " << content << " , ";
        
        content = href->GetBinContent(bin);
        error = href->GetBinError(bin);
        href_rebinned->Fill(x, content);
        href_rebinned->SetBinError(href_rebinned->FindBin(x), error);
        if (bin % 100 == 0) cout << content << " ], ";
    }
    cout << endl;
    // Calculate ratio
    TH1* ratio = (TH1*)h_rebinned->Clone("ratio");
    ratio->Divide(href_rebinned);

    // Apply cuts if needed
    if (usecut) {
        int bin_low = ratio->GetXaxis()->FindBin(low_cut);
        int bin_high = ratio->GetXaxis()->FindBin(high_cut);
        for (int bin = 1; bin <= ratio->GetNbinsX(); ++bin) {
            if (bin < bin_low || bin > bin_high) {
                ratio->SetBinContent(bin, 0); // Set bins outside the range to 0
                ratio->SetBinError(bin, 0);   // Set bin errors to 0
            }
        }
    }
    
    ratio->GetXaxis()->SetRangeUser(0.2, 10.0);
    ratio->SetXTitle("p_{T} (GeV/c)");

    return ratio;
}

//-----------------------------------------------------

double Parameterization::P_n(const double* x, const double* par) {
    bool debug{false};
    double n = x[0];            // n (x-axis value)
    double n_knee = par[0];     // Fit parameter: n_knee
    double a1 = par[1];         // Fit parameter: a1
    double a2 = par[2];         // Fit parameter: a2
    double a3 = par[3];         // Fit parameter: a3
    double sigma = par[4];      // Fit parameter: sigma

    double integral = 0.0;
    int n_steps = 1000;   //n number of integration steps
    double step = 1.0 / n_steps;

    for (int i = 0; i <= n_steps; ++i) {
        double c_b = i * step; // integrating variable
        double exponent = n_knee * TMath::Exp(-a1 * c_b - a2 * c_b * c_b - a3 * c_b * c_b * c_b);
        double arg = (n - exponent) / sigma;
        integral += TMath::Exp(-0.5 * arg * arg) / (sigma * TMath::Sqrt(2 * TMath::Pi()));
        if (i == n_steps && debug){
            
            std::cout << "Step " << i
            << ": n = " << n
            << ", n_knee = " << n_knee
            << ", a1 = " << a1
            << ", a2 = " << a2
            << ", a3 = " << a3
            << ", sigma = " << sigma
            << ", c_b = " << c_b
            << ", exponent = " << exponent
            << ", arg = " << arg
            << ", term = " << TMath::Exp(-0.5 * arg * arg) / (sigma * TMath::Sqrt(2 * TMath::Pi()))
            << ", integral*step = " << integral * step << std::endl;
        }
    }
    return integral * step; // multiplying step size to approximate the integral
}

//-----------------------------------------------------

TH1D* Parameterization::GetKneeParams(TH1D *h, int ref, bool norm, bool fit, bool ratio) {
    //Make boolean in arguments for returning scaled x-axis or not scaled x-axis
    Parameterization params;
    // https://arxiv.org/pdf/1708.00081 using scenario (B) sigma(c_b) = sigma(0) = sigma
    TH1D* h_uncorr = (TH1D*)h->Clone("h_uncorr");
    h_uncorr->Scale(1./h_uncorr->Integral());
    int nBins = h->GetNbinsX();
    h_uncorr->GetYaxis()->SetTitle("Event fraction (a.u.)");
    
    //Modulate
    double n_min = 500; // minimum Ncharged cutoff ( as recommended by https://arxiv.org/pdf/1708.00081 )
//    double n_max = 800;
    double n_max = 734.38;
    
    TRandom3 rand;
    //Unccorelated errors using Guassian input:
        for (int i = 1; i <= nBins; ++i) {
            double content = h_uncorr->GetBinContent(i);
            double error = h_uncorr->GetBinError(i);
            // Generate new bin content with Gaussian fluctuations
            if (error < 0 || std::isnan(error)) {
                    std::cerr << "Invalid error in bin " << i << std::endl;
                }
            double new_content = rand.Gaus(content, error);
            h_uncorr->SetBinContent(i, new_content);
            h_uncorr->SetBinError(i, error); // Preserve errors
        }

    if (norm){
        cout << "*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*" << endl;
        cout << "Getting knee parameters (normalized fit)..." << endl;
        TH1D* h_fitvals = (TH1D*)h_uncorr->Clone("h_fitvals");
        TF1* fitFunc = new TF1("fitFunc", Parameterization::P_n, n_min,
                               n_max, 5); // 5 parameters
//        TF1* fitFunc = new TF1("fitFunc", Parameterization::P_n, static_cast<double>(n_min / ref),
//                               static_cast<double>(n_max / ref), 5); // 5 parameters
        
        fitFunc->SetParName(0, "N_{ch,knee}");
        fitFunc->SetParName(1, "a_{1}");
        fitFunc->SetParName(2, "a_{2}");
        fitFunc->SetParName(3, "a_{3}");
        fitFunc->SetParName(4, "#sigma_{0}");
        
        h_uncorr->GetXaxis()->SetTitle(Form("#LTdN_{ch}/d#eta#GT / #LTdN_{ch}/d#eta#GT_{0-5%%}"));
        
        // Initial guesses: n_knee, a1, a2, a3, sigma_0
            fitFunc->SetParameters(560, 4, 0.95, -0.3, 35);
        
        h_uncorr->Fit(fitFunc, "S");
        
        for (int bin = 1; bin <= nBins; ++bin) {
            double x = h_fitvals->GetBinCenter(bin); // Get the bin center
            double y = fitFunc->Eval(x);              // Evaluate the fit function at x
            if (x > n_min) {h_fitvals->SetBinContent(bin, y); h_fitvals->SetBinError(bin, 0);}
        }
        
        h_uncorr->GetXaxis()->SetLimits(h_uncorr->GetXaxis()->GetXmin() / ref,
                                        h_uncorr->GetXaxis()->GetXmax() / ref);
        h_fitvals->GetXaxis()->SetLimits(h_fitvals->GetXaxis()->GetXmin() / ref, h_fitvals->GetXaxis()->GetXmax() / ref);
        
        double limit = 1.415; // Set the x-axis limit
        int nBins = h_fitvals->GetNbinsX();

        for (int i = 1; i <= nBins; ++i) {
            double binCenter = h_fitvals->GetXaxis()->GetBinCenter(i);
            if (binCenter > limit) {
                h_fitvals->SetBinContent(i, 0); // Nullify content
                h_fitvals->SetBinError(i, 0);   // Nullify error
            }
        }

        std::cout << "h_fitvals X-axis limits: ["
        << h_fitvals->GetXaxis()->GetXmin() << ", "
        << h_fitvals->GetXaxis()->GetXmax() << "]" << std::endl;
        
        std::cout << "Histogram X-axis limits: ["
        << h_uncorr->GetXaxis()->GetXmin() << ", "
        << h_uncorr->GetXaxis()->GetXmax() << "]" << std::endl;
        
        double chi2 = fitFunc->GetChisquare();
        double ndf = fitFunc->GetNDF();
        std::cout << "Chi^2 / NDF: " << chi2 << " / " << ndf << " = "<<  chi2 / ndf << std::endl;
        cout << "Knee parameter acquisition complete..." << endl;
        cout << "*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*" << endl;
        if (ratio) {
            TH1D* ratio = new TH1D("Data/Fit (norm)", "Data/Fit (norm)",
                                      h_uncorr->GetNbinsX(), h_uncorr->GetXaxis()->GetXmin(), h_uncorr->GetXaxis()->GetXmax());
                ratio->GetYaxis()->SetTitle("Data/Fit");
                ratio->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT");
                ratio->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT / #LTdN_{ch}/d#eta#GT_{0-5%}");
                
                for (int i = 1; i <= h_uncorr->GetNbinsX(); ++i) {
                    double data = h_uncorr->GetBinContent(i);
                    double fitValue = h_fitvals->GetBinContent(i);
                    if (fitValue != 0 && data != 0) {
                        double ratioValue = data / fitValue;
                        ratio->SetBinContent(i, ratioValue);
                    } else {
                        ratio->SetBinContent(i, 0);
                    }
                }
            return ratio;
        }
        if (fit) return h_fitvals;
    }
    else{
        cout << "*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*" << endl;
        cout << "Getting knee parameters (regular fit)..." << endl;
        TF1* fitFunc = new TF1("fitFunc", Parameterization::P_n, n_min,
                               n_max, 5); // 5 parameters
        
        // Set parameter names and initial values
        fitFunc->SetParName(0, "N_{ch,knee}");
        fitFunc->SetParName(1, "a_{1}");
        fitFunc->SetParName(2, "a_{2}");
        fitFunc->SetParName(3, "a_{3}");
        fitFunc->SetParName(4, "#sigma_{0}");
        
        h_uncorr->GetXaxis()->SetTitle(Form("#LTdN_{ch}/d#eta#GT"));
        
        
        fitFunc->SetParameters(560, 4, 0.95, -0.3, 35);
        
        h_uncorr->Fit(fitFunc, "S");
        
        double chi2 = fitFunc->GetChisquare();
        double ndf = fitFunc->GetNDF();
        std::cout << "Chi^2 / NDF: " << chi2 << " / " << ndf << " = "<<  chi2 / ndf << std::endl;
        cout << "Knee parameter acquisition complete..." << endl;
        cout << "*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*---*" << endl;
        if (ratio) return static_cast<TH1D*>(params.GetDataFitRatio(h_uncorr, fitFunc, "D/F Knee Fitting" ));
    }
    return h_uncorr;
}

//-----------------------------------------------------

double Parameterization::f(double n){
    double sigma = 0.0592; //par2
//    double sigma = 0.025;
    
//    double knee = 1.06131;
    double knee = 1.1168; //par1

    double norm = sigma * TMath::Sqrt(2 / TMath::Pi());
    double exponent = TMath::Exp(-TMath::Power(n - knee, 2) / (2 * TMath::Power(sigma, 2)));
    double erfc = TMath::Erfc((n - knee) / (TMath::Sqrt(2) * sigma));

    return n - (norm * exponent / erfc);
}

TH1D* Parameterization::GetBlackCurve(){
    //https://arxiv.org/pdf/1909.11609
    double n_min = 0.8;
    double n_max = 1.3;
    int n_bins = 1000;

    TH1D* h_black = new TH1D("h_black", "f(N_{ch}/N_{ch}_{0-5%})", n_bins, n_min, n_max);
    h_black->GetYaxis()->SetTitle("f(N_{ch}/N_{ch}_{0-5%})");
    h_black->GetXaxis()->SetTitle("N_{ch}/N_{ch}_{0-5%}");
    h_black->SetLineColor(kBlack);
    h_black->SetLineStyle(1);
    h_black->SetLineWidth(5);

    for (int i = 1; i <= n_bins; ++i) {
        double n = h_black->GetXaxis()->GetBinCenter(i);
        double value = Parameterization::f(n);
        h_black->SetBinContent(i, value);
    }
    return h_black;
}

TH1D* Parameterization::GetRedCurve(){
    //https://arxiv.org/pdf/1909.11609 eq.
    double n_min = 0.5;
    double n_max = 1.3;
    int n_bins = 1000;
    double cs2 = 0.08; //Setting cs2 to a const value (qualitative assesment)
//    double cs2 = 0.33;
    
    TH1D* h_red = new TH1D(Form("Analytic Model"), Form("Analytic Model"), n_bins, n_min, n_max);
//    TH1D* h_red = new TH1D(Form("Analytic Model: c_{s}^{2} = %0.2f", cs2), Form("Analytic Model: c_{s}^{2} = %0.2f", cs2), n_bins, n_min, n_max);
    h_red->GetYaxis()->SetTitle("#LTp_{T}#GT/#LTp_{T}#GT_{0-5%}");
    h_red->GetXaxis()->SetTitle("N_{ch}/N_{ch}_{0-5%}");
    h_red->SetLineColor(kRed);
    h_red->SetLineStyle(1);
    h_red->SetLineWidth(5);
    
    for (int i = 1; i <= n_bins; ++i) {
        double n = h_red->GetXaxis()->GetBinCenter(i);
        double value = pow( ( n / Parameterization::f(n) ), cs2);
        h_red->SetBinContent(i, value);
    }
    
    return h_red;
}

TH1D* Parameterization::GetBlueCurve(){
    //https://arxiv.org/pdf/1909.11609 eq. (10), Transverse Radius Ratio R/R_0
    double n_min = 0.8;
    double n_max = 1.3;
    double knee = 1.1168;
    int n_bins = 1000;
    
    TH1D* h_blue = new TH1D("h_blue", "R/R_{0-5%}", n_bins, n_min, n_max);
    h_blue->GetYaxis()->SetTitle("R/R_{0-5%}");
    h_blue->GetXaxis()->SetTitle("N_{ch}/N_{ch}_{0-5%}");
    h_blue->SetLineColor(kBlue);
    h_blue->SetLineStyle(1);
    h_blue->SetLineWidth(5);
    
    
    cout << "========-========-========-========-========-========-========-========-" << endl;
    cout << "Computing blue curve values... " << endl;
    cout << "n \t| value | f(n) | knee | f(n)/knee " << endl;
    for (int i = 1; i <= n_bins; ++i) {
        double n = h_blue->GetXaxis()->GetBinCenter(i);
        double value = pow( ( Parameterization::f(n) / knee ), (1.0 / 3.0) );
        if (i % 1000 == 0) cout << n << " | " << value << " | " << Parameterization::f(n) << " | "
            << knee << " | " << ( Parameterization::f(n) / knee ) << endl;
        h_blue->SetBinContent(i, value);
    }
    cout << "========-========-========-========-========-========-========-========-" << endl;
    return h_blue;
}

TH1D* Parameterization::GetOrangeCurve(){
    //https://arxiv.org/pdf/1909.11609 eq. (11), entropy density ratio s/s_0
    double n_min = 0.8;
    double n_max = 1.3;
    int n_bins = 1000;
    
    TH1D* h_or = new TH1D("h_or", "s/s_{0-5%}", n_bins, n_min, n_max);
    h_or->GetYaxis()->SetTitle("s/s_{0-5%}");
    h_or->GetXaxis()->SetTitle("N_{ch}/N_{ch}_{0-5%}");
    h_or->SetLineColor(kOrange);
    h_or->SetLineStyle(1);
    h_or->SetLineWidth(5);
    
    for (int i = 1; i <= n_bins; ++i) {
        double n = h_or->GetXaxis()->GetBinCenter(i);
        double value = ( n / Parameterization::f(n) );
        h_or->SetBinContent(i, value);
    }
    
    return h_or;
}

//-----------------------------------------------------

TGraphErrors* Parameterization::GetSoundFit(TGraphErrors* g, double low_cut, double high_cut, TList* outputList, double pTlow, double pThigh){
    bool debug = false;
    TGraphErrors* h_fit = (TGraphErrors*)g->Clone();
    
    TString fitName = TString::Format("PowerLaw_fit");
    TF1 *PowerLaw_fit = new TF1(fitName,Parameterization::PowerLaw,low_cut, high_cut,2); // 2 parameters (cs2 = par0, normalization = par1)
    
    PowerLaw_fit->SetParameter(0, 0.2);     // cs2
    PowerLaw_fit->SetParameter(1, 1.1);    // normalization (in accordance with knee)
    PowerLaw_fit->SetParLimits(0, 0.0, 0.33);
    PowerLaw_fit->SetParLimits(1, 0.8, 1.5);
    
    
    double linewidth = 5.0;
    PowerLaw_fit->SetLineWidth(linewidth);
    
    
    PowerLaw_fit->SetNpx(1000);
    PowerLaw_fit->SetLineColor(kBlack);
    if (debug) std::cout << "Placeholder" << std::endl;
    TFitResultPtr fitResult = nullptr;
    try {
        fitResult = h_fit->Fit(PowerLaw_fit, "SMERQ");
    } catch (...) {
        std::cerr << "Fit threw an exception!" << std::endl;
        return h_fit;  // Only valid if your function returns a pointer like TH1D*
    }
    if (debug) std::cout << "Placeholder" << std::endl;
       
    
        double cs2 = PowerLaw_fit->GetParameter(0);
        double cs2_err = PowerLaw_fit->GetParError(0);
    if (debug) std::cout << "Placeholder" << std::endl;
    
//    if (debug) {
//        for (int i = 0; i < g->GetN(); ++i) {
//            std::cout << "Y error = " << g->GetErrorY(i) << std::endl;
//        }
//        std::cout << "cs2 = " << cs2 << ", knee = " << knee << ", sigma = " << sigma << std::endl;
//    }
//    for (double n = 1.0; n < 1.25; n += 0.02) {
//        double norm = sigma * TMath::Sqrt(2 / TMath::Pi());
//        double exponent = TMath::Exp(-TMath::Power(n - knee, 2) / (2 * TMath::Power(sigma, 2)));
//        double erfc = TMath::Erfc((n - knee) / (TMath::Sqrt(2) * sigma));
//        double f = n - (norm * exponent / erfc);
//        double value = pow( ( n / f ), cs2);
//        double value = exp(cs2 * log(n / f));
//        if (debug) {
//            std::cout << "n=" << n << ", f=" << f << ", value=" << value << ", norm = " << norm
//            << ", exponent = " << exponent << ", erfc = " << erfc << std::endl;
//        }
//    }
    
    if (debug) std::cout << "Placeholder" << std::endl;
    
    double fitchi2ndf = 0.0;
    if (debug) std::cout << "computing fit result" << std::endl;
    if (fitResult.Get()) {
        if (debug) std::cout << "fit result" << std::endl;
        if (fitResult->IsValid()){
            if (debug) std::cout << "fit result" << std::endl;
            if (fitResult->Ndf() > 0) fitchi2ndf = fitResult->Chi2() / fitResult->Ndf();
        }
    }
    else return h_fit;
    
//        TString cs2Info = TString::Format("c_{s}^{2} = %.4f #pm %.6f (stat) (#chi^{2} / NDF =  %.2f)", cs2, cs2_err, fitchi2ndf);
    
    TString cs2Info = TString::Format("c_{s}^{2} = %.3f #pm %.4f (stat)", cs2, cs2_err);
    
    if (debug) std::cout << "Placeholder" << std::endl;
    
    
        TLegend* legend = new TLegend(0.2, 0.7, 0.3, 0.9);
        legend->SetTextFont(42);
        legend->SetTextSize(0.035);

    // cosmetics
//        legend->SetBorderSize(0);
//        legend->SetFillStyle(0);   // transparent
//        legend->SetMargin(0.0);    // since you're not showing a marker/line
//        legend->SetEntrySeparation(0.20); // vertical spacing between lines
    
        legend->AddEntry(h_fit, "STAR Au+Au #sqrt{s_{NN}} = 200 GeV","");
        legend->AddEntry(h_fit, Form("%.1f #leq p_{T} #leq %.1f GeV/c, |#eta| #leq 0.5", pTlow, pThigh), "");
        
    if (debug) std::cout << "Placeholder" << std::endl;
        legend->AddEntry(PowerLaw_fit, Form("Fit to extract c_{s}^{2} [%.3f,%.3f]", low_cut, high_cut), "L");
        legend->AddEntry(h_fit, cs2Info.Data(), "");
        legend->SetTextSize(0.03);
        legend->SetBorderSize(0); 
        legend->SetFillStyle(0);
        
        h_fit->GetListOfFunctions()->Add(legend);
    
        outputList->Add(h_fit);
    
    struct stat info;
    if (stat("SOS_Fitlogs", &info) != 0) {
        mkdir("SOS_Fitlogs", 0777);
    }
    if (debug) std::cout << "Placeholder" << std::endl;
    TString logFileName = Form("SOS_Fitlogs/cs2fit_pt_%.1f-%.1f_dndeta_%.2f-%.2f.txt", pTlow, pThigh, low_cut, high_cut);
    std::ofstream fitLog(logFileName.Data());
    if (!fitLog.is_open()) {
        std::cerr << "Error opening " << logFileName << " for writing." << std::endl;
        return h_fit;
    }
    if (debug) std::cout << "Placeholder" << std::endl;
    fitLog << "****************************************" << std::endl;
        fitLog << "Minimizer is Minuit2 / Migrad" << std::endl;
        fitLog << "Chi2                      = " << fitResult->Chi2() << std::endl;
        fitLog << "NDf                       = " << fitResult->Ndf() << std::endl;
        if (fitResult->Ndf() != 0)
            fitLog << "Chi2 / NDF                = " << fitResult->Chi2() / fitResult->Ndf() << std::endl;
        else
            fitLog << "Chi2 / NDF                = N/A (NDF = 0)" << std::endl;

        // Log Parameters
        for (int i = 0; i < PowerLaw_fit->GetNpar(); ++i) {
            double lowLimit, highLimit;
            PowerLaw_fit->GetParLimits(i, lowLimit, highLimit);

            fitLog << "p" << i << "                        = "
                   << std::setprecision(10) << PowerLaw_fit->GetParameter(i)
                   << "   +/-   " << std::setprecision(10) << PowerLaw_fit->GetParError(i);

            if (lowLimit != highLimit) {
                fitLog << "  (limits: [" << lowLimit << ", " << highLimit << "])";
            }
            fitLog << std::endl;
        }
    if (debug) std::cout << "Placeholder" << std::endl;
    fitLog.close();
    
    return h_fit;
}

//-----------------------------------------------------

int Parameterization::GetFineAmplitudeBinCutLow(int bin) const {
    int nBins{this->GetnFinePercentileBinning()};
    int valreturn{-1};
    
    if (bin < 0 || bin > nBins) {
        cout << "Error: Bin out of bounds." << '\n';
        valreturn = -1;
    }
    
    if (bin != nBins - 1) {
        if (this->GetSystem() == "AuAu_200") {
            if (this->GetMultEstimator() == "TPC") {
                valreturn = FineBinCutsTPC[bin + 1] + 1;
            } else if (this->GetMultEstimator() == "TracksEtaGapTPC") {
                valreturn = FineBinCutsEtaGapTPC[bin + 1] + 1;
            } else if (this->GetMultEstimator() == "Tracklets10") {
                valreturn = FineBinCutsTracklets10[bin + 1] + 1;
            } else if (this->GetMultEstimator() == "Tracklets14") {
                valreturn = FineBinCutsTracklets14[bin + 1] + 1;
            } else if (this->GetMultEstimator() == "TrackletsEtaGap") {
                valreturn = FineBinCutsTrackletsEtaGap[bin + 1] + 1;
            } else {
                cout << "fMultEstimator not defined!" << '\n';
            }
        }
        if (this->GetSystem() == "AuAu_3") {
            if (this->GetMultEstimator() == "TPC") {
                valreturn = FineBinCutsTPC[bin + 1] + 1;
            } else if (this->GetMultEstimator() == "TracksEtaGapTPC") {
                valreturn = FineBinCutsEtaGapTPC[bin + 1] + 1;
            } else if (this->GetMultEstimator() == "Tracklets10") {
                valreturn = FineBinCutsTracklets10[bin + 1] + 1;
            } else if (this->GetMultEstimator() == "Tracklets14") {
                valreturn = FineBinCutsTracklets14[bin + 1] + 1;
            } else if (this->GetMultEstimator() == "TrackletsEtaGap") {
                valreturn = FineBinCutsTrackletsEtaGap[bin + 1] + 1;
            } else {
                cout << "fMultEstimator not defined!" << '\n';
            }
        }
        
    }
    return valreturn;
}
//-----------------------------------------------------

int Parameterization::GetFineAmplitudeBinCutHigh(int bin) const {
  int nBins{this->GetnFinePercentileBinning()};
  int valreturn{-1};

  if (bin < 0 || bin > nBins) {
    cout << "Error: Bin out of bounds." << '\n';
    valreturn = -1;
  }

  if (this->GetSystem() == "AuAu_200") {
    if (this->GetMultEstimator() == "TPC") {
      valreturn = FineBinCutsTPC[bin];
    } else if (this->GetMultEstimator() == "TracksEtaGapTPC") {
      valreturn = FineBinCutsEtaGapTPC[bin];
    } else if (this->GetMultEstimator() == "Tracklets10") {
      valreturn = FineBinCutsTracklets10[bin];
    } else if (this->GetMultEstimator() == "Tracklets14") {
      valreturn = FineBinCutsTracklets14[bin];
    } else if (this->GetMultEstimator() == "TrackletsEtaGap") {
      valreturn = FineBinCutsTrackletsEtaGap[bin];
    } else {
      cout << "fMultEstimator not defined!" << '\n';
    }
  }
    if (this->GetSystem() == "AuAu_3") {
      if (this->GetMultEstimator() == "TPC") {
        valreturn = FineBinCutsTPC[bin];
      } else if (this->GetMultEstimator() == "TracksEtaGapTPC") {
        valreturn = FineBinCutsEtaGapTPC[bin];
      } else if (this->GetMultEstimator() == "Tracklets10") {
        valreturn = FineBinCutsTracklets10[bin];
      } else if (this->GetMultEstimator() == "Tracklets14") {
        valreturn = FineBinCutsTracklets14[bin];
      } else if (this->GetMultEstimator() == "TrackletsEtaGap") {
        valreturn = FineBinCutsTrackletsEtaGap[bin];
      } else {
        cout << "fMultEstimator not defined!" << '\n';
      }
    }
  return valreturn;
}

//-----------------------------------------------------

float Parameterization::GetCoarseAmplitudeCutHigh(int bin) const {
  float perc{1.0};
  perc = CoarsePercentileBinning[bin + 1];
  return perc;
}

//-----------------------------------------------------

float Parameterization::GetCoarseAmplitudeCutLow(int bin) const {
  float perc{1.0};
  perc = CoarsePercentileBinning[bin];
  return perc;
}

//-----------------------------------------------------

float Parameterization::GetFineAmplitudeCutLow(int bin) const {
  float perc{1.0};
    if (this->GetSystem() == "AuAu_200"){
        perc = CentralFinePercentileBinning[bin];
    }
    if (this->GetSystem() == "AuAu_3"){
        perc = CentralFinePercentileBinning[bin];
    }
  return perc;
}

//-----------------------------------------------------

float Parameterization::GetFineAmplitudeCutHigh(int bin) const {
  float perc{1.0};
    if (this->GetSystem() == "AuAu_200"){
        perc = CentralFinePercentileBinning[bin+1];
    }
    if (this->GetSystem() == "AuAu_3"){
        perc = CentralFinePercentileBinning[bin+1];
    }
    
  return perc;
}

//-----------------------------------------------------

int Parameterization::GetCoarseAmplitudeBinCutLowMC(int bin) const {

  int valreturn{-1};
  if (bin < 0 || bin > nCoarsePercentileBinning) {
    cout << "Error: Bin out of bounds." << '\n';
    valreturn = -1;
  }

  if (bin != nCoarsePercentileBinning - 1) {
    if (this->GetMultEstimator() == "TPC") {  //! Half TPC estimator
      if (fIsMCTrue) {
        valreturn = CoarseBinCutsTPCMCTrue[bin + 1] + 1;
      } else {
        valreturn = CoarseBinCutsTPCMCMeas[bin + 1] + 1;
      }
    } else if (this->GetMultEstimator() == "TrackletsEtaGap") {
      if (fIsMCTrue) {
        valreturn = CoarseBinCutsTrackletsEtaGapMCTrue[bin + 1] + 1;
      } else {
        valreturn = CoarseBinCutsTrackletsEtaGapMCMeas[bin + 1] + 1;
      }
    } else if (this->GetMultEstimator() == "TracksEtaGapTPC") {
      if (fIsMCTrue) {
        valreturn = CoarseBinCutsTracksEtaGapTPCMCTrue[bin + 1] + 1;
      } else {
        valreturn = CoarseBinCutsTracksEtaGapTPCMCMeas[bin + 1] + 1;
      }
    } else {
      cout << "Error: GetCoarseAmplitudeBinCutLowMC() fMultEstimator not "
              "defined!"
           << '\n';
    }
  } else {
    if (this->GetMultEstimator() == "TPC") {
      if (fIsMCTrue) {
        valreturn = CoarseBinCutsTPCMCTrue[bin + 1];
      } else {
        valreturn = CoarseBinCutsTPCMCMeas[bin + 1];
      }
    } else if (this->GetMultEstimator() == "TrackletsEtaGap") {
      if (fIsMCTrue) {
        valreturn = CoarseBinCutsTrackletsEtaGapMCTrue[bin + 1];
      } else {
        valreturn = CoarseBinCutsTrackletsEtaGapMCMeas[bin + 1];
      }
    } else if (this->GetMultEstimator() == "TracksEtaGapTPC") {
      if (fIsMCTrue) {
        valreturn = CoarseBinCutsTracksEtaGapTPCMCTrue[bin + 1];
      } else {
        valreturn = CoarseBinCutsTracksEtaGapTPCMCMeas[bin + 1];
      }
    } else {
      cout << "Error: GetCoarseAmplitudeBinCutLowMC() fMultEstimator not "
              "defined!"
           << '\n';
    }
  }

  return valreturn;
}

//-----------------------------------------------------

int Parameterization::GetCoarseAmplitudeBinCutHighMC(int bin) const {
  // Implements the bin cuts for Ntracklets/Nch to calculate efficiencies. This
  // is included for the SPD Tracklets and TPC Tracks with eta gap, and the Half
  // TPC estimator.

  int valreturn{-1};
  if (bin < 0 || bin > nCoarsePercentileBinning) {
    cout << "Error: Bin out of bounds." << '\n';
    valreturn = -1;
  }

  if (this->GetMultEstimator() == "TPC") {
    if (fIsMCTrue) {
      valreturn = CoarseBinCutsTPCMCTrue[bin];
    } else {
      valreturn = CoarseBinCutsTPCMCMeas[bin];
    }
  } else if (this->GetMultEstimator() == "TrackletsEtaGap") {
    if (fIsMCTrue) {
      valreturn = CoarseBinCutsTrackletsEtaGapMCTrue[bin];
    } else {
      valreturn = CoarseBinCutsTrackletsEtaGapMCMeas[bin];
    }
  } else if (this->GetMultEstimator() == "TracksEtaGapTPC") {
    if (fIsMCTrue) {
      valreturn = CoarseBinCutsTracksEtaGapTPCMCTrue[bin];
    } else {
      valreturn = CoarseBinCutsTracksEtaGapTPCMCMeas[bin];
    }
  } else {
    cout << "Err: GetCoarseAmplitudeBinCutHighMC() fMultEstimator not defined!"
         << '\n';
  }

  return valreturn;
}

//-----------------------------------------------------

int Parameterization::GetCoarseAmplitudeBinCutLow(int bin) const {
  int valreturn{-1};
  if (bin < 0 || bin > nCoarsePercentileBinning) {
    cout << "Error: Bin out of bounds." << '\n';
    valreturn = -1;
  }

    if (bin != nCoarsePercentileBinning - 1) {
        if (this->GetSystem() == "AuAu_200") {
            if (this->GetMultEstimator() == "TPC") {  //! Half TPC estimator
                valreturn = CoarseBinCutsTPC[bin + 1] + 1;
            } else if (this->GetMultEstimator() == "TracksEtaGapTPC") {
                valreturn = CoarseBinCutsEtaGapTPC[bin + 1] + 1;
            } else if (this->GetMultEstimator() == "Tracklets10") {
                valreturn = CoarseBinCutsTracklets10[bin + 1] + 1;
            } else if (this->GetMultEstimator() == "Tracklets14") {
                valreturn = CoarseBinCutsTracklets14[bin + 1] + 1;
            } else if (this->GetMultEstimator() == "TrackletsEtaGap") {
                valreturn = CoarseBinCutsEtaGapSPD[bin + 1] + 1;
            } else {
                cout << "fMultEstimator not defined!" << '\n';
            }
        }
        if (this->GetSystem() == "AuAu_3") {
            if (this->GetMultEstimator() == "TPC") {  //! Half TPC estimator
                valreturn = CoarseBinCutsTPC[bin + 1] + 1;
            } else if (this->GetMultEstimator() == "TracksEtaGapTPC") {
                valreturn = CoarseBinCutsEtaGapTPC[bin + 1] + 1;
            } else if (this->GetMultEstimator() == "Tracklets10") {
                valreturn = CoarseBinCutsTracklets10[bin + 1] + 1;
            } else if (this->GetMultEstimator() == "Tracklets14") {
                valreturn = CoarseBinCutsTracklets14[bin + 1] + 1;
            } else if (this->GetMultEstimator() == "TrackletsEtaGap") {
                valreturn = CoarseBinCutsEtaGapSPD[bin + 1] + 1;
            } else {
                cout << "fMultEstimator not defined!" << '\n';
            }
        }
    }
  return valreturn;
}

//-----------------------------------------------------

int Parameterization::GetCoarseAmplitudeBinCutHigh(int bin) const {
  int valreturn{-1};
  if (bin < 0 || bin > nCoarsePercentileBinning) {
    cout << "Error: Bin out of bounds." << '\n';
    valreturn = -1;
  }

  if (this->GetSystem() == "AuAu_200") {
      valreturn = 1;
     if (this->GetMultEstimator() == "TPC") {
      valreturn = CoarseBinCutsTPC[bin];
    } else if (this->GetMultEstimator() == "TracksEtaGapTPC") {
      valreturn = CoarseBinCutsEtaGapTPC[bin];
    } else if (this->GetMultEstimator() == "Tracklets10") {
      valreturn = CoarseBinCutsTracklets10[bin];
    } else if (this->GetMultEstimator() == "Tracklets14") {
      valreturn = CoarseBinCutsTracklets14[bin];
    } else if (this->GetMultEstimator() == "TrackletsEtaGap") {
      valreturn = CoarseBinCutsEtaGapSPD[bin];
    } else {
      cout << "fMultEstimator not defined!" << '\n';
    }
  }
  if (this->GetSystem() == "XeXe_544") {
    if (this->GetMultEstimator() == "V0") {
      valreturn = CoarseBinCutsXeXe544[bin];
    } else if (this->GetMultEstimator() == "TPC") {
      valreturn = CoarseBinCutsTPCXeXe544[bin];
    } else if (this->GetMultEstimator() == "TrackletsEtaGap") {
      valreturn = CoarseBinCutsTrackletsEtaGapXeXe544[bin];
    } else {
      cout << "fMultEstimator not defined!" << '\n';
    }
  }
  if (this->GetSystem() == "PbPb_276") {
    if (this->GetMultEstimator() == "V0") {
      valreturn = CoarseBinCutsPbPb276[bin];
    } else if (this->GetMultEstimator() == "TPC") {
      valreturn = CoarseBinCutsTPCPbPb276[bin];
    } else if (this->GetMultEstimator() == "TrackletsEtaGap") {
      valreturn = CoarseBinCutsTrackletsEtaGapPbPb276[bin];
    } else {
      cout << "fMultEstimator not defined!" << '\n';
    }
  }

  return valreturn;
}


//-----------------------------------------------------
const char *Parameterization::GetEstimatorLatex() const {
  const char *latex = "";
  if (this->GetMultEstimator() == "V0") {
    latex = EstimatorLatex[0];
  } else if (this->GetMultEstimator() == "TPC") {
    latex = EstimatorLatex[1];
  } else if (this->GetMultEstimator() == "Tracklets10") {
    latex = EstimatorLatex[2];
  } else if (this->GetMultEstimator() == "Tracklets14") {
    latex = EstimatorLatex[3];
  } else if (this->GetMultEstimator() == "TrackletsEtaGap") {
    latex = EstimatorLatex[4];
  } else {
    latex = EstimatorLatex[5];
  }

  return latex;
}

//-----------------------------------------------------

int Parameterization::GetnCoarsePercentileBinning() const {
  return nCoarsePercentileBinning;
}

//-----------------------------------------------------

int Parameterization::GetnFinePercentileBinning() const {
  int nn = 1;
  if (this->GetSystem() == "AuAu_200") {
    if (this->GetMultEstimator() == "TPC") {
      nn = nCentralFinePercentileBinning;
    }
  } 
    if (this->GetSystem() == "AuAu_3") {
      if (this->GetMultEstimator() == "TPC") {
        nn = nCentralFinePercentileBinning;
      }
    }
  else {
    nn = nCentralFinePercentileBinning;
  }
  return nn;
}

//-----------------------------------------------------

const char *Parameterization::GetFineAmplitudeBinLatex(int bin) const {
  const char *latex = "";

  if (this->GetSystem() == "PbPb_502") {
    if (fMultEstimator == "V0") {
      latex = FinePercentileBinningLatex[bin];
    } else {
      latex = CentralFinePercentileBinningLatex[bin];
    }
  } else {
    latex = FinePercentileBinningLatex[bin];
  }
  return latex;
}

//-----------------------------------------------------

string Parameterization::GetNameInputFile() const {
  string name = "AnalysisResults_";
  if (this->GetSystem() == "PbPb_276") {
    name = NameInputFile[0];
  } else if (this->GetSystem() == "PbPb_502") {
    name = NameInputFile[1];
  } else if (this->GetSystem() == "XeXe_544") {
    name = NameInputFile[2];
  } else {
    cout << "ERROR: Invalis fSystem, GetNameInputFile().";
  }

  return name;
}

int Parameterization::GetColor() const {
  if (this->GetSystem() == "PbPb_502") {
    return 1;
  } else if (this->GetSystem() == "PbPb_276") {
    return 2;
  } else {
    return 4;
  }

  return 0;
}

int Parameterization::GetMarker() const {
  if (this->GetSystem() == "AuAu_200") {
    return 30;
  } else if (this->GetSystem() == "PbPb_276") {
    return 27;
  } else {
    return 28;
  }

  return 0;
}

//-----------------------------------------------------

bool Parameterization::hasPointsInRange(TGraphErrors* graph, double xmin, double xmax) {
    if (!graph) return false;
    
    int nPoints = graph->GetN();
    double* x = graph->GetX();

    for (int i = 0; i < nPoints; ++i) {
        if (x[i] >= xmin && x[i] <= xmax) {
            return true; // found at least one point inside the range
        }
    }
    return false;
}

//-----------------------------------------------------

const std::vector<std::pair<double, double>>& Parameterization::GetPtRanges() {
    return pT_ranges;
}

//-----------------------------------------------------

const std::vector<std::pair<double, double>>& Parameterization::GetfitRanges() {
    return fit_ranges;
}

//-----------------------------------------------------

const std::vector<std::pair<double, double>>& Parameterization::GetsosRanges() {
    return sos_ranges;
}

//-----------------------------------------------------

const std::vector<int>& Parameterization::Getbadpoints() {
    return bad_points;
}

void Parameterization::Ratio_compute_pt_Nch(TH2D* Pt_Nch_hist, TH1D* Nch_dist, Parameterization params) {
    gROOT->SetBatch(kTRUE);
    std::string output_label = "./graphs_output.root";
    
    std::vector<Double_t> low_perc_bins, high_perc_bins;
    std::vector<int> low_bin_cuts, high_bin_cuts;
    int nperc_bins{0};
    
//    std::vector<Function> fitFunctions = {Function::LevyTsallis, Function::BW};
//    std::vector<Function> fitFunctions = {Function::BW, Function::LevyTsallis, Function::Hagedorn};
    std::vector<Function> fitFunctions = {Function::LevyTsallis};
//    std::vector<Function> fitFunctions = {Function::LevyTsallis, Function::Hagedorn};
    
    // -----**-----**-----**-----**-----**--(1) Graphs--**-----**-----**-----**-----**-----**-----**
    
    //Absolute quantities of ⟨pT⟩ vs ⟨dN/dη⟩:
    TList* lOutPt_Abs_Int = new TList(); // TList for <pT> calculation w/ pT dependence [sub-sub dir] - pT-Integrated Nch
    lOutPt_Abs_Int->SetOwner(true);
    TList* lOutPt_Abs_Ref = new TList(); // TList for <pT> calculation w/ pT dependence [sub-sub dir] - Referenced Nch
    lOutPt_Abs_Ref->SetOwner(true);
    //Normalized (to ref class) quantities of ⟨pT⟩/⟨pT⟩_{ref} vs ⟨dN/dη⟩/⟨dN/dη⟩_{ref}:
    TList* lOutPt_Norm_Int = new TList(); // TList for <pT> calculation w/ pT dependence [sub-sub dir] - pT-Integrated Nch
    lOutPt_Norm_Int->SetOwner(true);
//    TList* lOutPt_Norm_Ref = new TList(); // TList for <pT> calculation w/ pT dependence [sub-sub dir] - Referenced Nch
//    lOutPt_Norm_Ref->SetOwner(true);
    
//    TList* lOutSpectra = new TList(); // TList for spectra study
//    lOutSpectra->SetOwner(true);
    TList* lOutGen = new TList(); // TList for general elements (QA, etc)
    lOutGen->SetOwner(true);
    
    lOutGen->Add(Pt_Nch_hist);
    lOutGen->Add(Nch_dist);
    
    
    TGraphErrors* gAna = new TGraphErrors();
    gAna->SetName("Ana_vs_Ref_dNdeta");
    gAna->SetMarkerStyle(20);
    gAna->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT [StRefMultCorr]");
    gAna->GetYaxis()->SetTitle("#LTdN_{ch}/d#eta#GT [p_{T} Integrated]");
    TGraphErrors* gAnaRef = new TGraphErrors();
    gAnaRef->SetName("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%} [p_{T} Integrated]");
    gAnaRef->SetMarkerStyle(20);
    gAnaRef->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%} [StRefMultCorr]");
    gAnaRef->GetYaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%} [p_{T} Integrated]");
    TGraphErrors* gAnaDRef = new TGraphErrors();
    gAnaDRef->SetName("#LTdN_{ch}/d#eta#GT_{p_{T} Integrated}/#LTdN_{ch}/d#eta#GT_{StRefMultCorr}");
    gAnaDRef->SetMarkerStyle(20);
    gAnaDRef->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT [StRefMultCorr]");
    gAnaDRef->GetYaxis()->SetTitle("#LTdN_{ch}/d#eta#GT_{p_{T} Integrated}/#LTdN_{ch}/d#eta#GT_{StRefMultCorr}");
    TGraphErrors* gAnaDNRef = new TGraphErrors();
    gAnaDNRef->SetName("#LTdN_{ch}/d#eta#GT_{p_{T} Integrated}^{*}/#LTdN_{ch}/d#eta#GT_{StRefMultCorr}^{*}");
    gAnaDNRef->SetMarkerStyle(20);
    gAnaDNRef->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT [StRefMultCorr]");
    gAnaDNRef->GetYaxis()->SetTitle("#LTdN_{ch}/d#eta#GT_{p_{T} Integrated}^{*}/#LTdN_{ch}/d#eta#GT_{StRefMultCorr}^{*}");
    
    // -----**-----**-----**-----**--(2) Reference + knee caluclations--**-----**-----**-----**-----**-----**
    
    TFile* fOut = new TFile(output_label.c_str(),"recreate");
    
    TDirectory* FitRangesDir = fOut->mkdir("Spectra fit ranges results");
    TDirectory* ptStudyDir = FitRangesDir->mkdir("pT study");
    TDirectory* ptStudyDir_Abs = ptStudyDir->mkdir("Absolute quantities");
    TDirectory* ptStudyDir_Norm = ptStudyDir->mkdir("Normalized (to ref class) quantities");
    TDirectory* ptStudyDir_NormX = ptStudyDir->mkdir("X axis Normalized (to ref class) quantities");
    
    TDirectory* ptStudyDir_Abs_Int = ptStudyDir_Abs->mkdir("pT-Integrated Nch");
    TDirectory* ptStudyDir_Norm_Int = ptStudyDir_Norm->mkdir("pT-Integrated Nch");
    TDirectory* ptStudyDir_NormX_Int = ptStudyDir_NormX->mkdir("pT-Integrated Nch");
    
    TDirectory* spectraStudyDir = fOut->mkdir("spectra study");
    TDirectory* CentralityDir = fOut->mkdir("centrality study");
    TDirectory* MeanPTRefDir = fOut->mkdir("MeanPT_0 study");
    
    TDirectory* generalDir = fOut->mkdir("general");
    
    int fctr = 1;
    
    std::cout << std::endl; std::cout << std::endl; std::cout << std::endl;
    int total_fit_ranges = Parameterization::GetfitRanges().size();
    int fit_range_counter = 1;
    
    int eta_minutes{0};
    int eta_seconds{0};
    
    double best_fitting_score = 0.0;
    Function BestChi2fitfunctype;
    std::pair<double, double> BestChi2fitrange;
    std::pair<double, double> pTrange;
    
    
    for (const auto& fit_range : Parameterization::GetfitRanges()) { //**pT fit range loop start**
            auto start_one_fit_range = std::chrono::steady_clock::now();
            
            double fitting_score = 0.0;
            std::ostringstream subdir_name;
            subdir_name << Form("fit_range_%.2f-%.2f_GeV", fit_range.first, fit_range.second);
            TDirectory* fitRangeSubdir = ptStudyDir_Norm_Int->mkdir(subdir_name.str().c_str());
            TDirectory* fitRangeSubdirX = ptStudyDir_NormX_Int->mkdir(subdir_name.str().c_str());
            TDirectory* fitRangeSubdirAbs = ptStudyDir_Abs_Int->mkdir(subdir_name.str().c_str());
            TDirectory* spectrafitRangeSubdir = spectraStudyDir->mkdir(subdir_name.str().c_str());
        
        
            std::cout << "\033[2F";
            std::cout << "\r\033[2K";
            std::cout << "\r\033[1;34mProcessing spectra fit range: [" << fit_range.first << ", " << fit_range.second << "] (GeV/c)...\033[0m (" << fctr << "/" << Parameterization::GetfitRanges().size() << " fit ranges) ⚙️... | ETA 🕒: " << eta_minutes << " min " << eta_seconds << " sec ";
            std::cout << "\033[2B";
            std::cout << std::flush;
            fctr++;
        int ctr = 1;
        bool first = true; bool centrality_first = true;
        for (const auto& range : Parameterization::GetPtRanges()) {
            std::cout << "\033[1F";
            std::cout << "\r\033[2K";
            std::cout << "\r\033[1;32mProcessing pT range: [" << range.first << ", " << range.second << "] (GeV/c)...\033[0m (" << ctr << "/" << Parameterization::GetPtRanges().size() << " ranges) ⚙️";
            std::cout << "\033[1B";
            std::cout << std::flush;
            ctr++;
            double low_cut{range.first}; double high_cut{range.second}; // pT range selection
            
//            int ref_bin{18}; // <- you can try varying reference class here 5%->12 (OLD: 5%->17; 10%->23)
            int ref_bin{27}; // <- you can try varying reference class here 5%->12 (OLD: 5%->17; 10%->23)
//            int ref_bin{17}; // <- you can try varying reference class here 5%->12 (OLD: 5%->17; 10%->23)
            const int low_ref{params.GetFineAmplitudeBinCutLow(ref_bin)}; /// Obtaining reference bins ---> There are 49 elements (30 corresponds
            //    const int low_ref{( params.GetFineAmplitudeBinCutHigh(ref_bin) + params.GetFineAmplitudeBinCutLow(ref_bin) ) / 2}; /// Obtaining reference bins ---> There are 49 elements (30 corresponds to 5% bin):: Here is where the data needs to be pulled from .txt
//            const int high_ref{params.GetFineAmplitudeBinCutHigh(ref_bin)}; /// Obtaining reference bins:: Here is where the data needs to be pulled from .txt
            const int high_ref{params.GetFineAmplitudeBinCutHigh(0)}; /// Obtaining reference bins:: Here is where the data needs to be pulled from .txt
//            const int high_ref{params.GetFineAmplitudeBinCutHigh(ref_bin)}; /// Obtaining reference bins:: Here is where the data needs to be pulled from .txt
            TH1D* hpt_ref = (TH1D*)Pt_Nch_hist->ProjectionY("(reference) 0-5%",low_ref,high_ref);
            hpt_ref->Sumw2();
            hpt_ref->SetTitle("0-5%");
                hpt_ref->SetYTitle("1/N_{evt} d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-1}"); // pT spectra (Omar's version)
//            hpt_ref->SetYTitle("1/N_{evt} 1/(2#pi p_{T}) d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-2}");
            params.NormalizeByBinWidth(hpt_ref); //1/(2pi*pT) and 1/(dpT)
            hpt_ref->Scale(1./Nch_dist->Integral(low_ref,high_ref)); //<-- Scales spectra by 1/events
//            params.SetErrorYield(hpt_ref, 36); //Second argument is centrality index (the ref is stored in the last .txt file in yield_errors/)
//            params.SetErrorYield(hpt_ref, 29); //Second argument is centrality index (the ref is stored in the last .txt file in yield_errors/)
            params.SetErrorYield(hpt_ref, params.GetnFinePercentileBinning()); //Second argument is centrality index (the ref is stored in the last .txt file in yield_errors/)
            for (const auto& fitFuncType : fitFunctions) { //**fit function species loop start**
                std::ostringstream oss;
                TGraphErrors* gdN = new TGraphErrors();
                gdN->SetMarkerStyle(20);
                gdN->SetMarkerColor(kGreen);
                gdN->SetLineWidth(4);
                gdN->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
                gdN->GetYaxis()->SetTitle("#LTp_{T}#GT/#LTp_{T}#GT_{0-5%}");
                oss << low_cut << " #leq p_{T} #leq " << high_cut << " (GeV/c)";
                gdN->SetName(oss.str().c_str());
                oss.str(""); oss.clear();
                if (range.first == 0. && range.second == 10.){
                    oss << "STAR (Au+Au #sqrt{s_{NN}} = 200 GeV), " << "p_{T} #leq " << high_cut << " (GeV/c), |#eta| < 0.5";
                }
                else {oss << "STAR (Au+Au #sqrt{s_{NN}} = 200 GeV), " << low_cut << " #leq p_{T} #leq " << high_cut << " (GeV/c), |#eta| < 0.5";}
                gdN->SetTitle(oss.str().c_str());
                oss.str(""); oss.clear();
                
                TGraphErrors* gAN = new TGraphErrors(); //AbsY NormX (CHANGE THIS BELOW! AND FILL HISTO)
                gAN->SetMarkerStyle(20);
                gAN->SetMarkerColor(kGreen);
                gdN->SetLineWidth(2);
                gdN->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
                gdN->GetYaxis()->SetTitle("#LTp_{T}#GT/#LTp_{T}#GT_{0-5%}");
                oss << low_cut << " #leq p_{T} #leq " << high_cut << " (GeV/c)";
                gdN->SetName(oss.str().c_str());
                oss.str(""); oss.clear();
                if (range.first == 0. && range.second == 10.){
                    oss << "STAR (Au+Au #sqrt{s_{NN}} = 200 GeV), " << "p_{T} #leq " << high_cut << " (GeV/c), |#eta| < 0.5";
                }
                else {oss << "STAR (Au+Au #sqrt{s_{NN}} = 200 GeV), " << low_cut << " #leq p_{T} #leq " << high_cut << " (GeV/c), |#eta| < 0.5";}
                gdN->SetTitle(oss.str().c_str());
                oss.str(""); oss.clear();
                
                oss << low_cut << " #leq p_{T} #leq " << high_cut << " (GeV/c)";
                TGraphErrors* g = new TGraphErrors();
                g->SetMarkerStyle(20);
                g->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%} (Reference)");
                g->GetYaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%} (Integrated)");
//                g->GetYaxis()->SetTitle("#LTp_{T}#GT/#LTp_{T}#GT_{0-5%}");
                g->SetName(oss.str().c_str());
                
                TGraphErrors* gdNX = new TGraphErrors();
                gdNX->SetMarkerStyle(20);
                gdNX->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
                gdNX->GetYaxis()->SetTitle("#LTp_{T}#GT");
                gdNX->SetName(oss.str().c_str());
                
                TGraphErrors* ggAna = new TGraphErrors();
                ggAna->SetMarkerStyle(24);
                ggAna->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT");
                ggAna->GetYaxis()->SetTitle("#LTp_{T}#GT");
                ggAna->SetName(oss.str().c_str());
                
                TGraphErrors* gg = new TGraphErrors();
                gg->SetMarkerStyle(24);
                gg->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT (Reference)");
                gg->GetYaxis()->SetTitle("#LTp_{T}#GT");
                gg->SetName(oss.str().c_str());
                oss.str(""); oss.clear();
                
                TGraphErrors* gNch = new TGraphErrors();
                gNch->SetMarkerStyle(24);
                gNch->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
                gNch->GetYaxis()->SetTitle("Event fraction");
                oss << "Event fraction, " << low_cut << " #leq p_{T} #leq " << high_cut;
                gNch->SetName(oss.str().c_str());
                oss.str(""); oss.clear();
                
                TGraphErrors* gMeanPTRef = new TGraphErrors();
                gMeanPTRef->SetMarkerStyle(24);
                gMeanPTRef->GetXaxis()->SetTitle("#sqrt{s_{NN}}");
                gMeanPTRef->GetYaxis()->SetTitle("#LTp_{T}#GT_{0-5%}");
                oss << "MeanPT_0, " << Form("fit_range_%.2f-%.2f_GeV",fit_range.first, fit_range.second);
                gMeanPTRef->SetName(oss.str().c_str());
                oss.str(""); oss.clear();
                
                if (params.GetDebugSetting()){
                    cout << "Perc 0-5: " << Nch_dist->Integral(low_ref,high_ref) / Nch_dist->Integral(1,Nch_dist->GetNbinsX());
                    printf("\n\t low_ref = %d | high_ref = %d\n",low_ref,high_ref);
                    printf("\n\t Number of FINE bins: %d\n",params.GetnFinePercentileBinning());
                    cout << "\n\t Number of events: " << Nch_dist->Integral(low_ref,high_ref) << endl;
                }
                
                oss << low_cut << " #leq p_{T} #leq " << high_cut << " (GeV/c)" << Form("[%s]", Parameterization::FunctionToString(fitFuncType).c_str());
                if (params.GetDebugSetting()){ std::cout << "Processing " << Parameterization::FunctionToString(fitFuncType).c_str() << "..." << std::endl; }
                std::string name = oss.str();
                if (!gdN) {
                    std::cerr << "[ERROR] gdN is a null pointer before SetName. Skipping this entry.\n";
                    continue;
                }
                gdN->SetName(name.c_str());
                oss.str(""); oss.clear();
                
                oss << "dNdEta " << low_cut << " #leq p_{T} #leq " << high_cut << " (GeV/c) " << Form("fit_range_%.2f-%.2f_GeV",fit_range.first, fit_range.second) << Form(" [%s]", Parameterization::FunctionToString(fitFuncType).c_str());
                TGraphErrors* gCentrality = new TGraphErrors();
                gCentrality->SetMarkerStyle(24);
                gCentrality->GetXaxis()->SetTitle("Centrality (%)");
                gCentrality->GetYaxis()->SetTitle("#LTdN_{ch}/d#eta#GT");
                gCentrality->SetName(oss.str().c_str());
                oss.str(""); oss.clear();
                
                oss << "<pT> " << low_cut << " #leq p_{T} #leq " << high_cut << " (GeV/c) " << Form("fit_range_%.2f-%.2f_GeV",fit_range.first, fit_range.second) << Form(" [%s]", Parameterization::FunctionToString(fitFuncType).c_str());
                TGraphErrors* gpTCentrality = new TGraphErrors();
                gpTCentrality->SetMarkerStyle(24);
                gpTCentrality->GetXaxis()->SetTitle("Centrality (%)");
                gpTCentrality->GetYaxis()->SetTitle("#LT p_{T} #GT");
                gpTCentrality->SetName(oss.str().c_str());
                oss.str(""); oss.clear();
                
                
                if (params.GetDebugSetting()) {cout << "Fitting for 0-5%...\n";}
                TString tf1Name = Form("fit_%s", Parameterization::FunctionToString(fitFuncType).c_str());
                TF1* fit = params.GetFit(hpt_ref, fit_range.first, fit_range.second, fitFuncType, tf1Name);
                double chi2 = fit->GetChisquare();
                int ndf     = fit->GetNDF();
                double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;
//                if (chi2ndf > 1000. || chi2ndf < 1e-9) {
//                    if (params.GetDebugSetting()) std::cout << "Poor fit on reference bin. Skipping..." << std::endl;
//                    continue;
//                }
                fitting_score += params.ComputeFitScore(chi2ndf);
                
                TString funcStr = Parameterization::FunctionToString(fitFuncType);
                TString label = Form("D/F 0-5%% [%s]", funcStr.Data());
                hpt_ref->SetName(Form("(ref) 0-5%% [%s]", funcStr.Data()));
                
                if (params.GetDebugSetting()){ std::cout << "Attempting " << std::string(funcStr.Data()) << " fit..." << std::endl; }
                
                    
                const double meanpt_ref{params.GetMean(hpt_ref,fit,true,high_cut,low_cut)};
                const double meannch_ref{params.GetNch(Nch_dist,low_ref,high_ref)};
                const double dndetaAna_ref{params.GetdNdEta(hpt_ref,fit,true,high_cut,low_cut)};
                const double error_ref{params.GetErrorMeanpT(hpt_ref,fit,true,high_cut,params.GetDebugSetting())};
                
                if (params.GetDebugSetting()){
                    cout << "\n( Before event & pT normalizations ): " << params.GetdNdEta(hpt_ref,fit,true,0.2,0);
                    cout << '\n' << "<pT>_{0-5%}: " << meanpt_ref;
                    cout << '\n' << "<dn/deta>_{0-5%, Mean Nch}: " << meannch_ref;
                    cout << '\n' << "<dn/deta>_{0-5%, Integrated}: " << dndetaAna_ref;
                }
                
                if (first && params.GetDisplayFitRangesSetting()){
                    double RootS = 200;
                    gMeanPTRef->SetPoint(0, RootS, meanpt_ref);
                    gMeanPTRef->SetPointError(0, 0.0, error_ref);
                    gdN->GetXaxis()->SetRangeUser(20.0, 300.0);
                    gdN->GetYaxis()->SetRangeUser(meanpt_ref - 0.2, meanpt_ref + 0.2);
                    MeanPTRefDir->cd();
                    gMeanPTRef->Write();
                    
                    spectrafitRangeSubdir->cd();
//                    TH1* hFitRatio_ref = params.GetDataFitRatio(hpt_ref, fit, "D/F 0-5%");
//                    params.GetDataFitRatio(hpt_ref, fit, std::string(label.Data()));
//                    std::cout << "Writing reference spectrum: " << hpt_ref->GetName() << std::endl;
                    hpt_ref->Write();
                }
                
//                for (int i = 1; i < Nch_dist->GetNbinsX() - 5; i += 5) {
//                    
//                    if (params.GetDebugSetting()){ std::cout << "Computing general event fractions...\n";}
//                    
//                    TH1D* hpt_n = nullptr;
//                    
//                    if (Nch_dist->Integral( i, i + 5 ) < 10) continue;
//                    
//                    hpt_n = Pt_Nch_hist->ProjectionY(
//                    Form("hpt_n_%d_%d", i, i+5),i, i+5);
//                    
//                    fit = params.GetFit(hpt_n, fit_range.first, fit_range.second, fitFuncType, tf1Name, i);
//                    
//                    double dndetaAna_n{params.GetdNdEta(hpt_n, fit, true, high_cut,low_cut)};
//                    double error_dndetaAna_n{params.GetErrordNdEta(hpt_n, fit, true, 0.2)};
//                    double error_dndetaAna_ref{params.GetErrordNdEta(hpt_ref, fit, true, 0.2)};
//                    
//                    gNch->SetPoint(
//                    i, dndetaAna_n / dndetaAna_ref, Nch_dist->Integral( i, i + 5 ) / Nch_dist->Integral());
//                    
//                    double error_dndeta_ratio_n{( dndetaAna_n / dndetaAna_ref) * sqrt(pow(( (error_dndetaAna_n / dndetaAna_n) - (error_dndetaAna_ref / dndetaAna_ref) ), 2) )};
//                    
////                    gNch->SetPointError(i, error_dndeta_ratio_n, 1e-6);
//                    gNch->SetPointError(i, 0.0, 0.);
////                    gNch->GetXaxis()->SetRangeUser(0.8, 1.25);
////                    gNch->GetYaxis()->SetRangeUser(10e-6, 1);
//                }
//                CentralityDir->cd();
//                gNch->Write();
                
                // -----**-----**-----**-----**-----**--(3) Bin Calculations--**-----**-----**-----**-----**-----**-----**
                nperc_bins = params.GetnFinePercentileBinning();
                double iter = 0.0;
                for (int k = 0; k < nperc_bins; k++) {
                    iter += 1.;
                    bool isBad = false;
                    for (const auto& bad_point : Parameterization::Getbadpoints()) {
                        if (k == bad_point){
                            isBad = true;
                            break;
                        }
                    }
                    if (isBad) continue;
                    
                    if (params.GetDebugSetting()){
                        cout << "\nparams.GetFineAmplitudeCutLow(k)*100: " << params.GetFineAmplitudeCutLow(k)*100;
                        cout << "\nparams.GetFineAmplitudeCutHigh(k)*100: " << params.GetFineAmplitudeCutHigh(k)*100;
                        cout << "\nparams.GetFineAmplitudeBinCutLow(k): " << params.GetFineAmplitudeBinCutLow(k);
                        cout << "\nparams.GetFineAmplitudeBinCutHigh(k): " << params.GetFineAmplitudeBinCutHigh(k) << endl;
                    }
                    
                    TH1D* hpt = nullptr;
                    
                    if (params.GetCentralitySetting() == "I"){ //Inclusive Centrality Analysis
                        hpt = Pt_Nch_hist->ProjectionY(Form("%.3f-%.3f%% [%s]", params.GetFineAmplitudeCutLow(0)*100,
                                                            params.GetFineAmplitudeCutHigh(k)*100,
                                                            funcStr.Data()),
                                                       params.GetFineAmplitudeBinCutLow(k),
                                                       params.GetFineAmplitudeBinCutHigh(0));
                        hpt->SetTitle(Form("%.0f-%.4f%%", params.GetFineAmplitudeCutLow(0)*100,
                                           params.GetFineAmplitudeCutHigh(k)*100));
                    }
                    else{                                      //Exclusive Centrality Analysis
                        hpt = Pt_Nch_hist->ProjectionY(Form("%.3f-%.3f%% [%s]", params.GetFineAmplitudeCutLow(k)*100,
                                                            params.GetFineAmplitudeCutHigh(k)*100,
                                                            funcStr.Data()),
                                                       params.GetFineAmplitudeBinCutLow(k),
                                                       params.GetFineAmplitudeBinCutHigh(k));
                        hpt->SetTitle(Form("%.3f-%.3f%%",
                                           params.GetFineAmplitudeCutLow(k)*100,
                                           params.GetFineAmplitudeCutHigh(k)*100));
                    }
                    
                    hpt->Sumw2();
//                    hpt->SetYTitle("1/N_{evt} 1/(2#pi p_{T}) d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-2}"); // Rutik/Rene's version (differential yield)
                    hpt->SetYTitle("1/N_{evt} d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-1}"); // Omar's version (spectra)
                    params.NormalizeByBinWidth(hpt);
                    if (params.GetCentralitySetting() == "I"){ //Inclusive Centrality Analysis
                        hpt->Scale(1./Nch_dist->Integral(params.GetFineAmplitudeBinCutLow(k),
                                                         params.GetFineAmplitudeBinCutHigh(0))); //Inclusive
                    }
                    else{                                      //Exclusive Centrality Analysis
                        hpt->Scale(1./Nch_dist->Integral(params.GetFineAmplitudeBinCutLow(k),
                                                         params.GetFineAmplitudeBinCutHigh(k))); //Exclusive
                    }
                    params.SetErrorYield(hpt, k);
                    
                    if (params.GetDebugSetting()){
                        std::cout << "\n Setting error in pT spectra: " << std::endl;
                        for (int bin = 1; bin < 10; bin++) {
                            std::cout << " [" << bin << "] error: " << hpt->GetBinError(bin);
                        }
                        std::cout << std::endl;
                    }
                    
                    bool ifpt_cut = true; //Enforcing 10.0 GeV and 200 MeV cutoff
                    tf1Name = Form("hfit_%s", Parameterization::FunctionToString(fitFuncType).c_str());
                    if (params.GetDebugSetting()) std::cout << "Going into fit" << std::endl;
                    fit = params.GetFit(hpt, fit_range.first, fit_range.second, fitFuncType, tf1Name, k);
                    //                    fit = params.GetFit(hpt, fit_range.first, fit_range.second, fitFuncType);
                    if (params.GetDebugSetting()) std::cout << "Fit completed" << std::endl;
                    const double xMinEx = 0.0; const double xMaxEx = 0.2; const int nEx    = 4;
                    TH1D* hptEx = nullptr;
                    
                    chi2 = fit->GetChisquare();
                    ndf     = fit->GetNDF();
                    chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;
                    if (k < 20 && (chi2ndf > 250. || chi2ndf < 1e-9)) continue;
                    fitting_score += params.ComputeFitScore(chi2ndf);
                    
                    if (params.GetDebugSetting()) {
                        std::cout << "Writing hpt";
                        cout << "\n Index: " << k;
                        cout << "\n first: " << first << std::endl;
                    }
                    
                    if (first && params.GetDisplayFitRangesSetting()){
                        spectrafitRangeSubdir->cd();
                        label = Form("D/F %.4f-%.4f%% [%s]",
                                     params.GetFineAmplitudeCutLow(k)*100,
                                     params.GetFineAmplitudeCutHigh(k)*100,
                                     funcStr.Data());
//                        params.GetDataFitRatio(hpt, fit, std::string(label.Data()));
//                        hpt->SetMarkerStyle(4); hpt->SetMarkerColor(kBlack); hpt->SetMarkerSize(1.2);
                        hpt->SetMarkerStyle(20); hpt->SetMarkerColor(kBlack); hpt->SetMarkerSize(1.2);
//                        hptEx->SetLineColor(0); hptEx->SetLineWidth(0);
                        hptEx = new TH1D( Form("%s_extrap", hpt->GetName()), "Extrapolated points",
                          nEx, xMinEx, xMaxEx
                        );
                        for (int b = 1; b <= nEx; ++b) {
                          const double x = hptEx->GetXaxis()->GetBinCenter(b);
                          const double y = fit->Eval(x);
                          hptEx->SetBinContent(b, y);
                          hptEx->SetBinError(b, 0.0);
                        }
                        hptEx->SetMarkerStyle(24); hptEx->SetMarkerColor(kBlack); hptEx->SetMarkerSize(1.2);
                        hpt->Write();
                        hptEx->Write();
                    }
                    //                    lOutSpectra->Add(hpt);
                    //                lOutSpectra->Add(params.GetSpectraDifference(hpt, hpt_ref, Form("Relative Sp. Diff. ((%.2f-%.2f%%) - (0-5%%))/(0-5%%)",
                    //                                                                         params.GetFineAmplitudeCutLow(0)*100,
                    //                                                                         params.GetFineAmplitudeCutHigh(k)*100)));
                    //                lOutSpectra->Add(params.GetSpectraDifference(hpt, hpt_ref, Form("Relative Sp. Diff. ((%.3f-%.3f%%) - (0-5%%))/(0-5%%)", params.GetFineAmplitudeCutLow(k)*100, params.GetFineAmplitudeCutHigh(k)*100)));
                    //                    lOutSpectra->Add(params.GetDataFitRatio(hpt, fit, Form("D/F %.3f-%.3f%%",
                    //                                                                           params.GetFineAmplitudeCutLow(0)*100,
                    //                                                                           params.GetFineAmplitudeCutHigh(k)*100)));
                    //                lOutSpectra->Add(params.GetDataFitRatio(hpt, fit, Form("D/F %.2f-%.2f%%",
                    //                                                                            params.GetFineAmplitudeCutLow(k)*100,
                    //                                                                            params.GetFineAmplitudeCutHigh(k)*100)));
                    
                        
                    if (params.GetDebugSetting()) std::cout << "Calculating mean pT, etc." << std::endl;
                    double mean_pt{params.GetMean(hpt,fit,ifpt_cut,high_cut,low_cut)};
                    double dndetaAna{params.GetdNdEta(hpt, fit, ifpt_cut, high_cut,low_cut)};

                    double error_pt{params.GetErrorMeanpT(hpt,fit,true,high_cut, params.GetDebugSetting())};
//                    if (error_pt > 0.0025 || error_ref > 0.0025) { if (params.GetDebugSetting()) {std::cout << "SKIPPING from pT error";} continue; }
//                    if (error_pt > 0.005 || error_ref > 0.005) { if (params.GetDebugSetting()) {std::cout << "SKIPPING from pT error";} continue; }
                    //For reference, largest CMS stat error is 0.002
                    //            if (k == 2) error_pt = params.GetErrorMeanpT(hpt_ref,fit,true,low_cut, true);
                    double mean_nch;
                    if (params.GetCentralitySetting() == "I"){
                        mean_nch = params.GetNch(Nch_dist,params.GetFineAmplitudeBinCutLow(0),params.GetFineAmplitudeBinCutHigh(k));
                    }
                    else { mean_nch = params.GetNch(Nch_dist,params.GetFineAmplitudeBinCutLow(k),params.GetFineAmplitudeBinCutHigh(k)); }
                    
                    
                    double error_dndetaAna{params.GetErrordNdEta(hpt, fit, true, 0.2)};
                    double error_dndetaAna_ref{params.GetErrordNdEta(hpt_ref, fit, true, 0.2)};
//                    if (error_dndetaAna > 50 || error_dndetaAna_ref > 50) {std::cout << "SKIPPING from dNdeta error"; continue;}
                    double error_nch{params.GetRelativeError(k)};
                    double error_nch_ref{params.GetRelativeError(ref_bin)};
                    
//                    double error_ratio{(mean_pt / meanpt_ref) * sqrt(pow(( (error_pt / mean_pt) - (error_ref / meanpt_ref) ), 2.) )};
                    double error_ratio{(mean_pt / meanpt_ref) * sqrt( pow( (error_pt / mean_pt), 2.) + pow ((error_ref / meanpt_ref) , 2.) )};
                    // ^Error progation-> https://en.wikipedia.org/wiki/Propagation_of_uncertainty
                    //                    https://arxiv.org/pdf/1609.04150#page24
                    double error_nch_ratio{(mean_nch / meannch_ref) * sqrt(pow( (error_nch / mean_nch), 2.0 ) +  pow( (error_nch_ref / meannch_ref), 2.0 ) ) };
                    // ^Error progation-> https://en.wikipedia.org/wiki/Propagation_of_uncertainty
                    //                    https://arxiv.org/pdf/1609.04150#page24
                    double error_dndeta_ratio{( dndetaAna / dndetaAna_ref) * sqrt(pow(( (error_dndetaAna / dndetaAna) - (error_dndetaAna_ref / dndetaAna_ref) ), 2) )};
                    // ^Error progation-> https://en.wikipedia.org/wiki/Propagation_of_uncertainty
                    //                    https://arxiv.org/pdf/1609.04150#page24
                    
//                    error_ratio = error_ratio + 0.01*error_ratio; //sys uncertainty on mean pT ratio
//                     1% relative systematic estimated from /Users/calebbroodo/aux_sos_analysis/systematics/Code
//                    error_dndeta_ratio = error_dndeta_ratio + 0.008*error_dndeta_ratio; //sys uncertainty on dNdeta ratio
                    
                    if (params.GetDebugSetting()) {
                        cout << "\n Index: " << k;
                        cout << "\n Number of events: " << Nch_dist->Integral(params.GetFineAmplitudeBinCutLow(k),params.GetFineAmplitudeBinCutHigh(0));
                        cout << "\n params.GetFineAmplitudeBinCutHigh(k): " << params.GetFineAmplitudeBinCutHigh(0);
                        cout << "\n params.GetFineAmplitudeBinCutLow(k): " << params.GetFineAmplitudeBinCutLow(k);
                        cout << "\n dNdeta / dNdeta_0: " << dndetaAna / dndetaAna_ref;
                        cout << "\n mean_pt/meanpt_ref: " << mean_pt/meanpt_ref;
                        cout << "\n chi^2/ndf: " << mean_pt/meanpt_ref;
                        std::cout << "Setting points." << std::endl;
                        std::cout << "dndetaAna: " << dndetaAna << std::endl;
                        std::cout << "dndetaAna_ref: " << dndetaAna_ref << std::endl;
                        std::cout << "mean_pt: " << mean_pt << std::endl;
                        std::cout << "meanpt_ref: " << meanpt_ref << std::endl;
                        std::cout << "(mean_pt/meanpt_ref)" << (mean_pt/meanpt_ref);
                        std::cout << "(dndetaAna / dndetaAna_ref)" << (dndetaAna / dndetaAna_ref);
                        std::cout << "error_dndeta_ratio" << error_dndeta_ratio;
                        std::cout << "error_ratio" << error_ratio;
                    }
                    
                    if (k > 1 && k < 20){ //Checks if points are monotonic within errors
                        double prev_x, prev_y;
                        gdN->GetPoint(k - 1, prev_x, prev_y);
                        double prev_x_err = gdN->GetErrorX(k - 1);
                        double prev_y_err = gdN->GetErrorY(k - 1);
                        if (!params.IsMonotonicWithinError(prev_x, prev_y,
                                                    dndetaAna / dndetaAna_ref, mean_pt/meanpt_ref,
                                                    prev_x_err, prev_y_err,
                                                    error_dndeta_ratio, error_ratio)) continue;
                    }
                    
                    if (params.GetDebugSetting()) std::cout << "Setting points gdN." << std::endl;
                    gdN->SetPoint(k, dndetaAna / dndetaAna_ref, mean_pt/meanpt_ref);
                    gdN->SetPointError(k, error_dndeta_ratio, error_ratio);
                    //Ignore large error bars for high pT and high centrality bins
                    //            if (low_cut > 0.4) {gdN->SetPointError(k, 0, error_ratio);}
                    gdN->GetXaxis()->SetRangeUser(0.8, 1.25);
                    gdN->GetYaxis()->SetRangeUser(0.99, 1.025);
                    
                    
                    if (params.GetDebugSetting()) std::cout << "Setting points gdNX." << std::endl;
                    gdNX->SetPoint(k, dndetaAna / dndetaAna_ref, mean_pt);
                    gdNX->SetPointError(k, error_dndeta_ratio, error_ratio);
                    gdNX->GetXaxis()->SetRangeUser(0.8, 1.25);
                    gdNX->GetYaxis()->SetRangeUser(0.5, 0.7);
                    
                    
                    if (params.GetDebugSetting()) std::cout << "Setting points g." << std::endl;
//                    g->SetPoint(k, mean_nch / meannch_ref, mean_pt / meanpt_ref);
                    g->SetPoint(k, mean_nch / meannch_ref, dndetaAna / dndetaAna_ref);
                    g->SetPointError(k, error_nch_ratio, error_ratio);
                    g->GetXaxis()->SetRangeUser(0.8, 1.25);
                    g->GetYaxis()->SetRangeUser(0.8, 1.25);
//                    g->GetYaxis()->SetRangeUser(0.98, 1.025);
                    
                    
                    if (params.GetDebugSetting()) std::cout << "Setting points gg." << std::endl;
                    gg->SetPoint(k, mean_nch, mean_pt);
                    gg->SetPointError(k, error_nch, error_pt);
                    //            gg->GetXaxis()->SetRangeUser(400, 650);
                    //            gg->GetYaxis()->SetRangeUser(0.55, 0.57);
                    
                    
                    if (params.GetDebugSetting()) std::cout << "Setting points ggAna." << std::endl;
                    ggAna->SetPoint(k, dndetaAna, mean_pt);
                    ggAna->SetPointError(k, error_dndetaAna, error_pt);
                    //            ggAna->GetXaxis()->SetRangeUser(200, 1000);
                    //            ggAna->GetYaxis()->SetRangeUser(0.55, 0.57);
                    //            glog->SetPoint(k, log(mean_nch), log(mean_pt));
                    //            glog->SetPointError(k, error_nch_ratio/mean_nch, error_pt/mean_pt);
                    
                    if (params.GetDebugSetting()) std::cout << "Setting points gAna." << std::endl;
                    gAna->SetPoint(k, mean_nch, dndetaAna );
                    gAna->SetPointError(k, error_nch, error_dndetaAna);
                    
                    if (params.GetDebugSetting()) std::cout << "Setting points gAnaRef." << std::endl;
                    gAnaRef->SetPoint(k, mean_nch/meannch_ref, dndetaAna / dndetaAna_ref);
                    gAnaRef->SetPointError(k, error_nch_ratio, error_dndeta_ratio);
                    gAnaRef->GetXaxis()->SetRangeUser(0.2, 1.25);
                    
                    if (params.GetDebugSetting()) std::cout << "Setting points gAnaDRef." << std::endl;
                    gAnaDRef->SetPoint(k, mean_nch, dndetaAna / mean_nch);
                    gAnaDRef->SetPointError(k, error_nch_ratio, 0);
                    gAnaDRef->GetXaxis()->SetRangeUser(100, 800);
                    
                    if (params.GetDebugSetting()) std::cout << "Setting points gAnaDNRef." << std::endl;
                    gAnaDNRef->SetPoint(k, mean_nch, (dndetaAna / dndetaAna_ref) / (mean_nch / meannch_ref));
                    gAnaDNRef->SetPointError(k, 0, 0);
                    gAnaDNRef->GetXaxis()->SetRangeUser(100, 800);
                    
                    if(centrality_first){
                        if (params.GetDebugSetting()) std::cout << "Setting points Centrality." << std::endl;
                        double cent{params.GetCentrality(k, params, false)};
                        gCentrality->SetPoint(k, cent, dndetaAna);
                        gCentrality->SetPointError(k, 0, error_dndetaAna);
                        gpTCentrality->SetPoint(k, cent, mean_pt);
                        gpTCentrality->SetPointError(k, 0, error_pt);
                    }
                    
                    if (params.GetDebugSetting()) std::cout << "Finished setting points" << std::endl;
                    
                }
                
                if (params.GetDebugSetting()) {cout << "\n#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#\n" << endl;}
                
                if (centrality_first && params.GetDisplayFitRangesSetting()){
                    CentralityDir->cd();
                    gCentrality->Write();
//                    gNch->Write();
                    gpTCentrality->Write();
//                    delete gNch;
                    delete gCentrality;
                }
                centrality_first = false;
                
                fitRangeSubdir->cd();
                for (const auto& sos_range : Parameterization::GetsosRanges()){
                    if (gdN && gdN->GetN() > 0 && hasPointsInRange(gdN, sos_range.first, sos_range.second)) {
                        std::cout << "\r\033[2K";
                        std::cout << "\r\033[1;31mProcessing " << Parameterization::GetsosRanges().size()
                        << " speed of sound fit ranges: ["
                        << Parameterization::GetsosRanges().front().first << ", " << Parameterization::GetsosRanges().front().second
                        << "] - ["
                        << Parameterization::GetsosRanges().back().first << ", " << Parameterization::GetsosRanges().back().second
                        << "] (GeV/c)...\033[0m⚙️" << std::flush;
                        
                        if (params.GetDebugSetting()) std::cout << "Performing sound fit" << std::endl;
                        TGraphErrors* gdN_clone = params.GetSoundFit(gdN, sos_range.first, sos_range.second, lOutPt_Norm_Int, low_cut, high_cut);
                        if (params.GetDebugSetting()) std::cout << "cloning gdN" << std::endl;
                        if (gdN_clone && params.GetDisplayFitRangesSetting()){
                            gdN_clone->Write();
                            delete gdN_clone;
                        }
                    }
                } //sos fit range loop end
                
                if (params.GetDebugSetting()) std::cout << "Writing gdN" << std::endl;
                if (params.GetDisplayFitRangesSetting()) gdN->Write();
                delete gdN;
                
                if (params.GetDebugSetting()) std::cout << "Writing gdNX" << std::endl;
                fitRangeSubdirX->cd();
                if (params.GetDisplayFitRangesSetting()) gdNX->Write();
                delete gdNX;
                
                if (params.GetDebugSetting()) std::cout << "Writing ggAna" << std::endl;
                fitRangeSubdirAbs->cd();
                if (params.GetDisplayFitRangesSetting()) ggAna->Write();
                delete ggAna;
                
                if (params.GetDebugSetting()) std::cout << "Writing g." << std::endl;
                CentralityDir->cd();
                if (params.GetDisplayFitRangesSetting()) g->Write();
                delete g;
                
                fitting_score /= (iter + 1);
                if (fitting_score > best_fitting_score){
                    if (params.GetDebugSetting()) std::cout << "\n Fitting score: " << fitting_score << std::endl;
                    best_fitting_score = fitting_score;
                    BestChi2fitfunctype = fitFuncType;
                    BestChi2fitrange = fit_range;
                    pTrange = range;
                }
                
                if (params.GetDebugSetting()) std::cout << "End of fit function loop" << std::endl;
            } //fit function loop end
            first = false;
        } //pT fit range loop end
        
        if (params.GetDebugSetting()) std::cout << "the best fitting score was " << best_fitting_score << std::endl;
        
        auto end_one_fit_range = std::chrono::steady_clock::now();
        auto elapsed_one_fit_range = std::chrono::duration_cast<std::chrono::seconds>(end_one_fit_range - start_one_fit_range).count();

        int remaining_fit_ranges = total_fit_ranges - fit_range_counter;
        int estimated_remaining_time_sec = elapsed_one_fit_range * remaining_fit_ranges;

        eta_minutes = estimated_remaining_time_sec / 60;
        eta_seconds = estimated_remaining_time_sec % 60;

        fit_range_counter++;
    } //spectra fit range end

    // -----**-----**-----**-----**-----**--(5)--**-----**-----**-----**-----**-----**-----**
    std::cout << "\n Writing files for the spectra range with the best fitting score (" << best_fitting_score << ")..." << std::endl;
//    params.WriteBestChi2(BestChi2fitfunctype, BestChi2fitrange, pTrange, fOut, Pt_Nch_hist, Nch_dist, params);
    
    std::cout << "\nWriting final directories..." << std::endl;
    
    generalDir->cd();
    lOutGen->Write();
    fOut->cd();
    fOut->Write();
    
    std::cout << "deleting files..." << std::endl;
    delete fOut;
//    delete lOutPt_Norm_Int;
//    delete lOutPt_Abs_Int; delete lOutPt_Abs_Ref; delete lOutPt_Norm_Ref;
//    delete lOutSpectra; delete lOutGen;
    gROOT->SetBatch(kFALSE);
}

// -----------------------------------------------------


void Parameterization::saveToCSV(int nperc_bins, std::vector<int> high_bin_cuts, std::vector<int> low_bin_cuts, std::vector<double> low_perc_bins, std::vector<double> high_perc_bins, TH1* Nch_dist, Parameterization params) {
    // Open the CSV file
    std::ofstream csvFile("Event_Numerics.csv");

    // Check if the file is open
    if (!csvFile.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Write the headers to the file
    csvFile << "Bin Index,Bin Cuts,Bin Cut Range,Centrality,N Events, Relative Error (Evt fraction)\n";

    // Write the data in columns
    for (int k = 1; k < nperc_bins; k++) {
        // bin index
        std::string bin_index = std::to_string(k-1);
        
        // bin cuts
        std::string bin_cuts = std::to_string(high_bin_cuts[k-1]) + " - " + std::to_string(low_bin_cuts[k-1]);
        
        // bin cut range
        int range = abs (high_bin_cuts[k-1] - low_bin_cuts[k-1]);
        std::string bin_range = std::to_string( range );

        // centrality
        std::stringstream centralityStream;
                centralityStream << std::fixed << std::setprecision(3)
                                 << low_perc_bins[k-1] * 100 << "% - "
                                 << high_perc_bins[k-1] * 100 << "%";
                std::string centrality = centralityStream.str();
        
        // relative error
        double error = params.GetRelativeError(k-1);
        std::string bin_error = std::to_string( error );

        // N Events using the integral
        double n_events = Nch_dist->Integral(params.GetFineAmplitudeBinCutLow(k-1), params.GetFineAmplitudeBinCutHigh(k-1));

        // Write the row to the CSV file
        csvFile << bin_index << "," << bin_cuts << "," << bin_range << "," << centrality << "," << n_events << "," << error << "\n";
    }
    
    // Close the file
    csvFile.close();
}

// -----------------------------------------------------
void Parameterization::PrintFineBinCuts() {
    cout << "Reading stored Fine Bin Cut values... " << '\n';
    for (size_t i = 0; i < (this->GetnFinePercentileBinning() + 1); ++i) {
        std::cout << this->GetFineAmplitudeBinCutHigh(i) << " | ";
        std::cout << this->GetFineAmplitudeBinCutLow(i) << ", ";
//        std::cout << FineBinCutsTPC[i] << ", ";
    }
}

// -----------------------------------------------------

std::string Parameterization::FunctionToString(Function f) {
    switch (f) {
        case Function::LevyTsallis:       return "LevyTsallis";
        case Function::Tsallis:           return "Tsallis";
        case Function::Hagedorn:          return "Hagedorn";
        case Function::TBW:               return "TBW";
        case Function::BW:                return "BW";
        case Function::FokkerPlanck:      return "FokkerPlanck";
        case Function::PowerLaw:          return "PowerLaw";
        case Function::ExpDecay:          return "ExpDecay";
        case Function::RightKneeFunction: return "RightKneeFunction";
        case Function::RightKneeFunctionApproximation: return "RightKneeFunctionApproximation";
        default:                          return "UnknownFunction";
    }
}
// -----------------------------------------------------

bool Parameterization::IsMonotonicWithinError(double prev_x, double prev_y,
                            double curr_x, double curr_y,
                            double prev_x_err, double prev_y_err,
                            double curr_x_err, double curr_y_err) {
    double tolerance = 0.8;
    bool debug = false;
    double combined_x_err =  std::sqrt(prev_x_err * prev_x_err + curr_x_err * curr_x_err);
    double combined_y_err =  std::sqrt(prev_y_err * prev_y_err + curr_y_err * curr_y_err);

    bool x_increasing = curr_x > tolerance * (prev_x - combined_x_err);
    bool y_increasing = curr_y > tolerance * (prev_y - combined_y_err);

    if (debug){
            if (!x_increasing || !y_increasing) {
                std::cerr << "[Reject] Monotonic check failed:"
                          << " prev_x=" << prev_x << " ± " << prev_x_err
                          << ", curr_x=" << curr_x << " ± " << curr_x_err
                          << " | prev_y=" << prev_y << " ± " << prev_y_err
                          << ", curr_y=" << curr_y << " ± " << curr_y_err
                          << std::endl;
            }
    }


    return x_increasing && y_increasing;
}

// -----------------------------------------------------

double Parameterization::ComputeFitScore(double chi2ndf){
    double score = 0.0;
    double tolerance = 1000.0;
    
    if (chi2ndf > 10000 || chi2ndf < 1e-9) return score;
    
    score = exp( -1 * pow((chi2ndf - 1) / tolerance, 2) );
    
    return score;
}

// -----------------------------------------------------

void Parameterization::WriteBestChi2(Parameterization::Function BestChi2fitfunctype, std::pair<double, double> BestChi2fitrange, std::pair<double, double> pTrange, TFile* fOut, TH2D* Pt_Nch_hist, TH1D* Nch_dist, Parameterization params){
    std::cout << "\nSaving results from best fitting score on pT spectra..." << std::endl;
        gROOT->SetBatch(kTRUE);
        Function fitFuncType = BestChi2fitfunctype;
        std::pair<double, double> fit_range = BestChi2fitrange;
        int nperc_bins{0};
        
        // -----**-----**-----**-----**-----**--(1) Graphs--**-----**-----**-----**-----**-----**-----**
        
        TList* lOutPt_Norm_Int = new TList(); // TList for <pT> calculation w/ pT dependence [sub-sub dir] - pT-Integrated Nch
        lOutPt_Norm_Int->SetOwner(true);
        TList* lOutGen = new TList(); // TList for general elements (QA, etc)
        lOutGen->SetOwner(true);
        
        lOutGen->Add(Pt_Nch_hist);
        lOutGen->Add(Nch_dist);
        
        TGraphErrors* gAna = new TGraphErrors();
        gAna->SetName("Ana_vs_Ref_dNdeta");
        gAna->SetMarkerStyle(20);
        gAna->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT [StRefMultCorr]");
        gAna->GetYaxis()->SetTitle("#LTdN_{ch}/d#eta#GT [p_{T} Integrated]");
        TGraphErrors* gAnaRef = new TGraphErrors();
        gAnaRef->SetName("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%} [p_{T} Integrated]");
        gAnaRef->SetMarkerStyle(20);
        gAnaRef->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%} [StRefMultCorr]");
        gAnaRef->GetYaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%} [p_{T} Integrated]");
        TGraphErrors* gAnaDRef = new TGraphErrors();
        gAnaDRef->SetName("#LTdN_{ch}/d#eta#GT_{p_{T} Integrated}/#LTdN_{ch}/d#eta#GT_{StRefMultCorr}");
        gAnaDRef->SetMarkerStyle(20);
        gAnaDRef->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT [StRefMultCorr]");
        gAnaDRef->GetYaxis()->SetTitle("#LTdN_{ch}/d#eta#GT_{p_{T} Integrated}/#LTdN_{ch}/d#eta#GT_{StRefMultCorr}");
        TGraphErrors* gAnaDNRef = new TGraphErrors();
        gAnaDNRef->SetName("#LTdN_{ch}/d#eta#GT_{p_{T} Integrated}^{*}/#LTdN_{ch}/d#eta#GT_{StRefMultCorr}^{*}");
        gAnaDNRef->SetMarkerStyle(20);
        gAnaDNRef->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT [StRefMultCorr]");
        gAnaDNRef->GetYaxis()->SetTitle("#LTdN_{ch}/d#eta#GT_{p_{T} Integrated}^{*}/#LTdN_{ch}/d#eta#GT_{StRefMultCorr}^{*}");
        // -----**-----**-----**-----**--(2) Reference + knee caluclations--**-----**-----**-----**-----**-----**
                
                TDirectory* Bestchi2Dir = fOut->mkdir("Best Chi2Ndf on spectra");
                TDirectory* Bestchi2ptStudyDir = Bestchi2Dir->mkdir("pT study");
                TDirectory* Bestchi2ptStudyDir_Abs = Bestchi2ptStudyDir->mkdir("Absolute quantities");
                TDirectory* Bestchi2ptStudyDir_Norm = Bestchi2ptStudyDir->mkdir("Normalized (to ref class) quantities");
                TDirectory* Bestchi2ptStudyDir_NormX = Bestchi2ptStudyDir->mkdir("X axis Normalized (to ref class) quantities");
                TDirectory* Bestchi2spectraStudyDir = Bestchi2Dir->mkdir("spectra study");
                TDirectory* Bestchi2CentralityDir = Bestchi2Dir->mkdir("centrality study");
                std::ostringstream subdir_name;
    std::cout << "Processing fit range: " << fit_range.first << ", " << fit_range.second << std::endl;
                subdir_name << Form("fit_range_%.2f-%.2f_GeV", fit_range.first, fit_range.second);
                TDirectory* fitRangeSubdir = Bestchi2ptStudyDir_Norm->mkdir(subdir_name.str().c_str());
                TDirectory* fitRangeSubdirX = Bestchi2ptStudyDir_NormX->mkdir(subdir_name.str().c_str());
                TDirectory* fitRangeSubdirAbs = Bestchi2ptStudyDir_Abs->mkdir(subdir_name.str().c_str());
                TDirectory* spectrafitRangeSubdir = Bestchi2spectraStudyDir->mkdir(subdir_name.str().c_str());
                std::cout << "\033[2F";
                std::cout << "\r\033[2K";
                std::cout << "\r\033[1;34m[BEST fitting score] Processing spectra fit range: [" << fit_range.first << ", " << fit_range.second << "] (GeV/c)...\033[0m (";
                std::cout << "\033[2B";
            int ctr = 1;
            bool first = true; bool centrality_first = true;
            std::pair<double, double> range = pTrange;
//            std::pair<double, double> range = std::make_pair(0.0, 10.0);
                std::cout << "\033[1F";
                std::cout << "\r\033[2K";
                std::cout << "\r\033[1;32mProcessing pT range: [" << range.first << ", " << range.second << "] (GeV/c)...\033[0m (" << ctr << "/" << Parameterization::GetPtRanges().size() << " ranges) ⚙️";
                std::cout << "\033[1B";
                ctr++;
                double low_cut{range.first}; double high_cut{range.second}; // pT range selection
                
                int ref_bin{12}; // <- you can try varying reference class here 5%-> 13 (OLD: 5%->17; 10%->23)
                const int low_ref{params.GetFineAmplitudeBinCutLow(ref_bin)}; /// Obtaining reference bins ---> There are 49 elements (30 corresponds
                const int high_ref{params.GetFineAmplitudeBinCutHigh(0)}; /// Obtaining reference bins:: Here is where the data needs to be pulled from .txt
                TH1D* hpt_ref = (TH1D*)Pt_Nch_hist->ProjectionY("(reference) 0-5%",low_ref,high_ref);
                hpt_ref->Sumw2();
                hpt_ref->SetTitle("0-5%");
                hpt_ref->SetYTitle("1/N_{evt} 1/(2#pi p_{T}) d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-2}");
                params.NormalizeByBinWidth(hpt_ref); //1/(2pi*pT) and 1/(dpT)
                hpt_ref->Scale(1./Nch_dist->Integral(low_ref,high_ref)); //<-- Scales spectra by 1/events
                params.SetErrorYield(hpt_ref, 29); //Second argument is centrality index (the ref is stored in the last .txt file in yield_errors/)
//                params.SetErrorYield(hpt_ref, 36); //Second argument is centrality index (the ref is stored in the last .txt file in yield_errors/)
                    std::ostringstream oss;
                    TGraphErrors* gdN = new TGraphErrors();
                    gdN->SetMarkerStyle(20);
                    gdN->SetMarkerColor(kGreen);
                    gdN->SetLineWidth(2);
                    gdN->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
                    gdN->GetYaxis()->SetTitle("#LTp_{T}#GT/#LTp_{T}#GT_{0-5%}");
                    oss << low_cut << " #leq p_{T} #leq " << high_cut << " (GeV/c)";
                    gdN->SetName(oss.str().c_str());
                    oss.str(""); oss.clear();
                    if (range.first == 0. && range.second == 10.){
                        oss << "STAR (Au+Au #sqrt{s_{NN}} = 200 GeV), " << "p_{T} #leq " << high_cut << " (GeV/c), |#eta| < 0.5";
                    }
                    else {oss << "STAR (Au+Au #sqrt{s_{NN}} = 200 GeV), " << low_cut << " #leq p_{T} #leq " << high_cut << " (GeV/c), |#eta| < 0.5";}
                    gdN->SetTitle(oss.str().c_str());
                    oss.str(""); oss.clear();
                    
                    TGraphErrors* gAN = new TGraphErrors(); //AbsY NormX (CHANGE THIS BELOW! AND FILL HISTO)
                    gAN->SetMarkerStyle(20);
                    gAN->SetMarkerColor(kGreen);
                    gdN->SetLineWidth(0);
                    gdN->SetLineColor(0);
                    gdN->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
                    gdN->GetYaxis()->SetTitle("#LTp_{T}#GT/#LTp_{T}#GT_{0-5%}");
                    oss << low_cut << " #leq p_{T} #leq " << high_cut << " (GeV/c)";
                    gdN->SetName(oss.str().c_str());
                    oss.str(""); oss.clear();
                    if (range.first == 0. && range.second == 10.){
                        oss << "STAR (Au+Au #sqrt{s_{NN}} = 200 GeV), " << "p_{T} #leq " << high_cut << " (GeV/c), |#eta| < 0.5";
                    }
                    else {oss << "STAR (Au+Au #sqrt{s_{NN}} = 200 GeV), " << low_cut << " #leq p_{T} #leq " << high_cut << " (GeV/c), |#eta| < 0.5";}
                    gdN->SetTitle(oss.str().c_str());
                    oss.str(""); oss.clear();
                    
                    oss << low_cut << " #leq p_{T} #leq " << high_cut << " (GeV/c)";
                    TGraphErrors* g = new TGraphErrors();
                    g->SetMarkerStyle(20);
                    g->SetLineWidth(0);
                    g->SetLineColor(0);
                    g->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%} (Reference)");
                    g->GetYaxis()->SetTitle("#LTp_{T}#GT/#LTp_{T}#GT_{0-5%}");
                    g->SetName(oss.str().c_str());
                    
                    TGraphErrors* gdNX = new TGraphErrors();
                    gdNX->SetMarkerStyle(20);
                    gdNX->SetLineWidth(0);
                    gdNX->SetLineColor(0);
                    gdNX->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
                    gdNX->GetYaxis()->SetTitle("#LTp_{T}#GT");
                    gdNX->SetName(oss.str().c_str());
                    
                    TGraphErrors* ggAna = new TGraphErrors();
                    ggAna->SetMarkerStyle(24);
                    ggAna->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT");
                    ggAna->GetYaxis()->SetTitle("#LTp_{T}#GT");
                    ggAna->SetName(oss.str().c_str());
                    
                    TGraphErrors* gg = new TGraphErrors();
                    gg->SetMarkerStyle(24);
                    gg->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT (Reference)");
                    gg->GetYaxis()->SetTitle("#LTp_{T}#GT");
                    gg->SetName(oss.str().c_str());
    
                    TGraphErrors* gNch = new TGraphErrors();
                    gNch->SetMarkerStyle(24);
                    gNch->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
                    gNch->GetYaxis()->SetTitle("Event fraction");
                    gNch->SetName(oss.str().c_str());
                    oss.str(""); oss.clear();
                    
                    if (params.GetDebugSetting()){
                        cout << "Perc 0-5: " << Nch_dist->Integral(low_ref,high_ref) / Nch_dist->Integral(1,Nch_dist->GetNbinsX());
                        printf("\n\t low_ref = %d | high_ref = %d\n",low_ref,high_ref);
                        printf("\n\t Number of FINE bins: %d\n",params.GetnFinePercentileBinning());
                        cout << "\n\t Number of events: " << Nch_dist->Integral(low_ref,high_ref) << endl;
                    }
                    
                    
                    oss << low_cut << " #leq p_{T} #leq " << high_cut << " (GeV/c)" << Form("[%s]", Parameterization::FunctionToString(fitFuncType).c_str());
                    if (params.GetDebugSetting()){ std::cout << "Processing " << Parameterization::FunctionToString(fitFuncType).c_str() << "..." << std::endl; }
                    std::string name = oss.str();
                    if (!gdN) {
                        std::cerr << "[ERROR] gdN is a null pointer before SetName. Skipping this entry.\n";
                    }
                    gdN->SetName(name.c_str());
                    oss.str(""); oss.clear();
                    
                    oss << "dNdEta " << low_cut << " #leq p_{T} #leq " << high_cut << " (GeV/c) " << Form("fit_range_%.2f-%.2f_GeV",fit_range.first, fit_range.second) << Form(" [%s]", Parameterization::FunctionToString(fitFuncType).c_str());
                    TGraphErrors* gCentrality = new TGraphErrors();
                    gCentrality->SetMarkerStyle(24);
                    gCentrality->GetXaxis()->SetTitle("Centrality (%)");
                    gCentrality->GetYaxis()->SetTitle("#LTdN_{ch}/d#eta#GT");
                    gCentrality->SetName(oss.str().c_str());
                    oss.str(""); oss.clear();
                    
                    oss << "<pT> " << low_cut << " #leq p_{T} #leq " << high_cut << " (GeV/c) " << Form("fit_range_%.2f-%.2f_GeV",fit_range.first, fit_range.second) << Form(" [%s]", Parameterization::FunctionToString(fitFuncType).c_str());
                    TGraphErrors* gpTCentrality = new TGraphErrors();
                    gpTCentrality->SetMarkerStyle(24);
                    gpTCentrality->GetXaxis()->SetTitle("Centrality (%)");
                    gpTCentrality->GetYaxis()->SetTitle("#LT p_{T} #GT");
                    gpTCentrality->SetName(oss.str().c_str());
                    oss.str(""); oss.clear();
                    
                    
                    if (params.GetDebugSetting()) {cout << "Fitting for 0-5%...\n";}
                    TString tf1Name = Form("fit_%s", Parameterization::FunctionToString(fitFuncType).c_str());
                    TF1* fit = params.GetFit(hpt_ref, fit_range.first, fit_range.second, fitFuncType, tf1Name, -1);
                    double chi2 = fit->GetChisquare();
                    int ndf     = fit->GetNDF();
                    double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;
                    if (chi2ndf > 1000. || chi2ndf < 1e-9) {
                        if (params.GetDebugSetting()) std::cout << "Poor fit on reference bin. Skipping..." << std::endl;
                    }
                    
                    TString funcStr = Parameterization::FunctionToString(fitFuncType);
                    TString label = Form("D/F 0-5%% [%s]", funcStr.Data());
                    hpt_ref->SetName(Form("(ref) 0-5%% [%s]", funcStr.Data()));
                    
                    if (params.GetDebugSetting()){ std::cout << "Attempting " << std::string(funcStr.Data()) << " fit..." << std::endl; }
                    if (first){
                        spectrafitRangeSubdir->cd();
    //                    TH1* hFitRatio_ref = params.GetDataFitRatio(hpt_ref, fit, "D/F 0-5%");
    //                    params.GetDataFitRatio(hpt_ref, fit, std::string(label.Data()));
    //                    std::cout << "Writing reference spectrum: " << hpt_ref->GetName() << std::endl;
                        hpt_ref->Write();
                    }
                        
                    const double meanpt_ref{params.GetMean(hpt_ref,fit,true,high_cut,low_cut)};
                    const double meannch_ref{params.GetNch(Nch_dist,low_ref,high_ref)};
                    const double dndetaAna_ref{params.GetdNdEta(hpt_ref,fit,true,high_cut,low_cut)};
                    const double error_ref{params.GetErrorMeanpT(hpt_ref,fit,true,high_cut)};
                    
                    if (params.GetDebugSetting()){
                        cout << "\n( Before event & pT normalizations ): " << params.GetdNdEta(hpt_ref,fit,true,0.2,0);
                        cout << '\n' << "<pT>_{0-5%}: " << meanpt_ref;
                        cout << '\n' << "<dn/deta>_{0-5%, Mean Nch}: " << meannch_ref;
                        cout << '\n' << "<dn/deta>_{0-5%, Integrated}: " << dndetaAna_ref;
                    }
                    
                    // -----**-----**-----**-----**-----**--(3) Bin Calculations--**-----**-----**-----**-----**-----**-----**
                    nperc_bins = params.GetnFinePercentileBinning();
                    for (int k = 0; k < nperc_bins; k++) {
                        bool isBad = false;
                        for (const auto& bad_point : Parameterization::Getbadpoints()) {
                            if (k == bad_point){
                                isBad = true;
                                break;
                            }
                        }
                        if (isBad) continue;
                        
                        if (params.GetDebugSetting()){
                            cout << "\nparams.GetFineAmplitudeCutLow(k)*100: " << params.GetFineAmplitudeCutLow(k)*100;
                            cout << "\nparams.GetFineAmplitudeCutHigh(k)*100: " << params.GetFineAmplitudeCutHigh(k)*100;
                            cout << "\nparams.GetFineAmplitudeBinCutLow(k): " << params.GetFineAmplitudeBinCutLow(k);
                            cout << "\nparams.GetFineAmplitudeBinCutHigh(k): " << params.GetFineAmplitudeBinCutHigh(k) << endl;
                        }
                        
                        TH1D* hpt = nullptr;
                        
                        if (params.GetCentralitySetting() == "I"){ //Inclusive Centrality Analysis
                            hpt = Pt_Nch_hist->ProjectionY(Form("%.3f-%.3f%% [%s]", params.GetFineAmplitudeCutLow(0)*100,
                                                                params.GetFineAmplitudeCutHigh(k)*100,
                                                                funcStr.Data()),
                                                           params.GetFineAmplitudeBinCutLow(k),
                                                           params.GetFineAmplitudeBinCutHigh(0));
                            hpt->SetTitle(Form("%.0f-%.4f%%", params.GetFineAmplitudeCutLow(0)*100,
                                               params.GetFineAmplitudeCutHigh(k)*100));
                        }
                        else{                                      //Exclusive Centrality Analysis
                            hpt = Pt_Nch_hist->ProjectionY(Form("%.3f-%.3f%% [%s]", params.GetFineAmplitudeCutLow(k)*100,
                                                                params.GetFineAmplitudeCutHigh(k)*100,
                                                                funcStr.Data()),
                                                           params.GetFineAmplitudeBinCutLow(k),
                                                           params.GetFineAmplitudeBinCutHigh(k));
                            hpt->SetTitle(Form("%.3f-%.3f%%",
                                               params.GetFineAmplitudeCutLow(k)*100,
                                               params.GetFineAmplitudeCutHigh(k)*100));
                        }
                        
                        hpt->Sumw2();
                        hpt->SetYTitle("1/N_{evt} 1/(2#pi p_{T}) d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-2}"); // Rutik/Rene's version (differential yield)
                        //            hpt->SetYTitle("1/N_{evt} d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-1}"); // Omar's version (spectra)
                        params.NormalizeByBinWidth(hpt);
                        if (params.GetCentralitySetting() == "I"){ //Inclusive Centrality Analysis
                            hpt->Scale(1./Nch_dist->Integral(params.GetFineAmplitudeBinCutLow(k),
                                                             params.GetFineAmplitudeBinCutHigh(0))); //Inclusive
                        }
                        else{                                      //Exclusive Centrality Analysis
                            hpt->Scale(1./Nch_dist->Integral(params.GetFineAmplitudeBinCutLow(k),
                                                             params.GetFineAmplitudeBinCutHigh(k))); //Exclusive
                        }
                        params.SetErrorYield(hpt, k);
                        
                        if (params.GetDebugSetting()){
                            std::cout << "\n Setting error in pT spectra: " << std::endl;
                            for (int bin = 1; bin < 10; bin++) {
                                std::cout << " [" << bin << "] error: " << hpt->GetBinError(bin);
                            }
                            std::cout << std::endl;
                        }
                        
                        bool ifpt_cut = true; //Enforcing 10.0 GeV and 200 MeV cutoff
                        tf1Name = Form("hfit_%s", Parameterization::FunctionToString(fitFuncType).c_str());
                        if (params.GetDebugSetting()) std::cout << "Going into fit" << std::endl;
                        fit = params.GetFit(hpt, fit_range.first, fit_range.second, fitFuncType, tf1Name);
                        //                    fit = params.GetFit(hpt, fit_range.first, fit_range.second, fitFuncType);
                        if (params.GetDebugSetting()) std::cout << "Fit completed" << std::endl;
                        
                        chi2 = fit->GetChisquare();
                        ndf     = fit->GetNDF();
                        chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;
                        if (k < 20 && (chi2ndf > 1000. || chi2ndf < 1e-9)) continue;
                        
                        if (params.GetDebugSetting()) {
                            std::cout << "Writing hpt";
                            cout << "\n Index: " << k;
                            cout << "\n first: " << first << std::endl;
                        }
                        
                        if (first){
                            spectrafitRangeSubdir->cd();
                            label = Form("D/F %.4f-%.4f%% [%s]",
                                         params.GetFineAmplitudeCutLow(k)*100,
                                         params.GetFineAmplitudeCutHigh(k)*100,
                                         funcStr.Data());
    //                        params.GetDataFitRatio(hpt, fit, std::string(label.Data()));
                            hpt->SetMarkerStyle(4); hpt->SetMarkerColor(kBlack); hpt->SetMarkerSize(1.2);
                            hpt->Write();
                        }
                        
                        
                        double mean_nch;
                        if (params.GetCentralitySetting() == "I"){
                            mean_nch = params.GetNch(Nch_dist,params.GetFineAmplitudeBinCutLow(0),params.GetFineAmplitudeBinCutHigh(k));
                        }
                        else { mean_nch = params.GetNch(Nch_dist,params.GetFineAmplitudeBinCutLow(k),params.GetFineAmplitudeBinCutHigh(k)); }

                        if (params.GetDebugSetting()) std::cout << "Calculating mean pT, etc." << std::endl;
                        double mean_pt{params.GetMean(hpt,fit,ifpt_cut,high_cut,low_cut)};
                        double dndetaAna{params.GetdNdEta(hpt, fit, ifpt_cut, high_cut,low_cut)};
                        double error_pt{params.GetErrorMeanpT(hpt,fit,true,high_cut)};
//                        if (error_pt > 0.0025 || error_ref > 0.0025) { if (params.GetDebugSetting()) {std::cout << "SKIPPING from pT error";} continue; }
                        
        
                        
                        double error_dndetaAna{params.GetErrordNdEta(hpt, fit, true, 0.2)};
                        double error_dndetaAna_ref{params.GetErrordNdEta(hpt_ref, fit, true, 0.2)};
    //                    if (error_dndetaAna > 50 || error_dndetaAna_ref > 50) {std::cout << "SKIPPING from dNdeta error"; continue;}
                        double error_nch{params.GetRelativeError(k)};
                        double error_nch_ref{params.GetRelativeError(ref_bin)};
                        
                        double error_ratio{(mean_pt / meanpt_ref) * sqrt(pow(( (error_pt / mean_pt) - (error_ref / meanpt_ref) ), 2) )};
                        // ^Error progation-> https://en.wikipedia.org/wiki/Propagation_of_uncertainty
                        //                    https://arxiv.org/pdf/1609.04150#page24
                        double error_nch_ratio{(mean_nch / meannch_ref) * sqrt(pow( (error_nch / mean_nch), 2.0 ) +  pow( (error_nch_ref / meannch_ref), 2.0 ) ) };
                        // ^Error progation-> https://en.wikipedia.org/wiki/Propagation_of_uncertainty
                        //                    https://arxiv.org/pdf/1609.04150#page24
                        double error_dndeta_ratio{( dndetaAna / dndetaAna_ref) * sqrt(pow(( (error_dndetaAna / dndetaAna) - (error_dndetaAna_ref / dndetaAna_ref) ), 2) )};
                        // ^Error progation-> https://en.wikipedia.org/wiki/Propagation_of_uncertainty
                        //                    https://arxiv.org/pdf/1609.04150#page24
                        
                        if (params.GetDebugSetting()) {
                            cout << "\n Index: " << k;
                            cout << "\n Number of events: " << Nch_dist->Integral(params.GetFineAmplitudeBinCutLow(k),params.GetFineAmplitudeBinCutHigh(0));
                            cout << "\n params.GetFineAmplitudeBinCutHigh(k): " << params.GetFineAmplitudeBinCutHigh(0);
                            cout << "\n params.GetFineAmplitudeBinCutLow(k): " << params.GetFineAmplitudeBinCutLow(k);
                            cout << "\n dNdeta / dNdeta_0: " << dndetaAna / dndetaAna_ref;
                            cout << "\n mean_pt/meanpt_ref: " << mean_pt/meanpt_ref;
                            cout << "\n chi^2/ndf: " << mean_pt/meanpt_ref;
                            std::cout << "Setting points." << std::endl;
                            std::cout << "dndetaAna: " << dndetaAna << std::endl;
                            std::cout << "dndetaAna_ref: " << dndetaAna_ref << std::endl;
                            std::cout << "mean_pt: " << mean_pt << std::endl;
                            std::cout << "meanpt_ref: " << meanpt_ref << std::endl;
                            std::cout << "(mean_pt/meanpt_ref)" << (mean_pt/meanpt_ref);
                            std::cout << "(dndetaAna / dndetaAna_ref)" << (dndetaAna / dndetaAna_ref);
                            std::cout << "error_dndeta_ratio" << error_dndeta_ratio;
                            std::cout << "error_ratio" << error_ratio;
                        }
                        
//                        if (k > 1 && k < 20){ //Checks if points are monotonic within errors
//                            double prev_x, prev_y;
//                            gdN->GetPoint(k - 1, prev_x, prev_y);
//                            double prev_x_err = gdN->GetErrorX(k - 1);
//                            double prev_y_err = gdN->GetErrorY(k - 1);
//                            if (!params.IsMonotonicWithinError(prev_x, prev_y,
//                                                        dndetaAna / dndetaAna_ref, mean_pt/meanpt_ref,
//                                                        prev_x_err, prev_y_err,
//                                                        error_dndeta_ratio, error_ratio)) continue;
//                        }
                        
                        if (params.GetDebugSetting()) std::cout << "Setting points gdN." << std::endl;
                        if (params.GetDebugSetting()) std::cout << "Setting points k = " << k << "with x = " << dndetaAna / dndetaAna_ref << std::endl;
                        gdN->SetPoint(k, dndetaAna / dndetaAna_ref, mean_pt/meanpt_ref);
                        
                        gdN->SetPointError(k, error_dndeta_ratio, error_ratio);
                        gdN->GetXaxis()->SetRangeUser(0.8, 1.25);
                        gdN->GetYaxis()->SetRangeUser(0.99, 1.025);
                        
                        
                        if (params.GetDebugSetting()) std::cout << "Setting points gdNX." << std::endl;
                        gdNX->SetPoint(k, dndetaAna / dndetaAna_ref, mean_pt);
                        gdNX->SetPointError(k, error_dndeta_ratio, error_ratio);
                        gdNX->GetXaxis()->SetRangeUser(0.8, 1.25);
                        gdNX->GetYaxis()->SetRangeUser(0.5, 0.7);
                        
                        
                        if (params.GetDebugSetting()) std::cout << "Setting points g." << std::endl;
                        g->SetPoint(k, mean_nch / meannch_ref, mean_pt / meanpt_ref);
                        g->SetPointError(k, error_nch_ratio, error_ratio);
                        g->GetXaxis()->SetRangeUser(0.8, 1.25);
                        g->GetYaxis()->SetRangeUser(0.98, 1.025);
                        
                        
                        if (params.GetDebugSetting()) std::cout << "Setting points gg." << std::endl;
                        gg->SetPoint(k, mean_nch, mean_pt);
                        gg->SetPointError(k, error_nch, error_pt);
                        
                        if (params.GetDebugSetting()) std::cout << "Setting points ggAna." << std::endl;
                        ggAna->SetPoint(k, dndetaAna, mean_pt);
                        ggAna->SetPointError(k, error_dndetaAna, error_pt);
                        
                        if (params.GetDebugSetting()) std::cout << "Setting points gAna." << std::endl;
                        gAna->SetPoint(k, mean_nch, dndetaAna );
                        gAna->SetPointError(k, error_nch, error_dndetaAna);
                        
                        if (params.GetDebugSetting()) std::cout << "Setting points gAnaRef." << std::endl;
                        gAnaRef->SetPoint(k, mean_nch/meannch_ref, dndetaAna / dndetaAna_ref);
                        gAnaRef->SetPointError(k, error_nch_ratio, error_dndeta_ratio);
                        gAnaRef->GetXaxis()->SetRangeUser(0.2, 1.25);
                        
                        if (params.GetDebugSetting()) std::cout << "Setting points gAnaDRef." << std::endl;
                        gAnaDRef->SetPoint(k, mean_nch, dndetaAna / mean_nch);
                        gAnaDRef->SetPointError(k, error_nch_ratio, 0);
                        gAnaDRef->GetXaxis()->SetRangeUser(100, 800);
                        
                        if (params.GetDebugSetting()) std::cout << "Setting points gAnaDNRef." << std::endl;
                        gAnaDNRef->SetPoint(k, mean_nch, (dndetaAna / dndetaAna_ref) / (mean_nch / meannch_ref));
                        gAnaDNRef->SetPointError(k, 0, 0);
                        gAnaDNRef->GetXaxis()->SetRangeUser(100, 800);
                        
                        if(centrality_first){
                            if (params.GetDebugSetting()) std::cout << "Setting points Centrality." << std::endl;
                            double cent{params.GetCentrality(k, params, false)};
                            gCentrality->SetPoint(k, cent, dndetaAna);
                            gCentrality->SetPointError(k, 0, error_dndetaAna);
                            gpTCentrality->SetPoint(k, cent, mean_pt);
                            gpTCentrality->SetPointError(k, 0, error_pt);
                        }
                        
                        if (params.GetDebugSetting()) std::cout << "Finished setting points" << std::endl;
                        
                    }
                    
                    if (params.GetDebugSetting()) {cout << "\n#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#\n" << endl;}
                    
                    if (centrality_first ){
                        Bestchi2CentralityDir->cd();
                        gCentrality->Write();
                        gpTCentrality->Write();
                        gAnaRef->Write();
                        delete gCentrality;
                    }
                    centrality_first = false;
                    
                    fitRangeSubdir->cd();
                    
                    TGraphErrors* gdN_clone_out = nullptr;;
                    double sos_chi2ndf = 0.0; double sos_fit_score = 0.0; double fit_score_tmp = 0.0;
                    for (const auto& sos_range : Parameterization::GetsosRanges()){
                        if (gdN && gdN->GetN() > 0 && hasPointsInRange(gdN, sos_range.first, sos_range.second)) {
                            std::cout << "\r\033[2K";
                            std::cout << "\r\033[1;31mProcessing " << Parameterization::GetsosRanges().size()
                            << " speed of sound fit ranges: ["
                            << Parameterization::GetsosRanges().front().first << ", " << Parameterization::GetsosRanges().front().second
                            << "] - ["
                            << Parameterization::GetsosRanges().back().first << ", " << Parameterization::GetsosRanges().back().second
                            << "] (GeV/c)...\033[0m⚙️" << std::flush;
                            
                            if (params.GetDebugSetting()) std::cout << "Performing sound fit" << std::endl;
                            TGraphErrors* gdN_clone = params.GetSoundFit(gdN, sos_range.first, sos_range.second, lOutPt_Norm_Int, low_cut, high_cut);
                            
                            TF1* fitFunc = gdN_clone->GetFunction("PowerLaw_fit");
                            double cs2_err = fitFunc->GetParError(0);
                            
                            if (fitFunc && gdN_clone) {
                                double chi2 = fitFunc->GetChisquare();
                                double ndf  = fitFunc->GetNDF();
                                sos_chi2ndf = (ndf > 0) ? chi2 / ndf : -1;
                                sos_fit_score = params.ComputeFitScore(sos_chi2ndf);
                                if (sos_fit_score > fit_score_tmp && cs2_err < 0.04){
                                    if (params.GetDebugSetting()) std::cout << "[BestChi2 speed of sound] cloning gdN" << std::endl;
                                    if (gdN_clone_out) delete gdN_clone_out;
                                    if (!gdN_clone) {if (params.GetDebugSetting()) std::cout << "[BestChi2 speed of sound] NO gdN!" << std::endl;}
                                    if (params.GetDebugSetting()) std::cout << "[BestChi2 speed of sound] PLACEHOLDER" << std::endl;
                                    gdN_clone_out = (TGraphErrors*)gdN_clone->Clone();
                                    if (params.GetDebugSetting()) std::cout << "[BestChi2 speed of sound] PLACEHOLDER" << std::endl;
//                                    gdN_clone_out = gdN_clone;
                                }
                            }
                            delete gdN_clone;
                        }
                    } //sos fit range loop end
                    
                    if (params.GetDebugSetting()) std::cout << "[BestChi2 speed of sound] writing cloned gdN" << std::endl;
                    if (gdN_clone_out) gdN_clone_out->Write();
                    if (gdN_clone_out) delete gdN_clone_out;
    
                    if (params.GetDebugSetting()) std::cout << "Writing gdN" << std::endl;
                    gdN->Write();
                    delete gdN;
                    
                    if (params.GetDebugSetting()) std::cout << "Writing gdNX" << std::endl;
                    fitRangeSubdirX->cd();
                    gdNX->Write();
                    delete gdNX;
                    
                    if (params.GetDebugSetting()) std::cout << "Writing ggAna" << std::endl;
                    fitRangeSubdirAbs->cd();
                    ggAna->Write();
                    delete ggAna;
                    
                    if (params.GetDebugSetting()) std::cout << "End of fit function loop" << std::endl;
                
                first = false;

        // -----**-----**-----**-----**-----**--(5)--**-----**-----**-----**-----**-----**-----**
    
        gROOT->SetBatch(kFALSE);
}


// -----------------------------------------------------

void Parameterization::MomentAna(TH2D* Pt_Nch_hist, TH2D* Pt_Mult_hist, TH2D* Nch_Mult_hist, TH1D* Nch_dist, Parameterization params) {
        gROOT->SetBatch(kTRUE); std::string output_label = "./Momentgraphs_output.root";
        TH1::AddDirectory(kFALSE);

        int nperc_bins{0};
    
//        std::vector<Function> fitFunctions = {Function::LevyTsallis, Function::Hagedorn, Function::Tsallis};
        std::vector<Function> fitFunctions = {Function::LevyTsallis,};
        
        // -----**-----**-----**-----**-----**--(1)--**-----**-----**-----**-----**-----**-----**
        
        MomentHistogram* h = new MomentHistogram();
        TFile* fOut = new TFile(output_label.c_str(),"recreate");
        TList* lOutGen = new TList(); // TList for general elements (QA, etc)
        lOutGen->SetOwner(true);
        
        lOutGen->Add(Pt_Nch_hist);
        lOutGen->Add(Pt_Mult_hist);
        lOutGen->Add(Nch_Mult_hist);
        lOutGen->Add(Nch_dist);
        
        TDirectory* ptStudyDir = fOut->mkdir("[pT] study");
        TDirectory* ptStudyDir_Dist = ptStudyDir->mkdir("Distributions");
    
        TDirectory* ptStudyDir_Cumulants = ptStudyDir->mkdir("Cumulants");
        TDirectory* ptStudyDir_Cumulants_Abs = ptStudyDir_Cumulants->mkdir("Absolute quantities");
        TDirectory* ptStudyDir_Cumulants_Norm = ptStudyDir_Cumulants->mkdir("Normalized (to ref class) quantities");
    
        TDirectory* ptStudyDir_Moments = ptStudyDir->mkdir("Moments");
        TDirectory* ptStudyDir_Moments_Abs = ptStudyDir_Moments->mkdir("Absolute quantities");
        TDirectory* ptStudyDir_Moments_Norm = ptStudyDir_Moments->mkdir("Normalized (to ref class) quantities");
    
        TDirectory* NchStudyDir = fOut->mkdir("Nch study");
        TDirectory* NchStudyDir_Dist = NchStudyDir->mkdir("Distributions");
        
        TDirectory* NchStudyDir_Cumulants = NchStudyDir->mkdir("Cumulants");
        TDirectory* NchStudyDir_Cumulants_Abs = NchStudyDir_Cumulants->mkdir("Absolute quantities");
        TDirectory* NchStudyDir_Cumulants_Norm = NchStudyDir_Cumulants->mkdir("Normalized (to ref class) quantities");
    
        TDirectory* NchStudyDir_Moments = NchStudyDir->mkdir("Moments");
        TDirectory* NchStudyDir_Moments_Abs = NchStudyDir_Moments->mkdir("Absolute quantities");
        TDirectory* NchStudyDir_Moments_Norm = NchStudyDir_Moments->mkdir("Normalized (to ref class) quantities");
        
        TDirectory* spectraStudyDir = fOut->mkdir("spectra study");
        TDirectory* generalDir = fOut->mkdir("general");
        
        
        std::cout << std::endl; std::cout << std::endl; std::cout << std::endl;
        double best_fitting_score = 0.0;
        Function BestChi2fitfunctype;
        std::pair<double, double> BestChi2fitrange;
        std::pair<double, double> pTrange;
    
        for (const auto& fit_range : Parameterization::GetfitRanges()) { //**pT fit range loop start**
                double fitting_score = 0.0;
                std::ostringstream subdir_name;
                subdir_name << Form("fit_range_%.2f-%.2f_GeV", fit_range.first, fit_range.second);
//                TDirectory* spectrafitRangeSubdir = spectraStudyDir->mkdir(subdir_name.str().c_str());
//            int ctr = 1;
//            bool first = true; bool centrality_first = true;
            for (const auto& range : Parameterization::GetPtRanges()) {
//                ctr++;
//                double low_cut{range.first}; double high_cut{range.second}; // pT range selection
                
                int ref_bin{12}; // <- you can try varying reference class here 5%->13 (OLD: 5%->17; 10%->23)
                const int low_ref{params.GetFineAmplitudeBinCutLow(ref_bin)}; /// Obtaining reference bins ---> There are 49 elements (30 corresponds
                const int high_ref{params.GetFineAmplitudeBinCutHigh(0)}; /// Obtaining reference bins:: Here is where the data needs to be pulled from .txt
                TH1D* hpt_ref = (TH1D*)Pt_Nch_hist->ProjectionY("(reference) 0-5%",low_ref,high_ref);
                hpt_ref->Sumw2();
                hpt_ref->SetTitle("0-5%");
                hpt_ref->SetYTitle("1/N_{evt} 1/(2#pi p_{T}) d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-2}");
                params.NormalizeByBinWidth(hpt_ref); //1/(2pi*pT) and 1/(dpT)
                hpt_ref->Scale(1./Nch_dist->Integral(low_ref,high_ref)); //<-- Scales spectra by 1/events
                params.SetErrorYield(hpt_ref, 29); //Second argument is centrality index (the ref is stored in the last .txt file in yield_errors/)
//                params.SetErrorYield(hpt_ref, 36); //Second argument is centrality index (the ref is stored in the last .txt file in yield_errors/)
                for (const auto& fitFuncType : fitFunctions) { //**fit function species loop start**
                    
                    if (params.GetDebugSetting()) {cout << "Fitting for 0-5%...\n";}
                    TString tf1Name = Form("fit_%s", Parameterization::FunctionToString(fitFuncType).c_str());
                    TF1* fit = params.GetFit(hpt_ref, fit_range.first, fit_range.second, fitFuncType, tf1Name);
                    double chi2 = fit->GetChisquare();
                    int ndf     = fit->GetNDF();
                    double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;
                    if (chi2ndf > 1000. || chi2ndf < 1e-9) {
                        if (params.GetDebugSetting()) std::cout << "Poor fit on reference bin. Skipping..." << std::endl;
                        continue;
                    }
                    
                    fitting_score += params.ComputeFitScore(chi2ndf);
                    TString funcStr = Parameterization::FunctionToString(fitFuncType);
                    
                    // -----**-----**-----**-----**-----**--(3) Bin Calculations--**-----**-----**-----**-----**-----**-----**
                    nperc_bins = params.GetnFinePercentileBinning();
                    double iter = 0.0;
                    for (int k = 0; k < nperc_bins; k++) {
                        iter += 1.;
                        bool isBad = false;
                        for (const auto& bad_point : Parameterization::Getbadpoints()) {
                            if (k == bad_point){
                                isBad = true;
                                break;
                            }
                        }
                        if (isBad) continue;
                        TH1D* hpt = nullptr;
                        
                        if (params.GetCentralitySetting() == "I"){ //Inclusive Centrality Analysis
                            hpt = Pt_Nch_hist->ProjectionY(Form("%.3f-%.3f%% [%s]", params.GetFineAmplitudeCutLow(0)*100,
                                                                params.GetFineAmplitudeCutHigh(k)*100,
                                                                funcStr.Data()),
                                                           params.GetFineAmplitudeBinCutLow(k),
                                                           params.GetFineAmplitudeBinCutHigh(0));
                        }
                        else{                                      //Exclusive Centrality Analysis
                            hpt = Pt_Nch_hist->ProjectionY(Form("%.3f-%.3f%% [%s]", params.GetFineAmplitudeCutLow(k)*100,
                                                                params.GetFineAmplitudeCutHigh(k)*100,
                                                                funcStr.Data()),
                                                           params.GetFineAmplitudeBinCutLow(k),
                                                           params.GetFineAmplitudeBinCutHigh(k));
                        }
                        
                        hpt->Sumw2();
                        hpt->SetYTitle("1/N_{evt} 1/(2#pi p_{T}) d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-2}"); // Rutik/Rene's version (differential yield)
                        params.NormalizeByBinWidth(hpt);
                        if (params.GetCentralitySetting() == "I"){ //Inclusive Centrality Analysis
                            hpt->Scale(1./Nch_dist->Integral(params.GetFineAmplitudeBinCutLow(k),
                                                             params.GetFineAmplitudeBinCutHigh(0))); //Inclusive
                        }
                        else{                                      //Exclusive Centrality Analysis
                            hpt->Scale(1./Nch_dist->Integral(params.GetFineAmplitudeBinCutLow(k),
                                                             params.GetFineAmplitudeBinCutHigh(k))); //Exclusive
                        }
                        params.SetErrorYield(hpt, k);
                        tf1Name = Form("hfit_%s", Parameterization::FunctionToString(fitFuncType).c_str());
                        if (params.GetDebugSetting()) std::cout << "Going into fit" << std::endl;
                        fit = params.GetFit(hpt, fit_range.first, fit_range.second, fitFuncType, tf1Name);
                        if (params.GetDebugSetting()) std::cout << "Fit completed" << std::endl;
                        
                        chi2 = fit->GetChisquare();
                        ndf     = fit->GetNDF();
                        chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;
                        if (k < 20 && (chi2ndf > 1000. || chi2ndf < 1e-9)) continue;
                        fitting_score += params.ComputeFitScore(chi2ndf);
                        
                    }
                    
                    fitting_score /= (iter + 1);
                    if (fitting_score > best_fitting_score){
                        if (params.GetDebugSetting()) std::cout << "\n Fitting score: " << fitting_score << std::endl;
                        best_fitting_score = fitting_score;
                        BestChi2fitfunctype = fitFuncType;
                        BestChi2fitrange = fit_range;
                        pTrange = range;
                    }
                } //fit function loop end
//                first = false;
            } //pT fit range loop end
            
            if (params.GetDebugSetting()) std::cout << "the best fitting score was " << best_fitting_score << std::endl;
        } //spectra fit range end
    
            Function fitFuncType = BestChi2fitfunctype;
            std::pair<double, double> fit_range = BestChi2fitrange;

                    std::ostringstream subdir_name;
                    std::cout << "Processing fit range: " << fit_range.first << ", " << fit_range.second << std::endl;
                    subdir_name << Form("fit_range_%.2f-%.2f_GeV", fit_range.first, fit_range.second);
                    TDirectory* spectrafitRangeSubdir = spectraStudyDir->mkdir(subdir_name.str().c_str());
                    std::cout << "\033[2F";
                    std::cout << "\r\033[2K";
                    std::cout << "\r\033[1;34m[Moment Analysis + BEST fitting score] Processing spectra fit range: [" << fit_range.first << ", " << fit_range.second << "] (GeV/c)...\033[0m (";
                    std::cout << "\033[2B";
                int ctr = 1;
                bool first = true; //bool centrality_first = true;
                std::pair<double, double> range = pTrange;
                    std::cout << "\033[1F";
                    std::cout << "\r\033[2K";
                    std::cout << "\r\033[1;32mProcessing pT range: [" << range.first << ", " << range.second << "] (GeV/c)...\033[0m (" << ctr << "/" << Parameterization::GetPtRanges().size() << " ranges) ⚙️";
                    std::cout << "\033[1B";
                    ctr++;
                    double low_cut{range.first}; double high_cut{range.second}; // pT range selection
                    int ref_bin{12}; // <- you can try varying reference class here 5%->17; 10%->23
                    const int low_ref{params.GetFineAmplitudeBinCutLow(ref_bin)}; /// Obtaining reference bins ---> There are 49 elements (30 corresponds
                    const int high_ref{params.GetFineAmplitudeBinCutHigh(0)}; /// Obtaining reference bins:: Here is where the data needs to be pulled from .txt
    
                
//                    TH1D* hpt_ref = (TH1D*)Pt_Nch_hist->ProjectionY("(reference) 0-5%",low_ref,high_ref);
//                    TH1D* hMpt_ref = (TH1D*)Pt_Mult_hist->ProjectionY("[p_{T}] (reference) 0-5%",low_ref,high_ref);
//                    TH1D* hNch_ref = (TH1D*)Nch_Mult_hist->ProjectionY("Nch (reference) 0-5%",low_ref,high_ref);
                    TH1D* hpt_ref  = (TH1D*)Pt_Nch_hist->ProjectionY("(reference) 0-5%", low_ref, high_ref);
                    hpt_ref->SetDirectory(nullptr);

                    TH1D* hMpt_ref = (TH1D*)Pt_Mult_hist->ProjectionY("[p_{T}] (reference) 0-5%", low_ref, high_ref);
                    hMpt_ref->SetDirectory(nullptr);

                    TH1D* hNch_ref = (TH1D*)Nch_Mult_hist->ProjectionY("Nch (reference) 0-5%", low_ref, high_ref);
                    hNch_ref->SetDirectory(nullptr);

    
                    hpt_ref->Sumw2();
                    hpt_ref->SetTitle("0-5%");
                    hpt_ref->SetYTitle("1/N_{evt} 1/(2#pi p_{T}) d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-2}");
                    hMpt_ref->Sumw2();
                    hMpt_ref->SetTitle("0-5%");
                    hMpt_ref->SetYTitle("[p_{T}]");
                    hNch_ref->Sumw2();
                    hNch_ref->SetTitle("0-5%");
                    hNch_ref->SetYTitle("N_{ch}");
                    
                    params.NormalizeByBinWidth(hpt_ref);
                    hpt_ref->Scale(1./Nch_dist->Integral(low_ref,high_ref)); 
                    hMpt_ref->Scale(1./Nch_dist->Integral(low_ref,high_ref));
                    hNch_ref->Scale(1./Nch_dist->Integral(low_ref,high_ref));
//                    params.SetErrorYield(hpt_ref, 36); //Second argument is centrality index (the ref is stored in the last .txt file in yield_errors/)
                    params.SetErrorYield(hpt_ref, 29); //Second argument is centrality index (the ref is stored in the last .txt file in yield_errors/)
                    
                    if (params.GetDebugSetting()){ std::cout << "Calculating cumulants and moments" << std::endl; }
                    h->pTC1_ref = params.GetCumulant(hMpt_ref,1);
                    h->pTC2_ref = params.GetCumulant(hMpt_ref,2);
                    h->pTC3_ref = params.GetCumulant(hMpt_ref,3);
                    h->pTC4_ref = params.GetCumulant(hMpt_ref,4);
                    h->pTM_ref   = ( h->pTC1_ref );
                    h->pTVar_ref = ( h->pTC2_ref );
                    h->pTS_ref   = ( h->pTC3_ref ) / std::pow( h->pTC2_ref, 3./2. );
                    h->pTkap_ref = ( h->pTC4_ref ) / std::pow( h->pTC2_ref, 2. );
    
                    h->NchC1_ref = params.GetCumulant(hNch_ref,1);
                    h->NchC2_ref = params.GetCumulant(hNch_ref,2);
                    h->NchC3_ref = params.GetCumulant(hNch_ref,3);
                    h->NchC4_ref = params.GetCumulant(hNch_ref,4);
                    h->NchM_ref   = ( h->NchC1_ref );
                    h->NchVar_ref = ( h->NchC2_ref );
                    h->NchS_ref   = ( h->NchC3_ref ) / std::pow( h->NchC2_ref, 3./2. );
                    h->Nchkap_ref = ( h->NchC4_ref ) / std::pow( h->NchC2_ref, 2. );
    
                    h->pTC1_ref_error = params.GetCumulantError(hMpt_ref, 1);
                    h->pTC2_ref_error = params.GetCumulantError(hMpt_ref, 2);
                    h->pTC3_ref_error = params.GetCumulantError(hMpt_ref, 3);
                    h->pTC4_ref_error = params.GetCumulantError(hMpt_ref, 4);
                    h->pTM_ref_error   = h->pTC1_ref_error;
                    h->pTVar_ref_error = h->pTC2_ref_error;
                    h->pTS_ref_error   = std::sqrt(
                        std::pow(h->pTC3_ref_error / std::pow(h->pTC2_ref, 1.5), 2) +
                        std::pow(1.5 * h->pTC3_ref * h->pTC2_ref_error / std::pow(h->pTC2_ref, 2.5), 2)
                    );
                    h->pTkap_ref_error = std::sqrt(
                        std::pow(h->pTC4_ref_error / std::pow(h->pTC2_ref, 2), 2) +
                        std::pow(2 * h->pTC4_ref * h->pTC2_ref_error / std::pow(h->pTC2_ref, 3), 2)
                    );

                    h->NchC1_ref_error = params.GetCumulantError(hNch_ref, 1);
                    h->NchC2_ref_error = params.GetCumulantError(hNch_ref, 2);
                    h->NchC3_ref_error = params.GetCumulantError(hNch_ref, 3);
                    h->NchC4_ref_error = params.GetCumulantError(hNch_ref, 4);
                    h->NchM_ref_error   = h->NchC1_ref_error;
                    h->NchVar_ref_error = h->NchC2_ref_error;
                    h->NchS_ref_error   = std::sqrt(
                        std::pow(h->NchC3_ref_error / std::pow(h->NchC2_ref, 1.5), 2) +
                        std::pow(1.5 * h->NchC3_ref * h->NchC2_ref_error / std::pow(h->NchC2_ref, 2.5), 2)
                    );
                    h->Nchkap_ref_error = std::sqrt(
                        std::pow(h->NchC4_ref_error / std::pow(h->NchC2_ref, 2), 2) +
                        std::pow(2 * h->NchC4_ref * h->NchC2_ref_error / std::pow(h->NchC2_ref, 3), 2)
                    );

                if (params.GetDebugSetting()){ std::cout << "Finished calculating cumulants and moments" << std::endl; }
                        if (params.GetDebugSetting()){
                            cout << "Perc 0-5: " << Nch_dist->Integral(low_ref,high_ref) / Nch_dist->Integral(1,Nch_dist->GetNbinsX());
                            printf("\n\t low_ref = %d | high_ref = %d\n",low_ref,high_ref);
                            printf("\n\t Number of FINE bins: %d\n",params.GetnFinePercentileBinning());
                            cout << "\n\t Number of events: " << Nch_dist->Integral(low_ref,high_ref) << endl;
                            std::cout << "Processing " << Parameterization::FunctionToString(fitFuncType).c_str() << "..." << std::endl;
                            std::cout << "Fitting for 0-5%...\n";
                        }
    
                        TString tf1Name = Form("fit_%s", Parameterization::FunctionToString(fitFuncType).c_str());
                        TF1* fit = params.GetFit(hpt_ref, fit_range.first, fit_range.second, fitFuncType, tf1Name);
                        double chi2 = fit->GetChisquare();
                        int ndf     = fit->GetNDF();
                        double chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;
                        if (chi2ndf > 1000. || chi2ndf < 1e-9) {
                            if (params.GetDebugSetting()) std::cout << "Poor fit on reference bin. Skipping..." << std::endl;
                        }
                        const double dndetaAna_ref{params.GetdNdEta(hpt_ref,fit,true,high_cut,low_cut)};
    
                        TString funcStr = Parameterization::FunctionToString(fitFuncType);
                        TString label = Form("D/F 0-5%% [%s]", funcStr.Data());
                        hpt_ref->SetName(Form("(ref) 0-5%% [%s]", funcStr.Data()));
                        
                        if (params.GetDebugSetting()){ std::cout << "Attempting " << std::string(funcStr.Data()) << " fit..." << std::endl; }
                        if (first){
                            spectrafitRangeSubdir->cd();
                            hpt_ref->Write();
                            ptStudyDir_Dist->cd();
                            hMpt_ref->Write();
                            NchStudyDir_Dist->cd();
                            hNch_ref->Write();
                        }
                        
                        // -----**-----**-----**-----**-----**--(3) Bin Calculations--**-----**-----**-----**-----**-----**-----**
                        nperc_bins = params.GetnFinePercentileBinning();
                        for (int k = 0; k < nperc_bins; k++) {
                            bool isBad = false;
                            for (const auto& bad_point : Parameterization::Getbadpoints()) {
                                if (k == bad_point){
                                    isBad = true;
                                    break;
                                }
                            }
                            if (isBad) continue;
                            
                            if (params.GetDebugSetting()){
                                cout << "\nparams.GetFineAmplitudeCutLow(k)*100: " << params.GetFineAmplitudeCutLow(k)*100;
                                cout << "\nparams.GetFineAmplitudeCutHigh(k)*100: " << params.GetFineAmplitudeCutHigh(k)*100;
                                cout << "\nparams.GetFineAmplitudeBinCutLow(k): " << params.GetFineAmplitudeBinCutLow(k);
                                cout << "\nparams.GetFineAmplitudeBinCutHigh(k): " << params.GetFineAmplitudeBinCutHigh(k) << endl;
                            }
                            
                            TH1D* hpt = nullptr;
                            TH1D* hMpt = nullptr;
                            TH1D* hNch = nullptr;
                            
                            if (params.GetCentralitySetting() == "I"){ //Inclusive Centrality Analysis
                                hpt = Pt_Nch_hist->ProjectionY(Form("%.3f-%.3f%% [%s]", params.GetFineAmplitudeCutLow(0)*100,
                                                                    params.GetFineAmplitudeCutHigh(k)*100,
                                                                    funcStr.Data()),
                                                               params.GetFineAmplitudeBinCutLow(k),
                                                               params.GetFineAmplitudeBinCutHigh(0));
                                hpt->SetTitle(Form("%.0f-%.4f%%", params.GetFineAmplitudeCutLow(0)*100,
                                                   params.GetFineAmplitudeCutHigh(k)*100));
                                hMpt = Pt_Mult_hist->ProjectionY(Form("%.3f-%.3f%%", params.GetFineAmplitudeCutLow(0)*100,
                                                                      params.GetFineAmplitudeCutHigh(k)*100),
                                                                 params.GetFineAmplitudeBinCutLow(k),
                                                                 params.GetFineAmplitudeBinCutHigh(0));
                                hMpt->SetTitle(Form("%.0f-%.4f%%", params.GetFineAmplitudeCutLow(0)*100,
                                                   params.GetFineAmplitudeCutHigh(k)*100));
                                hNch = Nch_Mult_hist->ProjectionY(Form("%.3f-%.3f%%", params.GetFineAmplitudeCutLow(0)*100,
                                                                      params.GetFineAmplitudeCutHigh(k)*100),
                                                                 params.GetFineAmplitudeBinCutLow(k),
                                                                 params.GetFineAmplitudeBinCutHigh(0));
                                hNch->SetTitle(Form("%.0f-%.4f%%", params.GetFineAmplitudeCutLow(0)*100,
                                                   params.GetFineAmplitudeCutHigh(k)*100));
                            }
                            else{                                      //Exclusive Centrality Analysis
                                hpt = Pt_Nch_hist->ProjectionY(Form("%.3f-%.3f%% [%s]", params.GetFineAmplitudeCutLow(k)*100,
                                                                    params.GetFineAmplitudeCutHigh(k)*100,
                                                                    funcStr.Data()),
                                                               params.GetFineAmplitudeBinCutLow(k),
                                                               params.GetFineAmplitudeBinCutHigh(k));
                                hpt->SetTitle(Form("%.3f-%.3f%%",
                                                   params.GetFineAmplitudeCutLow(k)*100,
                                                   params.GetFineAmplitudeCutHigh(k)*100));
                                hMpt = Pt_Mult_hist->ProjectionY(Form("%.3f-%.3f%%", params.GetFineAmplitudeCutLow(k)*100,
                                                                      params.GetFineAmplitudeCutHigh(k)*100),
                                                                 params.GetFineAmplitudeBinCutLow(k),
                                                                 params.GetFineAmplitudeBinCutHigh(k));
                                hMpt->SetTitle(Form("%.0f-%.4f%%", params.GetFineAmplitudeCutLow(0)*100,
                                                   params.GetFineAmplitudeCutHigh(k)*100));
                                hNch = Nch_Mult_hist->ProjectionY(Form("%.3f-%.3f%%", params.GetFineAmplitudeCutLow(k)*100,
                                                                      params.GetFineAmplitudeCutHigh(k)*100),
                                                                 params.GetFineAmplitudeBinCutLow(k),
                                                                 params.GetFineAmplitudeBinCutHigh(k));
                                hNch->SetTitle(Form("%.0f-%.4f%%", params.GetFineAmplitudeCutLow(k)*100,
                                                   params.GetFineAmplitudeCutHigh(k)*100));
                            }
                            
                            hpt->Sumw2();
                            hpt->SetYTitle("1/N_{evt} 1/(2#pi p_{T}) d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-2}"); // Rutik/Rene's version (differential yield)
                            hMpt->Sumw2();
                            hMpt->SetYTitle("Normalized Number of Events");
                            hMpt->SetXTitle("[p_{T}]");
                            hNch->Sumw2();
                            hNch->SetYTitle("Normalized Number of Events");
                            hNch->SetXTitle("N_{ch}");
                            
                            params.NormalizeByBinWidth(hpt);
                            if (params.GetCentralitySetting() == "I"){ //Inclusive Centrality Analysis
                                hpt->Scale(1./Nch_dist->Integral(params.GetFineAmplitudeBinCutLow(k),
                                                                 params.GetFineAmplitudeBinCutHigh(0))); //Inclusive
                                hMpt->Scale(1./Nch_dist->Integral(params.GetFineAmplitudeBinCutLow(k),
                                                                 params.GetFineAmplitudeBinCutHigh(0)));
                                hNch->Scale(1./Nch_dist->Integral(params.GetFineAmplitudeBinCutLow(k),
                                                                 params.GetFineAmplitudeBinCutHigh(0)));
                            }
                            else{                                      //Exclusive Centrality Analysis
                                hpt->Scale(1./Nch_dist->Integral(params.GetFineAmplitudeBinCutLow(k),
                                                                 params.GetFineAmplitudeBinCutHigh(k))); //Exclusive
                                hMpt->Scale(1./Nch_dist->Integral(params.GetFineAmplitudeBinCutLow(k),
                                                                 params.GetFineAmplitudeBinCutHigh(k)));
                                hNch->Scale(1./Nch_dist->Integral(params.GetFineAmplitudeBinCutLow(k),
                                                                 params.GetFineAmplitudeBinCutHigh(k)));
                            }
                            params.SetErrorYield(hpt, k);
                            
                            bool ifpt_cut = true; //Enforcing 10.0 GeV and 200 MeV cutoff
                            tf1Name = Form("hfit_%s", Parameterization::FunctionToString(fitFuncType).c_str());
                            fit = params.GetFit(hpt, fit_range.first, fit_range.second, fitFuncType, tf1Name);
                            chi2 = fit->GetChisquare();
                            ndf     = fit->GetNDF();
                            chi2ndf = (ndf > 0) ? chi2 / ndf : 0.0;
                            if (k < 20 && (chi2ndf > 1000. || chi2ndf < 1e-9)) continue;
                            if (params.GetDebugSetting()) {
                                std::cout << "Writing hpt";
                                cout << "\n Index: " << k;
                                cout << "\n first: " << first << std::endl;
                            }
                            if (first){
                                spectrafitRangeSubdir->cd();
                                hpt->SetMarkerStyle(4); hpt->SetMarkerColor(kBlack); hpt->SetMarkerSize(1.2);
                                hpt->Write();
                                ptStudyDir_Dist->cd();
                                hMpt->Write();
                                NchStudyDir_Dist->cd();
                                hNch->Write();
                            }
                            
                            //-----Quantity calculation-----
                            double dndetaAna{params.GetdNdEta(hpt, fit, ifpt_cut, high_cut,low_cut)};
                            double error_dndetaAna{params.GetErrordNdEta(hpt, fit, true, 0.2)};
                            double error_dndetaAna_ref{params.GetErrordNdEta(hpt_ref, fit, true, 0.2)};
                            double error_dndeta_ratio{( dndetaAna / dndetaAna_ref) * sqrt(pow(( (error_dndetaAna / dndetaAna) - (error_dndetaAna_ref / dndetaAna_ref) ), 2) )};
                            
                            h->pTC1 = params.GetCumulant(hMpt, 1);
                            h->pTC2 = params.GetCumulant(hMpt, 2);
                            h->pTC3 = params.GetCumulant(hMpt, 3);
                            h->pTC4 = params.GetCumulant(hMpt, 4);
                            h->pTM = h->pTC1;
                            h->pTVar = h->pTC2;
                            h->pTS = h->pTC3 / std::pow(h->pTC2, 3. / 2.);
                            h->pTkap = h->pTC4 / std::pow(h->pTC2, 2.);

                            h->NchC1 = params.GetCumulant(hNch, 1);
                            h->NchC2 = params.GetCumulant(hNch, 2);
                            h->NchC3 = params.GetCumulant(hNch, 3);
                            h->NchC4 = params.GetCumulant(hNch, 4);
                            h->NchM = h->NchC1;
                            h->NchVar = h->NchC2;
                            h->NchS = h->NchC3 / std::pow(h->NchC2, 3. / 2.);
                            h->Nchkap = h->NchC4 / std::pow(h->NchC2, 2.);
                            
                            h->pTC1_error = params.GetCumulantError(hMpt, 1);
                            h->pTC2_error = params.GetCumulantError(hMpt, 2);
                            h->pTC3_error = params.GetCumulantError(hMpt, 3);
                            h->pTC4_error = params.GetCumulantError(hMpt, 4);
                            h->pTM_error   = h->pTC1_error;
                            h->pTVar_error = h->pTC2_error;
                            h->pTS_error   = std::sqrt( std::pow(h->pTC3_error / std::pow(h->pTC2, 1.5), 2) +
                                                        std::pow(1.5 * h->pTC3 * h->pTC2_error / std::pow(h->pTC2, 2.5), 2) );
                            h->pTkap_error = std::sqrt( std::pow(h->pTC4_error / std::pow(h->pTC2, 2), 2) +
                                                        std::pow(2 * h->pTC4 * h->pTC2_error / std::pow(h->pTC2, 3), 2) );

                            h->NchC1_error = params.GetCumulantError(hNch, 1);
                            h->NchC2_error = params.GetCumulantError(hNch, 2);
                            h->NchC3_error = params.GetCumulantError(hNch, 3);
                            h->NchC4_error = params.GetCumulantError(hNch, 4);
                            h->NchM_error   = h->NchC1_error;
                            h->NchVar_error = h->NchC2_error;
                            h->NchS_error   = std::sqrt( std::pow(h->NchC3_error / std::pow(h->NchC2, 1.5), 2) +
                                                         std::pow(1.5 * h->NchC3 * h->NchC2_error / std::pow(h->NchC2, 2.5), 2) );
                            h->Nchkap_error = std::sqrt( std::pow(h->NchC4_error / std::pow(h->NchC2, 2), 2) +
                                                         std::pow(2 * h->NchC4 * h->NchC2_error / std::pow(h->NchC2, 3), 2) );

                            
                            if (params.GetDebugSetting()) std::cout << "Setting points k = " << k << std::endl;
                            // X-axis
                            double xval = dndetaAna / dndetaAna_ref;
                            double xerr = error_dndeta_ratio;

                            // ----------------- Cumulants (Absolute) -----------------
                            h->gpTC1->SetPoint(k, xval, h->pTC1);
                            h->gpTC1->SetPointError(k, xerr, h->pTC1_error);

                            h->gpTC2->SetPoint(k, xval, h->pTC2);
                            h->gpTC2->SetPointError(k, xerr, h->pTC2_error);

                            h->gpTC3->SetPoint(k, xval, h->pTC3);
                            h->gpTC3->SetPointError(k, xerr, h->pTC3_error);

                            h->gpTC4->SetPoint(k, xval, h->pTC4);
                            h->gpTC4->SetPointError(k, xerr, h->pTC4_error);

                            // ----------------- Cumulants (Normalized to 0-5%) -----------------
                            h->gpTC1N->SetPoint(k, xval, h->pTC1 / h->pTC1_ref);
                            h->gpTC1N->SetPointError(k, xerr,
                                std::sqrt(
                                    std::pow(h->pTC1_error / h->pTC1_ref, 2) +
                                    std::pow(h->pTC1 * h->pTC1_ref_error / std::pow(h->pTC1_ref, 2), 2)
                                )
                            );

                            h->gpTC2N->SetPoint(k, xval, h->pTC2 / h->pTC2_ref);
                            h->gpTC2N->SetPointError(k, xerr,
                                std::sqrt(
                                    std::pow(h->pTC2_error / h->pTC2_ref, 2) +
                                    std::pow(h->pTC2 * h->pTC2_ref_error / std::pow(h->pTC2_ref, 2), 2)
                                )
                            );

                            h->gpTC3N->SetPoint(k, xval, h->pTC3 / h->pTC3_ref);
                            h->gpTC3N->SetPointError(k, xerr,
                                std::sqrt(
                                    std::pow(h->pTC3_error / h->pTC3_ref, 2) +
                                    std::pow(h->pTC3 * h->pTC3_ref_error / std::pow(h->pTC3_ref, 2), 2)
                                )
                            );

                            h->gpTC4N->SetPoint(k, xval, h->pTC4 / h->pTC4_ref);
                            h->gpTC4N->SetPointError(k, xerr,
                                std::sqrt(
                                    std::pow(h->pTC4_error / h->pTC4_ref, 2) +
                                    std::pow(h->pTC4 * h->pTC4_ref_error / std::pow(h->pTC4_ref, 2), 2)
                                )
                            );

                            // ----------------- Moments (Absolute) -----------------
                            h->gpTM->SetPoint(k, xval, h->pTM);
                            h->gpTM->SetPointError(k, xerr, h->pTM_error);

                            h->gpTVar->SetPoint(k, xval, h->pTVar);
                            h->gpTVar->SetPointError(k, xerr, h->pTVar_error);

                            h->gpTS->SetPoint(k, xval, h->pTS);
                            h->gpTS->SetPointError(k, xerr, h->pTS_error);

                            h->gpTkap->SetPoint(k, xval, h->pTkap);
                            h->gpTkap->SetPointError(k, xerr, h->pTkap_error);

                            // ----------------- Moments (Normalized to 0-5%) -----------------
                            h->gpTMN->SetPoint(k, xval, h->pTM / h->pTM_ref);
                            h->gpTMN->SetPointError(k, xerr,
                                std::sqrt(
                                    std::pow(h->pTM_error / h->pTM_ref, 2) +
                                    std::pow(h->pTM * h->pTM_ref_error / std::pow(h->pTM_ref, 2), 2)
                                )
                            );

                            h->gpTVarN->SetPoint(k, xval, h->pTVar / h->pTVar_ref);
                            h->gpTVarN->SetPointError(k, xerr,
                                std::sqrt(
                                    std::pow(h->pTVar_error / h->pTVar_ref, 2) +
                                    std::pow(h->pTVar * h->pTVar_ref_error / std::pow(h->pTVar_ref, 2), 2)
                                )
                            );

                            h->gpTSN->SetPoint(k, xval, h->pTS / h->pTS_ref);
                            h->gpTSN->SetPointError(k, xerr,
                                std::sqrt(
                                    std::pow(h->pTS_error / h->pTS_ref, 2) +
                                    std::pow(h->pTS * h->pTS_ref_error / std::pow(h->pTS_ref, 2), 2)
                                )
                            );

                            h->gpTkapN->SetPoint(k, xval, h->pTkap / h->pTkap_ref);
                            h->gpTkapN->SetPointError(k, xerr,
                                std::sqrt(
                                    std::pow(h->pTkap_error / h->pTkap_ref, 2) +
                                    std::pow(h->pTkap * h->pTkap_ref_error / std::pow(h->pTkap_ref, 2), 2)
                                )
                            );
                            
                            // ----------------- Cumulants (Absolute) -----------------
                            h->gNchC1->SetPoint(k, xval, h->NchC1);
                            h->gNchC1->SetPointError(k, xerr, h->NchC1_error);

                            h->gNchC2->SetPoint(k, xval, h->NchC2);
                            h->gNchC2->SetPointError(k, xerr, h->NchC2_error);

                            h->gNchC3->SetPoint(k, xval, h->NchC3);
                            h->gNchC3->SetPointError(k, xerr, h->NchC3_error);

                            h->gNchC4->SetPoint(k, xval, h->NchC4);
                            h->gNchC4->SetPointError(k, xerr, h->NchC4_error);

                            // ----------------- Cumulants (Normalized to 0-5%) -----------------
                            h->gNchC1N->SetPoint(k, xval, h->NchC1 / h->NchC1_ref);
                            h->gNchC1N->SetPointError(k, xerr,
                                std::sqrt(
                                    std::pow(h->NchC1_error / h->NchC1_ref, 2) +
                                    std::pow(h->NchC1 * h->NchC1_ref_error / std::pow(h->NchC1_ref, 2), 2)
                                )
                            );

                            h->gNchC2N->SetPoint(k, xval, h->NchC2 / h->NchC2_ref);
                            h->gNchC2N->SetPointError(k, xerr,
                                std::sqrt(
                                    std::pow(h->NchC2_error / h->NchC2_ref, 2) +
                                    std::pow(h->NchC2 * h->NchC2_ref_error / std::pow(h->NchC2_ref, 2), 2)
                                )
                            );

                            h->gNchC3N->SetPoint(k, xval, h->NchC3 / h->NchC3_ref);
                            h->gNchC3N->SetPointError(k, xerr,
                                std::sqrt(
                                    std::pow(h->NchC3_error / h->NchC3_ref, 2) +
                                    std::pow(h->NchC3 * h->NchC3_ref_error / std::pow(h->NchC3_ref, 2), 2)
                                )
                            );

                            h->gNchC4N->SetPoint(k, xval, h->NchC4 / h->NchC4_ref);
                            h->gNchC4N->SetPointError(k, xerr,
                                std::sqrt(
                                    std::pow(h->NchC4_error / h->NchC4_ref, 2) +
                                    std::pow(h->NchC4 * h->NchC4_ref_error / std::pow(h->NchC4_ref, 2), 2)
                                )
                            );

                            // ----------------- Moments (Absolute) -----------------
                            h->gNchM->SetPoint(k, xval, h->NchM);
                            h->gNchM->SetPointError(k, xerr, h->NchM_error);

                            h->gNchVar->SetPoint(k, xval, h->NchVar);
                            h->gNchVar->SetPointError(k, xerr, h->NchVar_error);

                            h->gNchS->SetPoint(k, xval, h->NchS);
                            h->gNchS->SetPointError(k, xerr, h->NchS_error);

                            h->gNchkap->SetPoint(k, xval, h->Nchkap);
                            h->gNchkap->SetPointError(k, xerr, h->Nchkap_error);

                            // ----------------- Moments (Normalized to 0-5%) -----------------
                            h->gNchMN->SetPoint(k, xval, h->NchM / h->NchM_ref);
                            h->gNchMN->SetPointError(k, xerr,
                                std::sqrt(
                                    std::pow(h->NchM_error / h->NchM_ref, 2) +
                                    std::pow(h->NchM * h->NchM_ref_error / std::pow(h->NchM_ref, 2), 2)
                                )
                            );

                            h->gNchVarN->SetPoint(k, xval, h->NchVar / h->NchVar_ref);
                            h->gNchVarN->SetPointError(k, xerr,
                                std::sqrt(
                                    std::pow(h->NchVar_error / h->NchVar_ref, 2) +
                                    std::pow(h->NchVar * h->NchVar_ref_error / std::pow(h->NchVar_ref, 2), 2)
                                )
                            );

                            h->gNchSN->SetPoint(k, xval, h->NchS / h->NchS_ref);
                            h->gNchSN->SetPointError(k, xerr,
                                std::sqrt(
                                    std::pow(h->NchS_error / h->NchS_ref, 2) +
                                    std::pow(h->NchS * h->NchS_ref_error / std::pow(h->NchS_ref, 2), 2)
                                )
                            );

                            h->gNchkapN->SetPoint(k, xval, h->Nchkap / h->Nchkap_ref);
                            h->gNchkapN->SetPointError(k, xerr,
                                std::sqrt(
                                    std::pow(h->Nchkap_error / h->Nchkap_ref, 2) +
                                    std::pow(h->Nchkap * h->Nchkap_ref_error / std::pow(h->Nchkap_ref, 2), 2)
                                )
                            );

                            if (params.GetDebugSetting()) std::cout << "Finished setting points" << std::endl;
                        } //End of index loop
        // -----**-----**-----**-----**-----**--(5)--**-----**-----**-----**-----**-----**-----**
        if (params.GetDebugSetting()) {cout << "\n#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#\n" << endl;}

        if (params.GetDebugSetting()) std::cout << "Writing absolute cumulants" << std::endl;
        ptStudyDir_Cumulants_Abs->cd();
        h->gpTC1->Write();
        h->gpTC2->Write();
        h->gpTC3->Write();
        h->gpTC4->Write();
    
        if (params.GetDebugSetting()) std::cout << "Writing normalized cumulants" << std::endl;
        ptStudyDir_Cumulants_Norm->cd();
        h->gpTC1N->Write();
        h->gpTC2N->Write();
        h->gpTC3N->Write();
        h->gpTC4N->Write();
    
        if (params.GetDebugSetting()) std::cout << "Writing absolute moments" << std::endl;
        ptStudyDir_Moments_Abs->cd();
        h->gpTM->Write();
        h->gpTVar->Write();
        h->gpTS->Write();
        h->gpTkap->Write();
    
        if (params.GetDebugSetting()) std::cout << "Writing normalized moments" << std::endl;
        ptStudyDir_Moments_Norm->cd();
        h->gpTMN->Write();
        h->gpTVarN->Write();
        h->gpTSN->Write();
        h->gpTkapN->Write();
                     
        if (params.GetDebugSetting()) std::cout << "Writing absolute cumulants" << std::endl;
        NchStudyDir_Cumulants_Abs->cd();
        h->gNchC1->Write();
        h->gNchC2->Write();
        h->gNchC3->Write();
        h->gNchC4->Write();

        if (params.GetDebugSetting()) std::cout << "Writing normalized cumulants" << std::endl;
        NchStudyDir_Cumulants_Norm->cd();
        h->gNchC1N->Write();
        h->gNchC2N->Write();
        h->gNchC3N->Write();
        h->gNchC4N->Write();

        if (params.GetDebugSetting()) std::cout << "Writing absolute moments" << std::endl;
        NchStudyDir_Moments_Abs->cd();
        h->gNchM->Write();
        h->gNchVar->Write();
        h->gNchS->Write();
        h->gNchkap->Write();

        if (params.GetDebugSetting()) std::cout << "Writing normalized moments" << std::endl;
        NchStudyDir_Moments_Norm->cd();
        h->gNchMN->Write();
        h->gNchVarN->Write();
        h->gNchSN->Write();
        h->gNchkapN->Write();

        gROOT->SetBatch(kFALSE);
        std::cout << "\nWriting final directories..." << std::endl;
        
        generalDir->cd();
        lOutGen->Write();
        fOut->cd();
        fOut->Write();
        
        std::cout << "deleting files..." << std::endl;
        delete h;
        delete fOut;
        gROOT->SetBatch(kFALSE);
    
}

double Parameterization::GetCumulant(TH1D* dist, int order){
    bool debug = false;
    if (debug) { std::cout << "Inside GetCumulant()" << std::endl; }
    //https://arxiv.org/pdf/2101.12413 page 9
    if (!dist || order < 1 || order > 4) return 0.0;
    if (debug) { std::cout << "Inside GetCumulant()" << std::endl; }
        double mean = dist->GetMean(); // ⟨N⟩
        double cumulant = 0.0;

    if (debug) { std::cout << "Inside GetCumulant()" << std::endl; }
        switch (order) {
            case 1: // C1 = ⟨N⟩
                if (debug) { std::cout << "Inside GetCumulant()" << std::endl; }
                cumulant = mean;
                break;

            case 2: // C2 = ⟨(δN)^2⟩ = μ2 = variance
                if (debug) { std::cout << "Inside GetCumulant()" << std::endl; }
                cumulant = dist->GetRMS() * dist->GetRMS();
                break;

            case 3: { // C3 = μ3 = ⟨(δN)^3⟩
                if (debug) { std::cout << "Inside GetCumulant()" << std::endl; }
                double mu3 = 0.0;
                int nBins = dist->GetNbinsX();
                for (int i = 1; i <= nBins; ++i) {
                    double x = dist->GetBinCenter(i);
                    double y = dist->GetBinContent(i);
                    mu3 += y * std::pow(x - mean, 3);
                }
                cumulant = mu3 / dist->GetEntries();
                break;
            }

            case 4: { // C4 = μ4 − 3μ2^2
                if (debug) { std::cout << "Inside GetCumulant()" << std::endl; }
                double mu2 = 0.0, mu4 = 0.0;
                int nBins = dist->GetNbinsX();
                for (int i = 1; i <= nBins; ++i) {
                    double x = dist->GetBinCenter(i);
                    double y = dist->GetBinContent(i);
                    mu2 += y * std::pow(x - mean, 2);
                    mu4 += y * std::pow(x - mean, 4);
                }
                mu2 /= dist->GetEntries();
                mu4 /= dist->GetEntries();
                cumulant = mu4 - 3.0 * std::pow(mu2, 2);
                break;
            }

            default:
                if (debug) { std::cout << "Inside GetCumulant()" << std::endl; }
                cumulant = 0.0;
                break;
        }
        if (debug) { std::cout << "Returning cumulant" << std::endl; }
        return cumulant;
}

double Parameterization::GetCumulantError(TH1D* dist, int order) {
    if (!dist || order < 1 || order > 4) return 0.0;

    double mean = dist->GetMean();
    double entries = dist->GetEntries();
    if (entries <= 1) return 0.0; // Avoid division by 0

    double mu2 = 0.0, mu3 = 0.0, mu4 = 0.0, mu6 = 0.0;
    int nBins = dist->GetNbinsX();

    for (int i = 1; i <= nBins; ++i) {
        double x = dist->GetBinCenter(i);
        double y = dist->GetBinContent(i);
        double dx = x - mean;
        mu2 += y * std::pow(dx, 2);
        mu3 += y * std::pow(dx, 3);
        mu4 += y * std::pow(dx, 4);
        mu6 += y * std::pow(dx, 6);  // Needed for C4 error
    }

    mu2 /= entries;
    mu3 /= entries;
    mu4 /= entries;
    mu6 /= entries;

    double error = 0.0;

    switch (order) {
        case 1: {
            double meanError = dist->GetMeanError();
            error = meanError;
            break;
        }
        case 2: {
            // Var(C2) ≈ (mu4 - mu2^2) / N
            error = std::sqrt((mu4 - mu2 * mu2) / entries);
            break;
        }
        case 3: {
            // Var(C3) ≈ (mu6 - mu3^2) / N
            error = std::sqrt((mu6 - mu3 * mu3) / entries);
            break;
        }
        case 4: {
            // Var(C4) = Var(mu4 - 3 * mu2^2)
            //         = Var(mu4) + 9 * Var(mu2^2) - 6 * Cov(mu4, mu2^2)
            // Approximate: Var(C4) ≈ (mu8 - (mu4 - 3*mu2^2)^2) / N
            // But that's complex; use:
            double var_mu4 = (mu6 - mu4 * mu4) / entries;
            double var_mu2 = (mu4 - mu2 * mu2) / entries;
            error = std::sqrt(var_mu4 + 36.0 * mu2 * mu2 * var_mu2);
            break;
        }
        default:
            error = 0.0;
    }

    return error;
}
