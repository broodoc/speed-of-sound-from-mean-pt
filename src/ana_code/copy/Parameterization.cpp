//NOTE: 9/15 - NormalizebyBinWidth function: I added pT normalization (by center)
//      9/15 - Changed hpT_ref from 0-5% to MB to compare to PID curves
#include "Parameterization.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

//-----------Analysis Headers--------------
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
//-----------Analysis Headers--------------

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
    /// NOTES:
    /// - EVENTUALLY NEED to use definition from event shuffling from Govert Nijs (maincollect.cpp in Trajectum?) to get exact centrality bin cuts (decimal, not estimated integer cuts). ALICE has more statistics so they can afford to get away with approximating centrality via integers.
    /// - NEED to make sure the number of bins (nCentralFinePercentileBinning) defined in Parameterization.h correctly matches the number of elements in the referenced array (CentralFinePercentileBinning)! Otherwise centrality definition will be off
    /// - NEED TO PERFORM GLAUBER FIT TO hin FOR EVENT SELECTION (THIS WILL AVOID PILEUP - WHEN THE TPC MISIDENTIFIES EVENTS THAT HAPPEN CLOSE TO EACH OTHER AS HIGHLY CENTRAL EVENT) AND THEN PERFORM THE REST OF FINDPERCENTILES WITH A GLAUBER TH1
    
    string name = hin->GetName();
    if (name.length() >= 4) {
        // Replace the last four characters with "bin_cuts.txt"
        name.replace(name.length() - 4, 4, "bin_cuts");
    }
    std::string output_label = name + ".txt";
    std::ofstream outFile(output_label);
    
  int nBins{10};
  if (isCoarseBinning) {
    nBins = this->GetnCoarsePercentileBinning();
      cout << "Getting coarse percentile binning for " << output_label << "... " << '\n';
  } else {
    nBins = this->GetnFinePercentileBinning();
      cout << "Getting fine percentile binning for " << output_label << "... " << '\n';
  }
  double newBin[nBins + 1];
  for (int bin_y_perc = 0; bin_y_perc < nBins + 1; ++bin_y_perc) {
    newBin[bin_y_perc] = 0;
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
  double counter{0.0};

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
      counter += 100 * abs(high_perc - low_perc);

        if (debug) {cout << " Real Percentage: " << 100 * abs(high_perc - low_perc) << "% |"
            << " Amplitude range: " << low_edge << " - " << up_edge << " |"
            << " Measured percentage: " << 100 * intergral_tmp << "% |"
            << " Number of events: " << events << '\n';
            
            cout << " Percentage window: " << 100 * low_perc << " - "
            << 100 * high_perc << " % |"
            << " Percentage covered: " << counter << " %" << '\n'
            << "========================================================" << '\n' ;
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
      for (int bin_y = firstbin; bin_y >= 0; bin_y--) {
        intergral_tmp += hin->GetBinContent(bin_y) / integral;
        events += hin->GetBinContent(bin_y);
        if (intergral_tmp > abs(high_perc - low_perc)) {
          foundbin = bin_y;
          break;
        }
      }

      if (bin_y_perc == (nBins - 1)) {
        foundbin = 1;
      }

      newBin[bin_y_perc] = firstbin;

      double up_edge{hin->GetXaxis()->GetBinUpEdge(firstbin)};
      double low_edge{hin->GetXaxis()->GetBinLowEdge(foundbin)};
      counter += 100 * abs(high_perc - low_perc);

        if (debug) {cout << " Real Percentage: " << 100 * abs(high_perc - low_perc) << "% |"
            << " N_{charged} range: " << low_edge << " - " << up_edge << " |"
            << " Measured percentage: " << 100 * intergral_tmp << "% |"
            << " Number of events: " << events << '\n';
            
            cout << " Percentage window: " << 100 * low_perc << " - "
            << 100 * high_perc << " % |"
            << " Percentage covered: " << counter << " %" << '\n'
            << "======================="
            << "======================="
            << "======================="
            << "======================="
            << "=======================" << '\n';
        }
    }
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
      return (xx / sqrt(par[0] * par[0] + xx * xx)) *
             pow((1 + xx / par[1]), -par[2]) * par[3];
      break;
    case Function::FokkerPlanck:
      // m - par[0]
      // A - par[1]
      // b - par[2]
      // c - par[3]
      // d - par[4]
      // T - par[5]
      return par[1] *
             exp(-1 * (par[2] / par[5]) *
                 atan((sqrt(xx * xx + par[0] * par[0]) - par[0]) / par[2])) /
             pow(1 + pow((sqrt(xx * xx + par[0] * par[0]) - par[0]) / par[2],
                         par[4]),
                 par[3]);
      break;
    case Function::PowerLaw:
      // n - par[0]
      return (par[1] / pow(xx, par[0]));
      break;
    case Function::ExpDecay:
      // n - par[0]
      return (par[0] + par[1] * xx * xx) / exp(par[2] * xx);
          break;
      case Function::RightKneeFunction:
          // p_0 = par[0]
          // p_1 (cs2) = par[1]
          // p_0^(p_1) = p_0^(cs2)
      return par[0] * TMath::Power(xx, par[1]);
          break;
      }
  }


//-----------------------------------------------------

TF1* Parameterization::GetFit(TH1 *h){
    //Fitting for <pT> spectra
    h->GetXaxis()->SetRangeUser(-0.1, 10);  // Set the fit region of h
//    Parameterization ExpDecay(Function::ExpDecay);
//      TF1 *fit = new TF1("fit", ExpDecay, 0.40, 10.0, 3);
//      fit->SetParameters(250, 20.0, 50.0);  //
//      fit->SetParLimits(0, 0.1, 300);  // parameter 0 // //
//      fit->SetParLimits(1, 0.001, 0.1);  // parameter 1 //   //
//      fit->SetParLimits(2, 0.05, 100.0);  // parameter 1 //   //
    //-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//        Parameterization PowerLaw(Function::PowerLaw);
//          TF1 *fit = new TF1("fit", PowerLaw, -0.1, 10.0, 2);
//          fit->SetParameters(1000, -2.0);  //
//          fit->SetParLimits(0, 0.1, 5000);  // parameter 0 // //
//          fit->SetParLimits(1, -100.0, 100.0);  // parameter 1 //   //
    //-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//    Parameterization Hagedorn(Function::Hagedorn); //Using Hagedorn for pT extrapolation
//      TF1 *fit = new TF1("fit", Hagedorn, -0.1, 2.5, 4);
//      fit->SetParameters(200, 1.0, 2.0,500);  //
//      fit->SetParLimits(0, 0.1, 5000);  // parameter 0 // // m - par[0]
//      fit->SetParLimits(1, 0.05, 50.0);  // parameter 1 //   // pT0 - par[1]
//      fit->SetParLimits(2, 0.1, 50.0);  // parameter 2 //  // b - par[2]
//      fit->SetParLimits(3, 250, 1000);  // parameter 3 // // constant - par[3]
//    fit->SetParameters(10, 0.5, 0.15,60);  //
//    fit->SetParLimits(0, 0.1, 500);  // parameter 0 // // m - par[0]
//    fit->SetParLimits(1, 0.05, 5.0);  // parameter 1 //   // pT0 - par[1]
//    fit->SetParLimits(2, 0.1, 50.0);  // parameter 2 //  // b - par[2]
//    fit->SetParLimits(3, 250, 1000);  // parameter 3 // // constant - par[3]
    //-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//    Parameterization Tsallis(Function::Tsallis); //Using Tsallis for pT extrapolation
//      TF1 *fit = new TF1("fit", Tsallis, -0.1, 3.0, 3);
//      fit->SetParameters(2e6, 3.50, 0.15);  //
//      fit->SetParLimits(0, 1e6, 5e6);  // parameter 0 //    V - par[0] (Volume - Tsallis Normalizing factor)
//      fit->SetParLimits(1, 1.0, 5.0);  // parameter 1 //   // q - par[1] (entropy index)
//      fit->SetParLimits(2, 0.08, 1.0);  // parameter 2 //     T - par[2] (Effective Temp)
    //-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    Parameterization levyTsallis(Function::LevyTsallis); //Using levyTsallis for pT extrapolation
      //Levy_Tsallis from 0 to 400 MeV
    TF1 *fit = new TF1("fit", levyTsallis, -0.1, 10, 4);
    fit->SetParameters(0.07, 0.107105, 55, 190);  //
    fit->SetParLimits(0, -0.1, 0.1);  // parameter 0 //    m - par[0]
    fit->SetParLimits(1, 0.1, 0.1);  // parameter 1 //   C - par[1]
    fit->SetParLimits(2, 50, 70);  // parameter 2 //     n - par[2]
    fit->SetParLimits(3, 150, 270);  // parameter 3 // dN/dy - par[3]
    //hpt_ref:
//    fit->SetParameters(0.0, 0.09, 5.125, 265);  //
//    fit->SetParLimits(0, -0.001, 0.001);  // parameter 0 //    m - par[0]
//    fit->SetParLimits(1, 0.1, 0.1);  // parameter 1 //   C - par[1]
//    fit->SetParLimits(2, 5.11, 5.14);  // parameter 2 //     n - par[2]
//    fit->SetParLimits(3, 260, 270);  // parameter 3 // dN/dy - par[3]
    //-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    // Performing and logging fit ---------------------------------------------------------
    fit->SetNpx(50000);
    TFitResultPtr fitResult = h->Fit(fit, "SQ");  // "S" ensures the fit result is returned
    
    //Fitting for <pT> Baryon Stopping (high pT peak) in dN/deta
    
    //Fitting for <pT> Particle Production (low pT peak) in dN/deta
    
    //Fitting for <pT> mixed (1/2(<pT>_B + <pT>_P)
    
    struct stat info;
    if (stat("Fitlogs", &info) != 0) {
        // Directory doesn't exist, create it
        mkdir("Fitlogs", 0777);  // 0777 gives full permissions to the folder
    }
    mkdir("Fitlogs", 0777);
    std::string histName = h->GetName();
    std::string logFileName = "Fitlogs/" + histName + "_fitlog.txt";
    std::ofstream fitLog(logFileName);
    if (!fitLog.is_open()) {
        std::cerr << "Error opening fitlog.txt for writing." << std::endl;
        return fit;
    }
    cout << "Performing pT extrapolation fit procedure for " << histName << ". Fit details are stored in " << logFileName << endl;
    fitLog << "****************************************" << std::endl;
    fitLog << "Minimizer is Minuit2 / Migrad" << std::endl;
//    fitLog << "Chi2                      = " << fitResult->Chi2() << std::endl;
//    fitLog << "NDf                       = " << fitResult->Ndf() << std::endl;
//    fitLog << "Edm                       = " << fitResult->Edm() << std::endl;
//    fitLog << "NCalls                    = " << fitResult->NCalls() << std::endl;

    // Log the parameters and their errors
    for (int i = 0; i < fit->GetNpar(); ++i) {
           double lowLimit, highLimit;
           fit->GetParLimits(i, lowLimit, highLimit);  // Get the limits of parameter i

           fitLog << "p" << i << "                        = "
                  << fit->GetParameter(i) << "   +/-   "
                  << fit->GetParError(i);

           // Check if the parameter has limits
           if (lowLimit != highLimit) {  // If limits are set, they won't be equal
               fitLog << "  \t (limited: [" << lowLimit << ", " << highLimit << "])";
           }

           fitLog << std::endl;
       }
    fitLog << "\n\n---------------FIT FUNCTION------------------------\n";
    for (double x = 0; x < 2.01; x += 0.01) {
        fitLog << "1/N_events 1/(2pi p_{T}) d^2N_{ch}/d#etadp_{T} f(" << x << ") = " << levyTsallis(&x, fit->GetParameters()) << endl;
//        fitLog << "1/N_events 1/(2pi p_{T}) d^2N_{ch}/d#etadp_{T} f(" << x << ") = " << ExpDecay(&x, fit->GetParameters()) << endl;
    }
    fitLog << "\n\n---------------DATASET------------------------\n";
    for (double x = 0; x < 2.01; x += 0.01) {
        int bin = h->FindBin(x);
        double y_value = h->GetBinContent(bin);
//        fitLog << "1/N_events 1/(2pi p_{T}) d^2N_{ch}/d#etadp_{T} f(" << x << ") = " << levyTsallis(&x, fit->GetParameters()) << endl;
        fitLog << "1/N_events 1/(2pi p_{T}) d^2N_{ch}/d#etadp_{T} f(" << x << ") = " << y_value << endl;
    }
    fitLog.close();
    // Performing and logging fit end ------------------------------------------------
    
    return fit;
}

//-----------------------------------------------------

void Parameterization::NormalizeByBinWidth(TH1 *h) {
  //! Normalize by bin width
    /// All other pT spectra normalizations folded in here (except 1/Nevts, this is from -> scale(1/Nevts) in Ratio plot func
    //Note: Also Normalizing by bin center (pT) here
  for (int bin = 1; bin <= h->GetNbinsX(); bin++) {
    double content{h->GetBinContent(bin)};
    double center{h->GetBinCenter(bin)};
    double error{h->GetBinError(bin)};
    double width{h->GetBinWidth(bin)};
    h->SetBinContent(bin, content / ( width * center * 2 * M_PI)); // 1/2pi (azimuthal normalization), 1/pT (pT = center, Invariance under Lorentz boost), 1/dpT (dpT = width)
//    h->SetBinContent(bin, content / width );
//    h->SetBinContent(bin, content / center); // Why 1/pT -> Ensures result is invariant under Lorentz trans.
//    h->SetBinContent(bin, content / ( 2 * 3.14159265)); //Why 1/2pi -> Azimuthal coverage from integration when translating invariant yeild expression to transverse quantities (pT and eta). We assume coverage over entire azimuthal plane
//    h->SetBinError( bin, error );
//    h->SetBinError(bin, error / ( width * center * 2 * M_PI)); //Here is where error needs to be reduced
//      h->SetBinError(bin, 0); //Here is where error needs to be reduced
//    h->SetBinError(bin, error / width);
//    h->SetBinError(bin, error / center);
//    h->SetBinError(bin, error / ( 2 * 3.14159265));
  }
}

//-----------------------------------------------------

double Parameterization::GetMean(TH1 *h, TF1 *fit, bool usecut,
                                 const double &high_cut, const double &low_cut) {
  double integral{0};
//  double mean{h->GetMean()};
//    cout << "mean: " << mean << endl;
    double mean{0};
//    return mean; //Remove later
//    double center{h->GetBinCenter(bin)};
//    double width{h->GetBinWidth(bin)};
//    double content{h->GetBinContent(bin)};r
  for (int bin = 0; bin <= h->GetNbinsX(); ++bin) {
    double center{h->GetBinCenter(bin)};
    double width{h->GetBinWidth(bin)};
//    double content{fit->Eval(center)};
    double content{h->GetBinContent(bin)};
    if (usecut && center > high_cut) {
      continue;
    }
    if (usecut && center < low_cut) {
        continue;
      }
    integral += width * center * content;
    mean += center * center * (width * content);
//      integral += width * content;
//    mean += center * (width * content);
  }
  mean /= integral;
    cout << "mean pT: " << mean << endl;
  return mean;
}

//-----------------------------------------------------

double Parameterization::GetdNdEta(TH1 *h, TF1 *fit, bool usecut,
                                   const double &high_cut, const double &low_cut) {
  double integral{0};
//    double center{h->GetBinCenter(bin)};
//    double width{h->GetBinWidth(bin)};
//    double content{h->GetBinContent(bin)};
    
  for (int bin = 0; bin <= h->GetNbinsX(); ++bin) {
    double center{h->GetBinCenter(bin)};
    double width{h->GetBinWidth(bin)};
    double content{h->GetBinContent(bin)};
//    double content{fit->Eval(center)};
    if (usecut && center > high_cut) {
      continue;
    }
      if (usecut && center < low_cut){
          continue;
      }
    integral += width * center * content;
//      integral += width * content;
  }
//    cout << h->GetName() << " N bins: " << h->GetNbinsX() << endl;
//    cout << "integral: " << integral << endl;
  return integral;
}

//-----------------------------------------------------

double Parameterization::GetNch(TH1 *h, const int low_ref, const int high_ref) {
//    Nch_dist,low_ref,high_ref
// This function is used for obtaining the mean Nch that corresponds to centrality class (as opposed to
// getdNdeta which produces the pT-integrated yield, in other words, the analysis corresponding
// Ncharged
//    h->Scale(1./h->Integral());s
    double integral{0};
  //  double mean{h->GetMean()};
  //    cout << "mean: " << mean << endl;
    double mean{0};
    int bin_low = h->GetXaxis()->FindBin(low_ref);
    int bin_high = h->GetXaxis()->FindBin(high_ref);
    if (bin_high < bin_low){return kFatal;}
    for (int bin = bin_low; bin <= bin_high; ++bin) {
      double center{h->GetBinCenter(bin)};
      double width{h->GetBinWidth(bin)};
  //    double content{fit->Eval(center)};
      double content{h->GetBinContent(bin)};
      integral += width * content;
      mean += center * (width * content);
  //      integral += width * content;
  //    mean += center * (width * content);
    }
    mean /= integral;
    cout << "Ncharged: " << mean << endl;
    return mean;
}

//-----------------------------------------------------

double Parameterization::GetErrordNdEta(const TH1 *h, bool usecut,
                                        const double &cut) {
  double uncertainty{0};
    // DO I NEED TO PERFORM FIT FOR ERROR? IF SO, HERE IS THE CODE:
//    Parameterization levyTsallis(Function::LevyTsallis); //Using levyTsallis for pT extrapolation
//      TF1 *fit = new TF1("fit", levyTsallis, 0, 10, 4);
//      fit->SetParameters(1.0, 0.1, 6.0, 10.0);  //
//      fit->SetParLimits(0, 0.01, 10);  // parameter 0
//      fit->SetParLimits(1, 0.001, 1);  // parameter 1
//      //save this into a GetMean_Spectra_Fit.txt:
//      for (double x = 0; x < 10; x += 0.1) {
//        std::cout << "GetMean: f(" << x << ") = " << levyTsallis(&x, fit->GetParameters()) << std::endl;
//      }
//    h->Fit(fit);
//    double center{h->GetBinCenter(bin)};
//    double width{h->GetBinWidth(bin)};
//    double error{fit->GetBinError(bin)};
    
  for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
    double center{h->GetBinCenter(bin)};
    double width{h->GetBinWidth(bin)};
    double error{h->GetBinError(bin)};
    if (usecut && center < cut) {
      continue;
    }
    uncertainty += pow(width, 2.0) * pow(error, 2.0);
  }
//    Remove later:
//    return 0;
    cout << "Error: " << uncertainty << endl;
  return sqrt(uncertainty);
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
    /*
  if (this->GetSystem() == "AuAu_200") {
    if (this->GetMultEstimator() == "V0") {
      perc = FinePercentileBinning[bin];
    } else {
      perc = CentralFinePercentileBinning[bin];
    }
  } else {
    perc = FinePercentileBinning[bin];
  }
     */
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
    /*
  if (this->GetSystem() == "AuAu_200") {
    if (this->GetMultEstimator() == "TPC") {
      perc = FinePercentileBinning[bin + 1];
    } else {
      perc = CentralFinePercentileBinning[bin + 1];
    }
  } else {
    perc = FinePercentileBinning[bin + 1];
  }
*/
  return perc;
}

//-----------------------------------------------------

int Parameterization::GetCoarseAmplitudeBinCutLowMC(int bin) const {
  // Implements the bin cuts for Ntracklets/Nch to calculate efficiencies. This
  // is included for the SPD Tracklets and TPC Tracks with eta gap, and the Half
  // TPC estimator.

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

void Parameterization::Ratio_compute_pt_Nch(TH2D* Pt_Nch_hist, TH1D* Nch_dist, bool isCoarseBinning, Parameterization params) {
    gROOT->SetBatch(kTRUE);
    string name = Pt_Nch_hist->GetName();
    if (name.length() >= 4) {
        // Replace the last four characters with "bin_cuts.txt"
        name.replace(name.length() - 4, 4, "Dec_Run21_from_raw_spectrum");
//        "Sep_25_from_spectrum_fit_LevyTsallis");
    }
    std::string output_label = "./graphs_" + name + ".root";
    gErrorIgnoreLevel = kError; /// Error handling
//    Parameterization params; //params
    params.SetMultEstimator("TPC"); /// params:Using the TPC as my mult estimator
    params.SetSystem("AuAu_200"); /// params:Setting system for Au+Au 200 GeV
    std::vector<Double_t> low_perc_bins; std::vector<Double_t> high_perc_bins;
    std::vector<int> low_bin_cuts; std::vector<int> high_bin_cuts;
    int nperc_bins = 0; bool debug = true;
    // -----------------------------
    
    TList* lOut = new TList();
    lOut->SetOwner(true);
    TGraphErrors* g = new TGraphErrors();
    g->SetName("norm_meanpt_vs_nch");
    g->SetMarkerStyle(20);
    g->GetXaxis()->SetTitle("<dN_{ch}/deta>/<dN_{ch}/deta>_{0-5%}");
    g->GetYaxis()->SetTitle("<p_{T}>/<p_{T}>_{0-5%}");
    TGraphErrors* gg = new TGraphErrors();
    gg->SetName("meanpt_vs_nch");
    gg->SetMarkerStyle(24);
    gg->GetXaxis()->SetTitle("N_{ch}");
    gg->GetYaxis()->SetTitle("<p_{T}>");
    TGraphErrors* glog = new TGraphErrors();
    glog->SetName("log_meanpt_vs_nch");
    glog->SetMarkerStyle(21);
    glog->GetXaxis()->SetTitle("log[dN_{ch}/d#eta]");
    glog->GetYaxis()->SetTitle("log[<p_{T}>]");
    lOut->Add(g);
    lOut->Add(gg);
    lOut->Add(glog);

    // -----------------------------
    
//    const int low_ref{params.GetFineAmplitudeBinCutLow(46)};
    ///Note: For low_ref 5% it's: elemental_index - 1 (ALWAYS)
    const int low_ref{params.GetFineAmplitudeBinCutLow(12)}; /// Obtaining reference bins ---> There are 49 elements (30 corresponds to 5% bin):: Here is where the data needs to be pulled from .txt
    const int high_ref{params.GetFineAmplitudeBinCutHigh(0)}; /// Obtaining reference bins:: Here is where the data needs to be pulled from .txt
    TH1D* hpt_ref = (TH1D*)Pt_Nch_hist->ProjectionY("hpt_ref",low_ref,high_ref);
    hpt_ref->Scale(1./Nch_dist->Integral(low_ref,high_ref)); //<-- Scales spectra by 1/events
//    hpt_ref->SetYTitle("1/N_{evt} d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-1}"); //
    hpt_ref->SetYTitle("1/N_{evt} 1/(2#pi p_{T}) d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-2}");
    //scale by pT binning and scale by 1/(2pi) here
    params.NormalizeByBinWidth(hpt_ref); //1/(2pi*pT) and 1/(dpT)
    lOut->Add(hpt_ref);
    
    cout << "Perc 0-5: " << Nch_dist->Integral(low_ref,high_ref) / Nch_dist->Integral(1,Nch_dist->GetNbinsX());
    if (debug){
        printf("\n\t low_ref = %d | high_ref = %d\n",low_ref,high_ref);
        printf("\n\t Number of FINE bins: %d\n",params.GetnFinePercentileBinning());
    }
    if (debug) {cout << "Fitting for 0-5%...\n";}
    TF1* fit{params.GetFit(hpt_ref)};
    const double meanpt_ref{params.GetMean(hpt_ref,fit,true,10.0,0.20)};
    const double dndeta_ref{params.GetNch(Nch_dist,low_ref,high_ref)};
//    const double dndeta_ref{params.GetdNdEta(hpt_ref,fit,true,10.0,0.20)};
    const double error_ref{params.GetErrordNdEta(hpt_ref,true,0.20)};
    
//    const double meanpt_ref{params.GetMean(hpt_ref,true,10.0)};
//    const double dndeta_ref{params.GetdNdEta(hpt_ref,true,10.0)};
//    const double error_ref{params.GetErrordNdEta(hpt_ref,true,10.0)};
    
    cout << '\n' << "<pT>_{0-5%}: " << meanpt_ref;
    cout << '\n' << "dn/deta_{0-5%}: " << dndeta_ref;
    
    /// Vector data acquisition for FINE and COARSE binning:
        //MAKE SURE TPC BIN CUTS ARE CORRECT (wrong cuts as of now)
        cout << '\n' << "===================================";
        cout << "===========================================" << '\n';
        cout << "\t RATIO PLOT WITH FINE BINNING PROCESSING... " << '\n';
        nperc_bins = params.GetnFinePercentileBinning();
        cout << "\t Number of bins: " << nperc_bins << '\n';
        cout << '\n' << "\t N events (MB): " << Nch_dist->Integral(0,high_ref) << '\n';
        
        if (debug) {cout << "\t Low Bin percentiles: ";}
        for (int i = 0; i < nperc_bins+1; i++){
            low_perc_bins.push_back(params.GetFineAmplitudeCutLow(i));
            low_bin_cuts.push_back(params.GetFineAmplitudeBinCutLow(i));
            if (debug) {cout << low_perc_bins[i]*100 << "%(" << low_bin_cuts[i] << "), ";}
        }
        if (debug) {cout << '\n' << "\t High Bin percentiles: ";}
        for (int j = 0; j < nperc_bins; j++){
            high_perc_bins.push_back(params.GetFineAmplitudeCutHigh(j));
            high_bin_cuts.push_back(params.GetFineAmplitudeBinCutHigh(j));
            if (debug) {cout << high_perc_bins[j]*100 << "%(" << high_bin_cuts[j] << "), "; }
        }
    cout << '\n';
        for (int k = 0; k < nperc_bins; k++) {
//            if (debug) {cout << "Fitting for " << Form("hpt_binindex_%d", k) << "...\n";}
            TH1D* hpt = Pt_Nch_hist->ProjectionY(Form("hpt_binindex_%d", k),params.GetFineAmplitudeBinCutLow(k),params.GetFineAmplitudeBinCutHigh(k));
//            hpt->SetYTitle("1/N_{evt} d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-1}"); //
            hpt->SetYTitle("1/N_{evt} 1/(2#pi p_{T}) d^{2}N_{ch}/d#eta dp_{T} (GeV/c)^{-2}");
            hpt->Scale(1./Nch_dist->Integral(params.GetFineAmplitudeBinCutLow(k),params.GetFineAmplitudeBinCutHigh(k)));
            params.NormalizeByBinWidth(hpt);
            lOut->Add(hpt);
            bool ifpt_cut = true; //Enforcing 10.0 GeV and 200 MeV cutoff
            fit = params.GetFit(hpt);
            
//            double mean_nch{params.GetdNdEta(hpt,fit,ifpt_cut,10.0,0.20)};
            double mean_nch{params.GetNch(Nch_dist,params.GetFineAmplitudeBinCutLow(k),params.GetFineAmplitudeBinCutHigh(k))};
            double mean_pt{params.GetMean(hpt,fit,ifpt_cut,10.0,0.20)};
            double error_pt{params.GetErrordNdEta(hpt,true,0.20)};
            double error_ratio{(mean_pt / meanpt_ref) * sqrt(pow(( (error_pt / mean_pt) - (error_ref / meanpt_ref) ), 2) )};
            // Check this: ^Error progation-> https://en.wikipedia.org/wiki/Propagation_of_uncertainty
//            double error_ratio{error_pt};
                
            g->SetPoint(k, mean_nch/dndeta_ref, mean_pt/meanpt_ref);
            g->SetPointError(k, 0.0, error_ratio);
            gg->SetPoint(k, mean_nch, mean_pt);
            gg->SetPointError(k, 0.0, error_pt);
            glog->SetPoint(k, log(mean_nch), log(mean_pt));
            glog->SetPointError(k, 0.0, error_pt/mean_pt);
        }
    params.saveToCSV(nperc_bins, high_bin_cuts, low_bin_cuts, low_perc_bins, high_perc_bins, Nch_dist, params);
//    Double_t cs2{params.GetFitcs2(g)};
//    TParameter<Double_t>* cs2Double = new TParameter<Double_t>("cs2", cs2);
//    TF1* gfit{params.GetFitcs2(g)};
    
//    if (gfit != nullptr) {lOut->Add(gfit); }
//    lOut->Add(gfit);
    // -----------------------------
        
    gErrorIgnoreLevel = kFatal; //Ignoring useless and nonFatal writing-based errors
    
    TFile* fOut = new TFile(output_label.c_str(),"recreate");
    
    fOut->Write();
    fOut->cd();
    lOut->Write();
    delete fOut;
    delete lOut;
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
    csvFile << "Bin Index,Bin Cuts,Centrality,N Events\n";

    // Write the data in columns
    for (int k = 1; k < nperc_bins; k++) {
        // bin index
        std::string bin_index = std::to_string(k-1);
        
        // bin cuts
        std::string bin_cuts = std::to_string(high_bin_cuts[k-1]) + " - " + std::to_string(low_bin_cuts[k-1]);

        // centrality
        std::stringstream centralityStream;
                centralityStream << std::fixed << std::setprecision(3)
                                 << low_perc_bins[k-1] * 100 << "% - "
                                 << high_perc_bins[k-1] * 100 << "%";
                std::string centrality = centralityStream.str();

        // N Events using the integral
        double n_events = Nch_dist->Integral(params.GetFineAmplitudeBinCutLow(k-1), params.GetFineAmplitudeBinCutHigh(k-1));

        // Write the row to the CSV file
        csvFile << bin_index << "," << bin_cuts << "," << centrality << "," << n_events << "\n";
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
TF1* Parameterization::GetFitcs2(TGraphErrors* g){
    double rkfx1{1.12}; double rkfx2{1.30};
    TString fitName = TString::Format("Right_knee_fit");
//    Parameterization RightKneeFunction(Function::RightKneeFunction); //Using levyTsallis for pT extrapolation
//      //Levy_Tsallis from 0 to 400 MeV
//    TF1 *fit = new TF1("fit", RightKneeFunction, rkfx1, rkfx2, 2);
//    fit->SetParameters(0.07, 0.107105, 55, 190);  //
//    fit->SetParLimits(0, -0.1, 0.1);  // parameter 0 //    m - par[0]
//    fit->SetParLimits(1, 0.1, 0.1);  // parameter 1 //   C - par[1]
//    
//    TF1 *Right_knee_fit = new TF1(fitName,Parameterization::RightKneeFunction,rkfx1, rkfx2,2); // 2 parameters (p_0, p_1) 0.90, 1.3,| 0.90, 1.4 |
//        Right_knee_fit->SetParLimits(0, 0.822, 1.5);
//        Right_knee_fit->SetParameter(1, 0.025);
//
//    double linewidth = 3.0;
//    Right_knee_fit->SetLineWidth(linewidth);
//    
//    TFitResultPtr fitResult = g->Fit(Right_knee_fit, "R0S", "");
//    g->Fit(Right_knee_fit, "R0", "");
//    g->SetLineColor(kRed);
//    Right_knee_fit->SetLineColor(kBlack);
//    Double_t cs2_1{0.0};
//    cs2_1 = Right_knee_fit->GetParameter(1);
//    Double_t cs2_error_1{0.0};
//    cs2_error_1 = Right_knee_fit->GetParError(1);
//    cout << "\t cs2: " << cs2_1 << " +\- " << cs2_error_1 <<endl;
//    
//    
//    std::string logFileName = "cs2log.txt";
//    std::ofstream fitLog(logFileName);
//    if (!fitLog.is_open()) {
//        std::cerr << "Error opening cs2log.txt for writing." << std::endl;
//    }
//    cout << "\tPerforming cs2 fit procedure. Fit details are stored in " << logFileName << endl;
//    fitLog << "****************************************" << std::endl;
//    fitLog << "Minimizer is Minuit2 / Migrad" << std::endl;
//    fitLog << "Chi2                      = " << fitResult->Chi2() << std::endl;
//    fitLog << "NDf                       = " << fitResult->Ndf() << std::endl;
//    fitLog << "Edm                       = " << fitResult->Edm() << std::endl;
//    fitLog << "NCalls                    = " << fitResult->NCalls() << std::endl;
//
//    // Log the parameters and their errors
//    for (int i = 0; i < Right_knee_fit->GetNpar(); ++i) {
//           double lowLimit, highLimit;
//        Right_knee_fit->GetParLimits(i, lowLimit, highLimit);  // Get the limits of parameter i
//
//           fitLog << "p" << i << "                        = "
//                  << Right_knee_fit->GetParameter(i) << "   +/-   "
//                  << Right_knee_fit->GetParError(i);
//
//           // Check if the parameter has limits
//           if (lowLimit != highLimit) {  // If limits are set, they won't be equal
//               fitLog << "  \t (limited: [" << lowLimit << ", " << highLimit << "])";
//           }
//
//           fitLog << std::endl;
//       }
//    fitLog.close();
//    
//    return Right_knee_fit;
    return nullptr;
}
// -----------------------------------------------------

// -----------------------------------------------------
