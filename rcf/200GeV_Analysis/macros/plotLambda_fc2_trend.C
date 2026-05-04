#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <Riostream.h>
#include <assert.h>

#include "TVector2.h"
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include "TF3.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TText.h"
#include "TRandom3.h"
#include "TArray.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TPaveStats.h"
#include "TASImage.h"

#define BohrR 1963.6885 // Bohr Radius for pions
#define FmToGeV 0.19733 // conversion to fm
#define mass_pion 0.13957 

using namespace std;



string fitName = "fit_01_ed0_kT0_M0";
string relPath = "../Results/";


//--------------------------------------------------------------------------------
void plotLambda_fc2_trend() 
{

  vector<pair<string,double>> list = {{"results_3p85_fxt_2018_v1.2.root",3.0}, {"results_4p59_fxt_2019_v1.2.root",3.2}, {"results_5p75_fxt_2020_v1.2.root",3.5}, {"results_7p3_fxt_2020_v1.2.root",3.9}, {"results_31p2_fxt_2020_v1.2.root",7.7}, {"results_200p0_col_2011_v1.3.root",200}};

  vector<double> ss_params;
  vector<double> ss_params_e;
  vector<double> ss_energies;
  
  vector<double> os_params;
  vector<double> os_params_e;
  vector<double> os_energies;

  // lambda
  // visually scanned from Grigory et al STAR paper proposal.
  vector<pair<double,double>> lambdas = {{3.0, 0.71},{3.2, 0.68},{3.5, 0.636},{3.9, 0.59},{7.7, 0.532},{11 ,0.508},{19,0.5},{27, 0.49},{39,0.49},{62, 0.47},{200, 0.45}};
  for( auto el : lambdas ) {
    ss_energies.push_back( el.first );
    ss_params.push_back( el.second );
    ss_params_e.push_back( 0 );
  }

  // fc2
  for( auto el : list ) {

    TFile *f = new TFile( (relPath + el.first).data(), "READ");
    TF1 *fit = (TF1*)f->Get( fitName.data() );
    os_params.push_back( fit->GetParameter(1) );
    os_params_e.push_back( fit->GetParError(1) );
    os_energies.push_back( el.second );
  }

  TGraphErrors *gr_lambda = new TGraphErrors( ss_energies.size(), ss_energies.data(), ss_params.data(), {}, ss_params_e.data() );
  gr_lambda->SetMarkerStyle(21);
  gr_lambda->SetMarkerColor(2);
  gr_lambda->SetMarkerSize(1.5);
  gr_lambda->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
  gr_lambda->GetXaxis()->SetTitleSize(0.06);
  gr_lambda->GetYaxis()->SetTitle("Dilution Parameter");
  gr_lambda->GetYaxis()->SetTitleSize(0.06);
  gr_lambda->SetTitle("");
  gr_lambda->SetMaximum(1.0);
  gr_lambda->SetMinimum(0.3);
  gr_lambda->Draw("ap");
 
  TGraphErrors *gr_fc2 = new TGraphErrors( os_energies.size(), os_energies.data(), os_params.data(), {}, os_params_e.data() );
  gr_fc2->SetMarkerStyle(20); // 29 (stars)
  gr_fc2->SetMarkerColor(4);
  gr_fc2->SetMarkerSize(1.5);
  gr_fc2->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
  gr_fc2->GetYaxis()->SetTitle("f_{c}^{2}");
  gr_fc2->SetMaximum(1.0);
  gr_fc2->SetMinimum(0.3);
  gr_fc2->Draw("p");

  gPad->SetLogx();
  gPad->SetGrid();

  TLegend *leg = new TLegend();
  leg->AddEntry(gr_lambda, "#lambda (#pi^{+} #pi^{+})", "p");
  leg->AddEntry(gr_fc2, "f_{c}^{2} (#pi^{-} #pi^{+})", "p");
  leg->Draw("same");

}
