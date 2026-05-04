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

using namespace std;

enum KtBins {
  KTBIN1,
  KTBIN2,
  KTBIN3,
  NKTBINS
};

int Mbins = 4;
int EDbins = 1;
double NormLowEndPercent = 0.15; // percentage below upper limit
double NormHighPercent = 0; // percentage below upper limit
string TwoPartDir = "TwoParticle/";

//--------------------------------------------------------------------------------
std::pair<int, int> GetKtBin( TH2D *hist, int kTchoice ) {
  
  double kTlow=0, kThigh=0;
  
  if( kTchoice == KtBins::KTBIN1 ) {
    kTlow = 0.15001;
    kThigh = 0.24999;
  }
  if( kTchoice == KtBins::KTBIN2 ) {
    kTlow = 0.25001;
    kThigh = 0.34999;
  }
  if( kTchoice == KtBins::KTBIN3 ) {
    kTlow = 0.34001;
    kThigh = 0.59999;
  }

  int kTlowBin = hist->GetXaxis()->FindBin( kTlow );
  int kThighBin = hist->GetXaxis()->FindBin( kThigh );

  return std::make_pair( kTlowBin, kThighBin );
}

//--------------------------------------------------------------------------------
void assembleRawCorrelations( string fname = "../farmOutput/data/jobsOut_3p0_2018_v1.2.root" ) 
{

  TFile *f = new TFile( fname.data(), "READ" );

  vector<TH1D*> C2;

  for( int m = 0; m < Mbins; m++ ) {
    for( int ed = 0; ed < EDbins; ed++ ) {

      for( int c1 = 0; c1 < 2; c1++ ) {
        for( int c2 = c1; c2 < 2; c2++ ) {
          string chargeComb;
          if( c1 == c2 ) { 
            if( c1 == 0 ) { chargeComb = "--"; }
            else { chargeComb = "++"; }
          }
          else { chargeComb = "+-"; }

          for( int kt = 0 ; kt < KtBins::NKTBINS; kt++ ) {

            string nameNum = TwoPartDir + "TwoParticle_Charge1_" + to_string( c1 ) + "_Charge2_" + to_string( c2 );
            nameNum += "_M_" + to_string( m ) + "_ED_" + to_string( ed ) + "_Term_1";
            TString nameDen = ( nameNum.data() );
            nameDen.ReplaceAll("Term_1", "Term_2");
            TH2D *num2D = (TH2D*)f->Get( nameNum.data() );
            TH2D *den2D = (TH2D*)f->Get( nameDen );

            int kTlowBin = 0, kThighBin = 0;
            std::tie( kTlowBin, kThighBin ) = GetKtBin( num2D, kt );

            string name1D = "C2_" + to_string(c1) + to_string(c2) + "_ed" + to_string(ed) + "_kT" + to_string(kt) + "_M" + to_string(m);
            string title = chargeComb + " centrality:" + to_string(m) + " kt:" + to_string(kt);
            TH1D *num1D = (TH1D*)(num2D->ProjectionY( name1D.data(), kTlowBin, kThighBin ));
            TH1D *den1D = (TH1D*)den2D->ProjectionY( nameDen + "_pro", kTlowBin, kThighBin );
            num1D->SetTitle( title.data() );

            double NormLowBin = num1D->GetXaxis()->FindBin( (1 - NormLowEndPercent) * num1D->GetXaxis()->GetXmax() );
            double NormHighBin = num1D->GetXaxis()->FindBin( (1 - NormHighEndPercent) * num1D->GetXaxis()->GetXmax() );

            den1D->Scale( num1D->Integral(NormLowBin, NormHighBin) / den1D->Integral(NormLowBin, NormHighBin) );
            num1D->Divide( den1D );

            C2.push_back( num1D );

          }
        }
      }
    }
  }
 
  TFile *fout = new TFile( "rawCorrelations.root", "RECREATE" );


  for( uint i = 0; i < C2.size(); i++ ) {
    C2[i]->Write();
  }

}
