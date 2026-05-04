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


int Npairs = 1e9;

double QinvMax = 0.1;

double Rmin = 3;
double Rmax = 10;
double Rstep = 1;
double fc2min = 0.3;
double fc2max = 0.8;
double fc2step = 0.1;


TH2D *mFSI_ss;
TH2D *mFSI_os;

void LoadFSIfile( TString filename );
double GetFemtoWeight( bool sameCharge, double fc2, double Rinv, double Qinv );
double GetFSI( bool sameCharge, double Rinv, double Qinv );

//--------------------------------------------------------------------------------
void makeMRC( string fname = "../farmOutput/embedding/embedding_3p85_fxt_2018_v1.0.root" ) 
{

  TRandom3 *rand = new TRandom3();
  rand->SetSeed(0);

  TF1 *LowQ_enhancer = new TF1("LowQ_enhancer","1 - [0]*x*x",0,0.1);
  LowQ_enhancer->SetParameter(0, 80);

  LoadFSIfile( "../corrections/KFile.root" );

  TFile *f = new TFile( fname.data(), "READ" );
  if( ! f ) {
    cout<<"source file not found"<<endl;
    return;
  }

  TH1D *Pt_plus = (TH1D*)f->Get("hMC_Pt_Plus");
  TH1D *Pt_minus = (TH1D*)f->Get("hMC_Pt_Minus");
  TH1D *Phi_plus = (TH1D*)f->Get("hMC_Phi_Plus");
  TH1D *Phi_minus = (TH1D*)f->Get("hMC_Phi_Minus");
  TH1D *Theta_plus = (TH1D*)f->Get("hMC_Theta_Plus");
  TH1D *Theta_minus = (TH1D*)f->Get("hMC_Theta_Minus");

  TH2D *Pt_Theta_plus = (TH2D*)f->Get("hMC_Pt_Theta_Plus");
  TH2D *Pt_Theta_minus = (TH2D*)f->Get("hMC_Pt_Theta_Minus");

  TH2D *dPt_plus2D = (TH2D*)f->Get("hMC_deltaPt_Plus");
  TH2D *dPt_minus2D = (TH2D*)f->Get("hMC_deltaPt_Minus");
  TH2D *dPhi_plus2D = (TH2D*)f->Get("hMC_deltaPhi_Plus");
  TH2D *dPhi_minus2D = (TH2D*)f->Get("hMC_deltaPhi_Minus");
  TH2D *dTheta_plus2D = (TH2D*)f->Get("hMC_deltaTheta_Plus");
  TH2D *dTheta_minus2D = (TH2D*)f->Get("hMC_deltaTheta_Minus");

  const int Nptbins = dPt_plus2D->GetNbinsX();
  TH1D *dPt[2][Nptbins];
  TH1D *dPhi[2][Nptbins];
  TH1D *dTheta[2][Nptbins];

  for( int i = 0; i < Nptbins; i++ ) {
    dPt[0][i] = (TH1D*)dPt_plus2D->ProjectionY(Form("dPt_plus%i",i), i+1, i+1);
    dPt[1][i] = (TH1D*)dPt_minus2D->ProjectionY(Form("dPt_minus%i",i), i+1, i+1);
    dPhi[0][i] = (TH1D*)dPhi_plus2D->ProjectionY(Form("dPhi_plus%i",i), i+1, i+1);
    dPhi[1][i] = (TH1D*)dPhi_minus2D->ProjectionY(Form("dPhi_minus%i",i), i+1, i+1);
    dTheta[0][i] = (TH1D*)dTheta_plus2D->ProjectionY(Form("dTheta_plus%i",i), i+1, i+1);
    dTheta[1][i] = (TH1D*)dTheta_minus2D->ProjectionY(Form("dTheta_minus%i",i), i+1, i+1);
  }

  double fc2HistMax = fc2max + fc2step/2.0;
  double fc2HistMin = fc2min - fc2step/2.0;
  int fc2Bins = round((fc2HistMax - fc2HistMin) / fc2step);
  double RHistMax = Rmax + Rstep/2.0;
  double RHistMin = Rmin - Rstep/2.0;
  int RBins = round((RHistMax - RHistMin) / Rstep);

  TH3D *ideal_ss = new TH3D("ideal_ss",";fc2;R_{inv} (fm);q_{inv} (GeV/c)", fc2Bins,fc2HistMin,fc2HistMax, RBins,RHistMin,RHistMax, 20,0,0.1);
  TH3D *ideal_os = new TH3D("ideal_os",";fc2;R_{inv} (fm);q_{inv} (GeV/c)", fc2Bins,fc2HistMin,fc2HistMax, RBins,RHistMin,RHistMax, 20,0,0.1);
  TH3D *ideal_noW = new TH3D("ideal_noW",";fc2;R_{inv} (fm);q_{inv} (GeV/c)", fc2Bins,fc2HistMin,fc2HistMax, RBins,RHistMin,RHistMax, 20,0,0.1);
  TH3D *smeared_ss = new TH3D("smeared_ss",";fc2;R_{inv} (fm);q_{inv} (GeV/c)", fc2Bins,fc2HistMin,fc2HistMax, RBins,RHistMin,RHistMax, 20,0,0.1);
  TH3D *smeared_os = new TH3D("smeared_os",";fc2;R_{inv} (fm);q_{inv} (GeV/c)", fc2Bins,fc2HistMin,fc2HistMax, RBins,RHistMin,RHistMax, 20,0,0.1);
  ideal_ss->Sumw2();
  smeared_ss->Sumw2();
  ideal_os->Sumw2();
  smeared_os->Sumw2();

  TH2D *Qcorrelation = new TH2D("Qcorrelation",";qTruth;qSmear", 20,0,0.1, 20,0,0.1);

  int seconds = time(0); 

  int pairN = 0;

  double pt_1 = 0;
  double theta_1 = 0;
  double phi_1 = 0;
  double pt_2 = 0;
  double theta_2 = 0;
  double phi_2 = 0;

  ////////////////////////////////////////////
  // pair loop
  while( pairN < Npairs ) {
  
    Pt_Theta_plus->GetRandom2( pt_1, theta_1 );
    phi_1 = Phi_plus->GetRandom();
    
    Pt_Theta_minus->GetRandom2( pt_2, theta_2 );
    phi_2 = Phi_minus->GetRandom();

    if( pt_1 < 0.15 || pt_2 < 0.15 ) { continue; }
    if( pt_1 > 0.55 || pt_2 > 0.55 ) { continue; }

    double px_1 = pt_1 * cos( phi_1 ); 
    double py_1 = pt_1 * sin( phi_1 );
    double pz_1 = pt_1 / tan( theta_1 );
    double px_2 = pt_2 * cos( phi_2 ); 
    double py_2 = pt_2 * sin( phi_2 );
    double pz_2 = pt_2 / tan( theta_2 );
    double E_1 = sqrt(pt_1*pt_1 + pz_1*pz_1 + mass_pion*mass_pion);
    double E_2 = sqrt(pt_2*pt_2 + pz_2*pz_2 + mass_pion*mass_pion);
    
    double qinvIdeal = sqrt( pow(px_1 - px_2, 2) + pow(py_1 - py_2, 2) + pow(pz_1 - pz_2, 2) - pow(E_1 - E_2, 2) );
    
    if( qinvIdeal > QinvMax ) { continue; }
    //if( LowQ_enhancer->Eval( qinvIdeal ) < rand->Rndm() ) { continue; }

    int ptIdx1 = dPt_plus2D->GetXaxis()->FindBin( pt_1 ) - 1;

    double dpx_1 = px_1 * dPt[0][ptIdx1]->GetRandom() - py_1 * dPhi[0][ptIdx1]->GetRandom();
    double dpy_1 = py_1 * dPt[0][ptIdx1]->GetRandom() + px_1 * dPhi[0][ptIdx1]->GetRandom();
    double dpz_1 = pz_1 * dPt[0][ptIdx1]->GetRandom() - pt_1/pow(sin(theta_1),2) * dTheta[0][ptIdx1]->GetRandom();
    
    int ptIdx2 = dPt_plus2D->GetXaxis()->FindBin( pt_2 ) - 1;
    
    int chIdx2 = 1; // use minus charge for second particle

    double dpx_2 = px_2 * dPt[chIdx2][ptIdx2]->GetRandom() - py_2 * dPhi[chIdx2][ptIdx2]->GetRandom();
    double dpy_2 = py_2 * dPt[chIdx2][ptIdx2]->GetRandom() + px_2 * dPhi[chIdx2][ptIdx2]->GetRandom();
    double dpz_2 = pz_2 * dPt[chIdx2][ptIdx2]->GetRandom() - pt_2/pow(sin(theta_2),2) * dTheta[chIdx2][ptIdx2]->GetRandom();


    double px1_smear = px_1 + dpx_1;
    double py1_smear = py_1 + dpy_1;
    double pz1_smear = pz_1 + dpz_1;
    double px2_smear = px_2 + dpx_2;
    double py2_smear = py_2 + dpy_2;
    double pz2_smear = pz_2 + dpz_2;
    double E1_smear = sqrt( pow(px1_smear,2) + pow(py1_smear,2) + pow(pz1_smear,2) + pow(mass_pion,2) );
    double E2_smear = sqrt( pow(px2_smear,2) + pow(py2_smear,2) + pow(pz2_smear,2) + pow(mass_pion,2) );
    
    double qinvSmeared = pow( px1_smear - px2_smear, 2);
    qinvSmeared += pow( py1_smear - py2_smear, 2);
    qinvSmeared += pow( pz1_smear - pz2_smear, 2);
    qinvSmeared -= pow( E1_smear - E2_smear, 2);
    qinvSmeared = sqrt(qinvSmeared);

    for( double fc2 = fc2min; fc2 < fc2max + fc2step/2.0; fc2 = fc2 + fc2step ) {
      for( double R = Rmin; R < Rmax + Rstep/2.0; R = R + Rstep ) {
        double weight_ss = GetFemtoWeight( true, fc2, R, qinvIdeal );
        double weight_os = GetFemtoWeight( false, fc2, R, qinvIdeal );

        ideal_ss->Fill(   fc2, R, qinvIdeal, weight_ss );
        ideal_os->Fill(   fc2, R, qinvIdeal, weight_os );
        ideal_noW->Fill(  fc2, R, qinvIdeal );
        smeared_ss->Fill( fc2, R, qinvSmeared, weight_ss );
        smeared_os->Fill( fc2, R, qinvSmeared, weight_os );
      }
    }

    Qcorrelation->Fill( qinvIdeal, qinvSmeared );

    pairN++;
  } // end pair loop

  // Construct corrections
  TFile *fout = new TFile("MRC_temp.root","RECREATE");

  TH3D *source_ss = (TH3D*)ideal_ss->Clone("source_ss");
  TH3D *source_os = (TH3D*)ideal_os->Clone("source_os");
  TH3D *mrc_ss = (TH3D*)ideal_ss->Clone("mrc_ss");
  TH3D *mrc_os = (TH3D*)ideal_os->Clone("mrc_os");

  for(int x = 1; x <= source_ss->GetNbinsX(); x++ ) {
    for(int y = 1; y <= source_ss->GetNbinsY(); y++ ) {
      for(int z = 1; z <= source_ss->GetNbinsZ(); z++ ) {

        if( ideal_noW->GetBinContent( x, y, z ) > 0 ) {
          double val = source_ss->GetBinContent( x, y, z ) / ideal_noW->GetBinContent( x, y, z );
          double err = source_ss->GetBinError( x, y, z ) / ideal_noW->GetBinContent( x, y, z );
          source_ss->SetBinContent( x, y, z, val );
          source_ss->SetBinError( x, y, z, err );
          //
          val = source_os->GetBinContent( x, y, z ) / ideal_noW->GetBinContent( x, y, z );
          err = source_os->GetBinError( x, y, z ) / ideal_noW->GetBinContent( x, y, z );
          source_os->SetBinContent( x, y, z, val );
          source_os->SetBinError( x, y, z, err );
        }
        if( smeared_ss->GetBinContent( x, y, z ) > 0 ) {
          double val = mrc_ss->GetBinContent( x, y, z ) / smeared_ss->GetBinContent( x, y, z );
          double err = mrc_ss->GetBinError( x, y, z ) / smeared_ss->GetBinContent( x, y, z );
          mrc_ss->SetBinContent( x, y, z, val );
          mrc_ss->SetBinError( x, y, z, err );
          //
          val = mrc_os->GetBinContent( x, y, z ) / smeared_os->GetBinContent( x, y, z );
          err = mrc_os->GetBinError( x, y, z ) / smeared_os->GetBinContent( x, y, z );
          mrc_os->SetBinContent( x, y, z, val );
          mrc_os->SetBinError( x, y, z, err );
        }
      }
    }
  }

  mrc_ss->SetMarkerStyle(20);
  mrc_os->SetMarkerStyle(24);

  source_ss->Write();
  source_os->Write();
  mrc_ss->Write();
  mrc_os->Write();
  Qcorrelation->Write();

  fout->Close();

  seconds = time(0) - seconds;
  cout<<"Minutes = "<<seconds/60.<<endl;
}

//------------------------------------------------------------------------
void LoadFSIfile( TString filename ) {

  TFile *f = new TFile( filename, "READ" );
  mFSI_ss = (TH2D*)f->Get("K2_ss");
  mFSI_os = (TH2D*)f->Get("K2_os");
  mFSI_ss->SetDirectory(0);
  mFSI_os->SetDirectory(0);

  f->Close();
}

//------------------------------------------------------------------------
double GetFemtoWeight( bool sameCharge, double fc2, double Rinv, double Qinv ) {

  double K2 = GetFSI( sameCharge, Rinv, Qinv );
 
  if( sameCharge) {
    return ( (1 - fc2) + fc2*K2*( 1 + exp(-pow(Rinv * Qinv/FmToGeV, 2) ) ) );
  }
  else {
    return ( (1 - fc2) + fc2*K2 );
  }

}

//------------------------------------------------------------------------
double GetFSI( bool sameCharge, double Rinv, double Qinv ) {

  // QbinLow, QbinHigh, RbinLow, RbinHigh
  // |-----c--x--|-----c-----|
  int QbinLow = mFSI_ss->GetXaxis()->FindBin( Qinv );
  int QbinHigh = QbinLow + 1;
  if( Qinv < mFSI_ss->GetXaxis()->GetBinCenter( Qinv ) ) {
    QbinLow--;
    QbinHigh--;
  }
  if( QbinLow < 1 ) { QbinLow = 1; QbinHigh = QbinLow; }
  if( QbinHigh > mFSI_ss->GetNbinsX() ) { QbinHigh = mFSI_ss->GetNbinsX(); QbinLow = QbinHigh; }
  double wQ = (Qinv - mFSI_ss->GetXaxis()->GetBinCenter( QbinLow )) / mFSI_ss->GetXaxis()->GetBinWidth(QbinLow );

  int RbinLow = mFSI_ss->GetYaxis()->FindBin( Rinv );
  int RbinHigh = RbinLow + 1;
  if( Qinv < mFSI_ss->GetYaxis()->GetBinCenter( Rinv ) ) {
    RbinLow--;
    RbinHigh--;
  }
  if( RbinLow < 1 ) { RbinLow = 1; RbinHigh = RbinLow; }
  if( RbinHigh > mFSI_ss->GetNbinsY() ) { RbinHigh = mFSI_ss->GetNbinsY(); RbinLow = RbinHigh; }
  double wR = (Rinv - mFSI_ss->GetYaxis()->GetBinCenter( RbinLow )) / mFSI_ss->GetYaxis()->GetBinWidth(RbinLow );

  TH2D *fsiHist;
  if( sameCharge ) { fsiHist = mFSI_ss; }
  else { fsiHist = mFSI_os; }

  // Bilinear interpolation
  // First along Qinv
  double a1 = fsiHist->GetBinContent( QbinLow, RbinLow ) * (1 - wQ) + fsiHist->GetBinContent( QbinHigh, RbinLow ) * wQ;
  double a2 = fsiHist->GetBinContent( QbinLow, RbinHigh ) * (1 - wQ) + fsiHist->GetBinContent( QbinHigh, RbinHigh ) * wQ;

  // Now along Rinv
  double fsi = a1 * (1 - wR) + a2 * wR;

  return fsi;
}

