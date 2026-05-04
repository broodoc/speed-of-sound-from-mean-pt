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
#include "TGraph.h"
#include "TSpline.h"
#include "TVirtualFitter.h"

#define BohrR 1963.6885 // Bohr Radius for pions
#define FmToGeV 0.19733 // conversion to fm

using namespace std;

enum KtBins {
  KTBIN1,
  KTBIN2,
  KTBIN3,
  NKTBINS
};

map<int,string> ccToString = { {0,"--"}, {1,"++"}, {2,"-+"} };

const int Mbins = 4;
const int EDbins = 1;
const double QFitMin = 0.005;
const double QFitMax = 0.1;
const double QMax_mrc = 0.07;

TH3D *gMRC_ss;
TH3D *gMRC_os;
TH2D *gFSI_ss;
TH2D *gFSI_os;


vector<TH1D*> gC2raw;
vector<TH1D*> gC2corr;
vector<TF1*> gC2fit;

TH1D *gC2;

double gRinv;
double gChi2 = 0;;
double gNFitPoints;
bool gSameCharge;

double GetFSI( bool sameCharge, double Rinv, double Qinv );
void fcn_c2qinv(int&, double*, double&, double[], int);
void ApplyMRC( bool sameCharge, TH1D *raw, TH1D *corr, double fc2, double Rinv );
bool LoadMRC( TString filename );
bool LoadFSIfile( TString filename );
Double_t fitFcn_forTF1(Double_t *x, Double_t *par);



//--------------------------------------------------------------------------------
void fitCorrelations( string fname = "rawCorrelations.root" ) 
{

  if( ! LoadMRC( "../corrections/MRC_3p0_2018.root" ) ) {
    cout<<"No MRC file"<<endl;
    return;
  }

  if( ! LoadFSIfile( "../corrections/KFile.root" ) ) {
    cout<<"No FSI file"<<endl;
    return;
  }

  const int npar_fit = 3;
  string parName[ npar_fit ] = {"N", "fc2", "Rinv"};
  double par_def[ npar_fit ] = {1.0, 0.5, 4.5};
  double par_steps[ npar_fit ] = {0.01, 0.01, 0.1};
  double par_mins[ npar_fit ] = {0.8, 0.3, 1.0};
  double par_maxs[ npar_fit ] = {1.1, 1.0, 10.0};

  double fitparams[ npar_fit ] = {0};
  double fitparams_e[ npar_fit ] = {0};
  
  TMinuit minuitFit_c2qinv( npar_fit );
  minuitFit_c2qinv.SetFCN( fcn_c2qinv );

  int ierflg = 0;
  double arglist = 0;
  minuitFit_c2qinv.mnexcm( "SET NOWarnings", &arglist, 1, ierflg );  
  arglist = -1;
  minuitFit_c2qinv.mnexcm( "SET PRint", &arglist, 1, ierflg );
  //arglist=2;// improve Minimization Strategy (1 is default)
  //minuitFit_c2qinv.mnexcm("SET STR",&arglist,1,ierflg);
  //arglist = 0;
  //minuitFit_c2qinv.mnexcm("SCAN", &arglist,1,ierflg);
  arglist = 5000;
  minuitFit_c2qinv.mnexcm( "MIGRAD", &arglist , 1, ierflg );

  
  TFile *f = new TFile( fname.data(), "READ" );

  
  for( int m = 0; m < Mbins; m++ ) {
    for( int ed = 0; ed < EDbins; ed++ ) {

      for( int kt = 0 ; kt < KtBins::NKTBINS; kt++ ) {

        double Rtemp[3] = {0}; // --, ++, -+

        for( int cc = 0; cc < 3; cc++ ) {

          string chargeComb;
          if( cc == 0 ) { chargeComb = "00"; }
          else if ( cc == 1) { chargeComb = "11"; }
          else { chargeComb = "01"; }

          string nameRaw = "C2_" + chargeComb + "_ed" + to_string(ed) + "_kT" + to_string(kt) + "_M" + to_string(m);
          TH1D *c2raw = (TH1D*)f->Get( nameRaw.data() );
          
          string nameC2corr = "C2_" + chargeComb + "_ed" + to_string(ed) + "_kT" + to_string(kt) + "_M" + to_string(m);
          gC2 = (TH1D*)c2raw->Clone( nameC2corr.data() );
          
          string rawName = "raw" + nameC2corr;
          c2raw->SetName( rawName.data() );
          string rawTitle = "raw " + string( c2raw->GetTitle() );
          c2raw->SetTitle( rawTitle.data() );

          // Initialize fit 
          for( int i = 0; i < npar_fit; i++ ) {
            minuitFit_c2qinv.DefineParameter( i, parName[i].c_str(), par_def[i], par_steps[i], par_mins[i], par_maxs[i] );
            minuitFit_c2qinv.GetParameter( i, fitparams[i], fitparams_e[i] );
          }

          bool satisfactoryFit = false;

          int ierr;
          if( cc != 2 ) { 
            gSameCharge = true;
            minuitFit_c2qinv.mnfree(1);
            gRinv = fitparams[2];
          }
          else {
            // Fix Rinv for os case from ss fit
            gRinv = (Rtemp[0] + Rtemp[1])/2.0;
            minuitFit_c2qinv.DefineParameter( 2, parName[2].c_str(), gRinv, par_steps[2], par_mins[2], par_maxs[2] );
            minuitFit_c2qinv.mnfixp(2, ierr); 
            gSameCharge = false;
          }
          
          //if( m != 0 || kt != 0 ) { continue; }
         
          ////////////////////////////
          while( ! satisfactoryFit ) {

            ApplyMRC( gSameCharge, c2raw, gC2, fitparams[1], gRinv );

            //----------------------------
            // Do the minimization!
            cout<<"Start fit M:"<<m<<"  ed:"<<ed<<"  cc:"<<ccToString[cc]<<"  kt:"<<kt<<endl;
            minuitFit_c2qinv.Migrad();// Minuit's best minimization algorithm
            minuitFit_c2qinv.mnexcm( "SHOw PARameters", &arglist, 1, ierflg );
            double chi2NDF = gChi2 / (gNFitPoints - minuitFit_c2qinv.GetNumFreePars());
            cout<<" Chi2/NDF = "<<chi2NDF<<endl;

            for( int i = 0; i < npar_fit; i++ ) {
              minuitFit_c2qinv.GetParameter( i, fitparams[i], fitparams_e[i] );
            }

            if( cc != 2 && fabs(gRinv - fitparams[2]) < 0.5 ) {
              satisfactoryFit = true;
            }
            if( cc == 2 ) { satisfactoryFit = true; }

            gRinv = fitparams[2];
          }
          Rtemp[ cc ] = fitparams[2];
          ////////////////////////////
          
          TH1D *c2corr = (TH1D*)gC2->Clone();

          TF1 *fit_temp = new TF1( "fit_temp", fitFcn_forTF1, 0.001,0.2, 3);
          fit_temp->SetParameter( 0, fitparams[0] );
          fit_temp->SetParameter( 1, fitparams[1] );
          fit_temp->SetParameter( 2, fitparams[2] );

          string fitname = "fit_" + chargeComb + "_ed" + to_string(ed) + "_kT" + to_string(kt) + "_M" + to_string(m);
          TF1 *fit = (TF1*)fit_temp->Clone( fitname.data() );
          fit->SetTitle( fitname.data() );
          
          // color on font
          if( gSameCharge ) {
            if( cc == 0 ) { 
              c2corr->SetMarkerStyle(24); 
              c2corr->SetMarkerColor(2);
              c2corr->SetLineColor(2);
            }
            else { 
              c2corr->SetMarkerStyle(25); 
              c2corr->SetMarkerColor(2);
              c2corr->SetLineColor(2);
              fit->SetLineStyle(2);
            }

            fit->SetLineColor(2);
          }
          else {
            c2corr->SetMarkerStyle(24);
            c2corr->SetMarkerColor(4);
            c2corr->SetLineColor(4);
            fit->SetLineColor(4);
          }

          gC2corr.push_back( c2corr );
          gC2raw.push_back( c2raw );
          gC2fit.push_back( fit );

        }
      }
    }
  }

  TH2D *FitSummary_R_00 = new TH2D("FitSummary_R_00","R_{inv} --;centrality;k_{T} (GeV/c);R_{inv} (fm)",Mbins,0.5,Mbins+0.5, KtBins::NKTBINS,0.5,KtBins::NKTBINS+0.5);
  TH2D *FitSummary_R_11 = new TH2D("FitSummary_R_11","R_{inv} --;centrality;k_{T} (GeV/c);R_{inv} (fm)",Mbins,0.5,Mbins+0.5, KtBins::NKTBINS,0.5,KtBins::NKTBINS+0.5);
  TH2D *FitSummary_lambda_00 = new TH2D("FitSummary_lambda_00","#lambda_{inv} --;centrality;k_{T} (GeV/c);#lambda_{inv}",Mbins,0.5,Mbins+0.5, KtBins::NKTBINS,0.5,KtBins::NKTBINS+0.5);
  TH2D *FitSummary_lambda_11 = new TH2D("FitSummary_lambda_11","#lambda_{inv} --;centrality;k_{T} (GeV/c);#lambda_{inv}",Mbins,0.5,Mbins+0.5, KtBins::NKTBINS,0.5,KtBins::NKTBINS+0.5);
  TH2D *FitSummary_fc2_01 = new TH2D("FitSummary_fc2_01","f_{c}^{2} --;centrality;k_{T} (GeV/c);f_{c}^{2}",Mbins,0.5,Mbins+0.5, KtBins::NKTBINS,0.5,KtBins::NKTBINS+0.5);

  for( uint i = 0; i < gC2fit.size(); i++ ) {

    TString title = gC2fit[ i ]->GetTitle();
    TString ktpart = title(title.Index("kT")+2, 1);
    int kTbin = ktpart.Atoi() + 1;
    TString Mpart = title(title.Index("_M")+2, 1);
    int Mbin = Mpart.Atoi() + 1;
    //cout<<title<<"  "<<kTbin<<"  "<<Mbin<<endl;
    if( title.Contains("00") ) {
      FitSummary_lambda_00->SetBinContent( Mbin, kTbin, gC2fit[ i ]->GetParameter(1) );
      FitSummary_lambda_00->SetBinError( Mbin, kTbin, gC2fit[ i ]->GetParError(1) );
      FitSummary_R_00->SetBinContent( Mbin, kTbin, gC2fit[ i ]->GetParameter(2) );
      FitSummary_R_00->SetBinError( Mbin, kTbin, gC2fit[ i ]->GetParError(2) );
    }
    else if( title.Contains("11") ) {
      FitSummary_lambda_11->SetBinContent( Mbin, kTbin, gC2fit[ i ]->GetParameter(1) );
      FitSummary_lambda_11->SetBinError( Mbin, kTbin, gC2fit[ i ]->GetParError(1) );
      FitSummary_R_11->SetBinContent( Mbin, kTbin, gC2fit[ i ]->GetParameter(2) );
      FitSummary_R_11->SetBinError( Mbin, kTbin, gC2fit[ i ]->GetParError(2) );
    }
    else {
      FitSummary_fc2_01->SetBinContent( Mbin, kTbin, gC2fit[ i ]->GetParameter(1) );
      FitSummary_fc2_01->SetBinError( Mbin, kTbin, gC2fit[ i ]->GetParError(1) );
    }
  }

  TFile *fout = new TFile( "results.root", "RECREATE" );

  for( uint i = 0; i < gC2corr.size(); i++ ) {
    gC2raw[i]->Write();
    gC2corr[i]->Write();
    gC2fit[i]->Write();
  }

  FitSummary_lambda_00->Write();
  FitSummary_lambda_11->Write();
  FitSummary_fc2_01->Write();
  FitSummary_R_00->Write();
  FitSummary_R_11->Write();

}

//------------------------------------------------------------------------
double GetFSI( bool sameCharge, double Rinv, double Qinv ) {

  // QbinLow, QbinHigh, RbinLow, RbinHigh
  // |-----c--x--|-----c-----|
  int QbinLow = gFSI_ss->GetXaxis()->FindBin( Qinv );
  int QbinHigh = QbinLow + 1;
  if( Qinv < gFSI_ss->GetXaxis()->GetBinCenter( Qinv ) ) {
    QbinLow--;
    QbinHigh--;
  }
  if( QbinLow < 1 ) { QbinLow = 1; QbinHigh = QbinLow; }
  if( QbinHigh > gFSI_ss->GetNbinsX() ) { QbinHigh = gFSI_ss->GetNbinsX(); QbinLow = QbinHigh; }
  double wQ = (Qinv - gFSI_ss->GetXaxis()->GetBinCenter( QbinLow )) / gFSI_ss->GetXaxis()->GetBinWidth(QbinLow );

  int RbinLow = gFSI_ss->GetYaxis()->FindBin( Rinv );
  int RbinHigh = RbinLow + 1;
  if( Qinv < gFSI_ss->GetYaxis()->GetBinCenter( Rinv ) ) {
    RbinLow--;
    RbinHigh--;
  }
  if( RbinLow < 1 ) { RbinLow = 1; RbinHigh = RbinLow; }
  if( RbinHigh > gFSI_ss->GetNbinsY() ) { RbinHigh = gFSI_ss->GetNbinsY(); RbinLow = RbinHigh; }
  double wR = (Rinv - gFSI_ss->GetYaxis()->GetBinCenter( RbinLow )) / gFSI_ss->GetYaxis()->GetBinWidth(RbinLow );

  TH2D *fsiHist;
  if( sameCharge ) { fsiHist = gFSI_ss; }
  else { fsiHist = gFSI_os; }

  // Bilinear interpolation
  // First along Qinv
  double a1 = fsiHist->GetBinContent( QbinLow, RbinLow ) * (1 - wQ) + fsiHist->GetBinContent( QbinHigh, RbinLow ) * wQ;
  double a2 = fsiHist->GetBinContent( QbinLow, RbinHigh ) * (1 - wQ) + fsiHist->GetBinContent( QbinHigh, RbinHigh ) * wQ;

  // Now along Rinv
  double fsi = a1 * (1 - wR) + a2 * wR;

  return fsi;
}

//--------------------------------------------------------------------------------
void fcn_c2qinv( int& npar, double* deriv, double& f, double par[], int flag ) {

  gChi2 = 0;
  gNFitPoints = 0;

  TH2D *K2 = gFSI_ss;
  if( ! gSameCharge ) { K2 = gFSI_os; }

  for( int qbin = 1; qbin <= gC2->GetNbinsX(); qbin++ ) {

    double qinv = gC2->GetXaxis()->GetBinCenter( qbin );

    if( qinv < QFitMin ) { continue; }
    if( qinv > QFitMax ) { continue; }

    double K2 = GetFSI( gSameCharge, gRinv, qinv );
    double val = gC2->GetBinContent( qbin );
    double val_e = gC2->GetBinError( qbin );
    //cout<<qinv<<"  "<<val_e<<endl;
    
    double func = 0;
    double err = 0;
    if( gSameCharge ) {
      func = par[0] * ( (1 - par[1]) + par[1]*K2*(1 + exp(-pow(qinv*par[2]/FmToGeV,2)) ) );
      err = val_e;
    }
    else {
      func = par[0] * ( (1 - par[1]) + par[1]*K2 );
      err = val_e;
    }
    
    gChi2 += pow( (val - func) / val_e, 2);
    gNFitPoints++;
  }

  f = gChi2;// Chi2 minimization
}

//--------------------------------------------------------------------------------
double fitFcn_forTF1(double *x, double *par)
{
   double qinv =x[0];
   double f = 0;
   if( gSameCharge ) {
     f = par[0] * ( (1 - par[1]) + par[1]*GetFSI( gSameCharge, par[2], qinv )*(1 + exp(-pow(qinv*par[2]/FmToGeV,2)) ) );
   }
   else {
     f = par[0] * ( (1 - par[1]) + par[1]*GetFSI( gSameCharge, gRinv, qinv ) );
   }
   return f;
}


//--------------------------------------------------------------------------------
void ApplyMRC( bool sameCharge, TH1D *raw, TH1D *corr, double fc2, double Rinv ) {

  TH3D *mrc = gMRC_ss;
  if( ! sameCharge ) { mrc = gMRC_os; }
    
  int fc2Bin = mrc->GetXaxis()->FindBin( fc2 );
  int RinvBin = mrc->GetYaxis()->FindBin( Rinv );

  for( int qbin = 1; qbin <= corr->GetNbinsX(); qbin++ ) {
    if( corr->GetXaxis()->GetBinUpEdge( qbin ) > QMax_mrc ) { continue; }

    corr->SetBinContent( qbin, raw->GetBinContent( qbin ) * mrc->GetBinContent( fc2Bin, RinvBin, qbin ) );
    corr->SetBinError( qbin, raw->GetBinError( qbin ) * mrc->GetBinContent( fc2Bin, RinvBin, qbin ) );
  }

}

//--------------------------------------------------------------------------------
bool LoadMRC( TString filename ) {

  TFile *f = new TFile( filename, "READ" );
  gMRC_ss = (TH3D*)f->Get( "mrc_ss" );
  gMRC_os = (TH3D*)f->Get( "mrc_os" );
  gMRC_ss->SetDirectory(0);
  gMRC_os->SetDirectory(0);

  f->Close();

  if( gMRC_ss && gMRC_os ) { return true; }
  else { return false; }
}

//--------------------------------------------------------------------------------
bool LoadFSIfile( TString filename ) {

  TFile *f = new TFile( filename, "READ" );
  gFSI_ss = (TH2D*)f->Get( "K2_ss" );
  gFSI_os = (TH2D*)f->Get( "K2_os" );
  gFSI_ss->SetDirectory(0);
  gFSI_os->SetDirectory(0);

  f->Close();

  if( gFSI_ss && gFSI_os ) { return true; }
  else { return false; }
}

