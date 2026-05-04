#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <Riostream.h>
#include <complex>

#include "TObject.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
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
#include "TPad.h"
#include "TMinuit.h"
#include "TChain.h"
#include "TLorentzVector.h"

#define PI 3.1415926
#define BohrR 1963.6885 // Mate's value(1963.6885) ~ 387.5/0.197327(1963.7455)
#define FmToGeV 0.197327 // 0.197327
#define fzero_aa 0.942597 // 0.186fm/FmToGeV (scattering length of pi+pi-) = 0.942597
#define fzero_ab -0.89192 // -0.176fm/FmtoGeV (scattering length (pi+pi-)-->(pi0pi0) = -0.89192
#define fzero_bb 0.3119 // fzero_aa + 1/sqrt(2)*fzero_ab = 0.0615/FmToGeV = 0.3119
#define dzero -50.6773 // -10/FmToGeV = -50.6773 (effective range)
#define EulerC 0.5772156649 // Euler constant
#define masspi0 0.1349766 // pi0 mass
#define masspiC 0.1395702 // pi+ mass
#define massKch 0.493677 // K+- mass
#define massK0 0.497614 // K0 mass

#define SF 1.0 // Scale Factor to test Riverside limit
#define RstarMax 405.4/SF // 405.4 (80 fm / FmToGeV), tried 202.7 as a variation
#define RstarMin 0.507/SF // 0.507 (0.1 fm / FmToGeV)
using namespace std;


bool StrongFSI   = true;
int NpairsPerQ   = 10000;
double massPOI   = masspiC;


// g and p used in Gamma function
const int g_gamma = 7;
double p_gamma[9] = {0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
// h coefficients for very small kstar (eta > 0.2)
double hcoef[15]={.8333333333e-1, .8333333333e-2, .396825396825e-2, .4166666667e-2, .7575757576e-2, .2109279609e-1, .8333333333e-1, .4432598039e0 , .305395433e1, .2645621212e2, .2814601449e3, .3607510546e4, .5482758333e5, .9749368235e6, .200526958e8};
// Scattering length (in GeV) coefficients in isospin representation from Lednicky's code from Colangelo (2001)
double fzero[2][6]={{.220, .268, -.0139, -.00139, 36.77, 15*0}, {-.0444, -.0857, -.00221, -.000129, -21.62, 15*0}};


double Gamov(double);
double fact(double);
double hFunction(double);
double KinverseFunction(double, int);
complex<double> gamma(complex<double> z);
complex<double> conf_hyp(complex<double> a, complex<double> b, complex<double> z);


//------------------------------------------------------------------------
void FSI( TString RValue_str = "5.0" ){

  // suffix for outfile name and K2 histo names 
  TString suffix = RValue_str;
  double RValue = RValue_str.Atof();

  TRandom3 *RandomNumber = new TRandom3();
  RandomNumber->SetSeed(0);

  TF1 *SingleParticleGaussian = new TF1("SingleParticleGaussian","exp(-0.5 * pow(x/[0],2) )", -100,100);
  SingleParticleGaussian->FixParameter( 0, RValue );

  TFile *outfile = new TFile(Form("output_fsi_%sfm.root",suffix.Data()), "RECREATE");

  TH1D *rstarPRF = new TH1D("rstarPRF","", 400,0,100);
  //
  int NQbins = 200;
  double QMin = 0.0005;
  double QMax = 0.2005;
  TH1D *Num_BEE_ss=new TH1D("Num_BEE_ss","",NQbins,QMin,QMax);
  TH1D *Num_FSI_ss = new TH1D("Num_FSI_ss","",NQbins,QMin,QMax);
  TH1D *Num_FSI_os = new TH1D("Num_FSI_os","",NQbins,QMin,QMax);
  //
  TH1D *PlaneWF_ss = new TH1D("PlaneWF_ss","",NQbins,QMin,QMax);
  TH1D *PlaneWF_os = new TH1D("PlaneWF_os","",NQbins,QMin,QMax);
  TH1D *Npairs_ss=new TH1D("Npairs_ss","",NQbins,QMin,QMax);
  TH1D *Npairs_os=new TH1D("Npairs_os","",NQbins,QMin,QMax);
  ////////////////////////

  complex<double> binput(1.0, 0.0); // this is fixed ("confluent" hypergeometric function used)
 
  // First Create B(rho,eta) and P(rho,eta) and h(eta) functions (for G hypergeometric functions (strong+Coulomb FSI))
  const int Btermslimit=500;// 500
  double rhoLimit = .05*100/FmToGeV;// 100 fm "limit"
  int rhoBins = rhoLimit/(0.001/FmToGeV);// this binning provides small changes between nearby bins
  double rhoStep = rhoLimit/rhoBins;
  double etaLimit = 1/(0.0005*BohrR);
  int etaBins = 2.0*etaLimit/(0.001);// this binning provides small changes between nearby bins
  double etaStep = (2*etaLimit)/etaBins;
  //cout<<rhoBins<<"  "<<etaBins<<"  "<<rhoLimit/rhoBins<<"  "<<etaLimit/etaBins<<endl;
  TH2D *B_histogram = new TH2D("B_histogram","",rhoBins,0,rhoLimit, etaBins,-etaLimit,etaLimit);
  TH2D *P_histogram = new TH2D("P_histogram","",rhoBins,0,rhoLimit, etaBins,-etaLimit,etaLimit);
  TH1D *h_histogram = new TH1D("h_histogram","",etaBins,-etaLimit,etaLimit);

  for(int i=1; i<=rhoBins; i++){
    for(int j=1; j<=etaBins; j++){
      double rho = (i+0.5)*rhoStep;
      double eta = (j+0.5)*etaStep - etaLimit;// starts from negative values
      if(fabs(eta) < 0.0001) continue;

      double Bfunc[Btermslimit]={0};
      double Pfunc[Btermslimit]={0};
      int StopPoint=Btermslimit;
      Bfunc[0]=1; Bfunc[1]=eta*rho;
      Pfunc[0]=1; Pfunc[1]=0;
      for(int ii=2; ii<Btermslimit; ii++) { 
        Bfunc[ii] = (2*eta*rho*Bfunc[ii-1] - rho*rho*Bfunc[ii-2])/(double(ii*(ii+1.)));
        Pfunc[ii] = (2*eta*rho*Pfunc[ii-1] - rho*rho*Pfunc[ii-2] - (2*(ii-1)+1)*2*eta*rho*Bfunc[ii-1])/(double((ii-1)*ii));
        if(fabs(Bfunc[ii])+fabs(Bfunc[ii-1])+fabs(Pfunc[ii])+fabs(Pfunc[ii-1]) < 1.0e-7) {StopPoint=ii; break;}
      }
      //cout<<StopPoint+1<<endl;
      double Bsum=0, Psum=0;
      for(int ii=0; ii<StopPoint+1; ii++) {Bsum += Bfunc[ii]; Psum += Pfunc[ii];}
      B_histogram->SetBinContent(i,j, Bsum);
      P_histogram->SetBinContent(i,j, Psum);
    }
  }
  cout<<"Done preparing B and P histograms"<<endl;
  for(int j=1; j<=etaBins; j++){
    double eta = (j+0.5)*etaStep - etaLimit;// starts from negative values
    if(fabs(eta) < 0.0001) continue;
    h_histogram->SetBinContent(j, hFunction(eta));
  }
  cout<<"Done preparing h histogram"<<endl;

  TVector3 P1_PRF;
  TVector3 P2_PRF;
  TVector3 X1_PRF;
  TVector3 X2_PRF;

  int seconds = time(0);

  double step_kstar = 0.0005;

  for( double kstar = 0.001; kstar < (QMax+step_kstar)/2.0; kstar = kstar + step_kstar ) {
    double E = sqrt( pow(kstar, 2) + pow(massPOI, 2) );
    P1_PRF[0] = kstar;
    P1_PRF[1] = 0;
    P1_PRF[2] = 0;
    P2_PRF[0] = -kstar;
    P2_PRF[1] = 0;
    P2_PRF[2] = 0;
    double qinv12 = (P1_PRF - P2_PRF).Mag();
    //cout<<"qinv: "<<qinv12<<endl;

    for( int pair = 0; pair < NpairsPerQ; pair++ ) {

      for( int i = 0; i< 3; i++) {
        X1_PRF[i] = SingleParticleGaussian->GetRandom() / FmToGeV / SF;
        X2_PRF[i] = SingleParticleGaussian->GetRandom() / FmToGeV / SF;
      }
         
      double rstar12_PRF = (X1_PRF - X2_PRF).Mag(); 
      if( rstar12_PRF < RstarMin ) { continue; }

      double rho = rstar12_PRF * qinv12/2.;
      double phase = (P1_PRF-P2_PRF)*(X1_PRF-X2_PRF);
      double Cosphase = phase/(rstar12_PRF)/(qinv12);

      // ++: chargeProd = 0    +-: chargeProd = 1
      for( int chargeProd = 0; chargeProd < 2; chargeProd++ ) { 
        bool sameCharge = (chargeProd == 0) ? true : false;

        double eta = 1/(qinv12/2.*BohrR);
        if( ! sameCharge ) {
          eta = -eta;
        }

        complex<double> ainput(0.0, -eta);
        complex<double> zinput12(0.0, rho + rho*Cosphase);
        complex<double> zinput21(0.0, rho - rho*Cosphase);
        complex<double> CHG12(1.0, 0);
        complex<double> CHG21(1.0, 0);
        double GamovFactor=1.0;
        double BEE = 0.;
        double WeightWF = 1.0;
        double WeightPlaneWF = 1.0;

        // Strong FSI for +- case
        complex<double> Psi_alpha(0., 0.);
        if( StrongFSI && ! sameCharge ) {
          double kb = sqrt(pow(massPOI,2)+pow(qinv12/2.,2) - pow(masspi0,2));
          double G2 = KinverseFunction(pow(qinv12/2.,2),1);
          double G1 = KinverseFunction(pow(qinv12/2.,2),0);
          double RK11 = 2/3.*G1 + 1/3.*G2;
          double RK22 = 2/3.*G2 + 1/3.*G1;
          double RK12 = -sqrt(2/9.)*(G1-G2);
          complex<double> Chi(h_histogram->GetBinContent(h_histogram->GetXaxis()->FindBin(eta)), Gamov(eta)/(2*eta));
          complex<double> C3 = RK11 - 2.0*Chi/BohrR;
          complex<double> C5(RK22, -kb);
          complex<double> C10 = C3*C5 - pow(RK12,2);
          complex<double> fc_aa = C5/C10;
          complex<double> G = P_histogram->GetBinContent(P_histogram->GetXaxis()->FindBin(rho), P_histogram->GetYaxis()->FindBin(eta));
          G += 2.0*eta*rho*( log(fabs(2.0*eta*rho)) + 2.0*EulerC - 1.0 + Chi )*B_histogram->GetBinContent(B_histogram->GetXaxis()->FindBin(rho), B_histogram->GetYaxis()->FindBin(eta));
          Psi_alpha = fc_aa*G/(rstar12_PRF);
        }

        CHG12= conf_hyp( ainput, binput, zinput12 );
        CHG21= conf_hyp( ainput, binput, zinput21 );
        GamovFactor = Gamov( eta );
        BEE = cos(phase);

        complex<double> planewave12(double(TMath::Cos(-rho*Cosphase)), double(TMath::Sin(-rho*Cosphase)));
        complex<double> planewave21(planewave12.real(), -planewave12.imag());
        //
        complex<double> Full12(0.0,0.0);
        Full12 += planewave12;
        Full12 *= CHG12;
        if( chargeProd == 1 ) { Full12 += Psi_alpha; }
        //
        complex<double> Full21(0.0,0.0);
        Full21 += planewave21;
        Full21 *= CHG21;
        //
        complex<double> FullWF = Full12;
        complex<double> FullPlaneWF = planewave12;
        if( sameCharge ) {
          FullWF += Full21;
          FullWF /= sqrt(2.);
          FullPlaneWF += planewave21;
          FullPlaneWF /= sqrt(2.);
        }

        WeightWF = GamovFactor*pow(abs(FullWF),2);
        WeightPlaneWF = pow(abs(FullPlaneWF),2);

        if( sameCharge ) { // same-charge
          rstarPRF->Fill( rstar12_PRF*FmToGeV );
          Num_BEE_ss->Fill( qinv12, BEE );
          Num_FSI_ss->Fill( qinv12, WeightWF );
          PlaneWF_ss->Fill( qinv12, WeightPlaneWF );
          Npairs_ss->Fill( qinv12 );
        }
        else { // mixed-charge
          Num_FSI_os->Fill( qinv12, WeightWF );
          PlaneWF_os->Fill( qinv12, WeightPlaneWF );
          Npairs_os->Fill( qinv12 );
        }
      } // chargeProd
    }// pair

  }// kstar

  outfile->cd();
  //
  rstarPRF->Write();
  Num_BEE_ss->Write();
  Num_FSI_ss->Write();
  Num_FSI_os->Write();
  //
  PlaneWF_ss->Write();
  PlaneWF_os->Write();
  Npairs_ss->Write();
  Npairs_os->Write();
  //
  TH1D *K2_ss = (TH1D*)Num_FSI_ss->Clone();
  K2_ss->SetName("K2_ss");
  K2_ss->SetTitle(Form("Gaussian source R=%s fm",suffix.Data()));
  K2_ss->Divide(PlaneWF_ss);
  K2_ss->SetMarkerStyle(20);
  K2_ss->SetMarkerColor(2);
  TH1D *K2_os = (TH1D*)Num_FSI_os->Clone();
  K2_os->SetName("K2_os");
  K2_os->SetTitle(Form("Gaussian source R=%s fm",suffix.Data()));
  K2_os->Divide(PlaneWF_os);
  K2_os->SetMarkerStyle(20);
  K2_os->SetMarkerColor(4);
  
  K2_ss->Draw("hist p");
  K2_os->Draw("hist p same");
  K2_ss->SetMaximum(2.2);
  gPad->SetGridx();
  gPad->SetGridy();

  K2_ss->Write();
  K2_os->Write();

  cout<<setprecision(8)<<"Qinv: "<<K2_ss->GetXaxis()->GetBinCenter(60)<<"  K2_ss:"<<K2_ss->GetBinContent(60)<<endl;
  //outfile->Close();

  //seconds = time(0) - seconds;
  //cout<<"Minutes = "<<seconds/60.<<endl;

}

//------------------------------------------------------------------------
double Gamov(double eta){
  return (2*PI*eta/(exp(2*PI*eta)-1));
}

//------------------------------------------------------------------------
double fact(double n){
  return (n < 1.00001) ? 1 : fact(n - 1) * n;
}

//------------------------------------------------------------------------
// h function from Lednicky, phys. of part. and nuclei 40:307 (2009)
double hFunction(double eta){
  double DS=1.0;
  double N=0;
  double S=0;
  double returnValue=0;
  if(eta < 3.0){
    while(DS > 1.0e-13){
      N++;
      DS = 1/(N*(pow(N/eta,2)+1.0));
      S += DS;
    }
    returnValue = (S - EulerC - log(fabs(eta)));
  }else{// small kstar, high eta
    double etaSquared=pow(eta,2);
    double etaPowered=etaSquared;
    for(int ii=0; ii<9; ii++){// 9 was maximum value for Lednicky's code
      returnValue += 1/etaPowered * hcoef[ii];
      etaPowered *= etaSquared;
    }
  }
  return returnValue; 
}

//------------------------------------------------------------------------
double KinverseFunction(double kstarSq, int J){
  double ESq = kstarSq + pow(massPOI,2);
  double E = sqrt(ESq);
  double xSq = kstarSq/pow(massPOI,2);
  double Kinverse = E*(4*ESq-fzero[J][4])/(4*pow(massPOI,2)-fzero[J][4]);
  Kinverse /= fzero[J][0] + fzero[J][1]*xSq + fzero[J][2]*pow(xSq,2) + fzero[J][3]*pow(xSq,3);
  return Kinverse;
}

//------------------------------------------------------------------------
complex<double> gamma(complex<double> z)
{// Lanczos approximation

  if(floor(abs(z)) == z)
    return fact(abs(z)-1.0);

  if(real(z)<(double)(0.5)) return PI/(sin(PI*z)*gamma((1.0-z)));
  z -= 1.0;
  complex<double> x = p_gamma[0];
  for(int i=1;i<g_gamma+2;i++) x = x+p_gamma[i]/(z+(double)i);
  complex<double> t = z+(double)g_gamma+0.5;
  return sqrt(2*PI)*pow(t,z+0.5)*exp(-t)*x;
}

//------------------------------------------------------------------------
// confluent hypergeometric function
complex<double> conf_hyp(complex<double> a, complex<double> b, complex<double> z)
{
  // b is not used (confluent)
  complex<double> S(1.0, 0.0);
  complex<double> F=S;
  complex<double> alf1 = a - 1.0;
  double J=0;
  double errorCutOff=0.00001;// series truncation point (0.00001, 1e-6 for strict case)
  if(abs(z) > 40.) {// high k*r* case (was 20., now 40. Improved accuracy for qinv>60 MeV)
    double eta = -a.imag();
    double RKS = z.imag();
    double D1 = log(RKS)*eta;
    double D2 = pow(eta,2)/RKS;
    complex<double> arg(1.0, eta);
    complex<double> EIDC = gamma(arg);
    EIDC /= abs(EIDC);
    complex<double> ZZ1(cos(D1), sin(D1));
    ZZ1 /= EIDC;
    complex<double> value(1.0, -D2);
    value *= ZZ1;
    complex<double> value2(cos(RKS), sin(RKS));
    F = value - value2*eta/RKS/ZZ1;
    return F;
  }else {// low k*r* case
    while(fabs(S.real())+fabs(S.imag()) > errorCutOff){
      J++;
      complex<double> A = (alf1 + J)/(pow(J,2));
      S *= A*z;
      F += S;
    }
    return F;
  }
}

