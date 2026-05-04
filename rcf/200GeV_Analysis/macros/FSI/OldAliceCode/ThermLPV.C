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
#include "TMinuit.h"
#include "TChain.h"

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
using namespace std;

// g and p used in Gamma function
const int g_gamma = 7;
double p_gamma[9] = {0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7};
// h coefficients for very small kstar (eta > 0.2)
double hcoef[15]={.8333333333e-1, .8333333333e-2, .396825396825e-2, .4166666667e-2, .7575757576e-2, .2109279609e-1, .8333333333e-1, .4432598039e0 , .305395433e1, .2645621212e2, .2814601449e3, .3607510546e4, .5482758333e5, .9749368235e6, .200526958e8};
// Scattering length (in GeV) coefficients in isospin representation from Lednicky's code from Colangelo (2001)
double fzero[2][6]={{.220, .268,-.0139,-.00139, 36.77, 15*0}, {-.0444,-.0857,-.00221,-.000129,-21.62, 15*0}};


float Gamov(float);
float fact(float);
double hFunction(double);
double KinverseFunction(double, int);
void Shuffle(int*, int, int);
void BoostPRF(float [4], float [4], float [4], float*);


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

// confluent hypergeometric function
complex<double> conf_hyp(complex<double> a, complex<double> b, complex<double> z)
{// b is assumed to be 1.0 here
  complex<double> S(1.0, 0.0);
  complex<double> F=S;
  complex<double> alf1 = a - 1.0;
  double J=0;
  double errorCutOff=0.00001;// series truncation point
  if(abs(z) > 20.) {// high k*r* case
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

void ThermLPV(){

  TRandom3 *RandomNumber = new TRandom3();
  RandomNumber->SetSeed(0);
  
  int Nfiles=4395;// 4395 for 60-70%, 1999 max for 40-50 and 50-60, 99 for 30-40, 19 for 0-5 and 5-10 
  bool IncludeFSI=kFALSE;

  TFile *outfile = new TFile("Therm_dist_temp.root","RECREATE");
  TH1D *r_dist=new TH1D("r_dist","",400,0,100);
  TH1D *x_dist=new TH1D("x_dist","",400,0,100);
  TH1D *Pt_dist=new TH1D("Pt_dist","",200,0,2);
  //
  TH2D *rstarPRF_dist_ss=new TH2D("rstarPRF_dist_ss","",20,0,1.0, 400,0,100);
  TH2D *tstar_dist_ss=new TH2D("tstar_dist_ss","",20,0,1.0, 400,0,100);
  //
  TH2D *rstarPRF_dist_os=new TH2D("rstarPRF_dist_os","",20,0,1.0, 400,0,100);
  TH2D *tstar_dist_os=new TH2D("tstar_dist_os","",20,0,1.0, 400,0,100);

  TH2D *Num_Cos_ss=new TH2D("Num_Cos_ss","",20,0,1, 40,0,0.2);
  TH2D *Num_Cos_os=new TH2D("Num_Cos_os","",20,0,1, 40,0,0.2);
  TH2D *Num_CosFSI_ss=new TH2D("Num_CosFSI_ss","",20,0,1, 40,0,0.2);
  TH2D *Num_CosFSI_os=new TH2D("Num_CosFSI_os","",20,0,1, 40,0,0.2);
  Num_CosFSI_ss->Sumw2();
  Num_CosFSI_os->Sumw2();
  TH2D *Num_PrimCosFSI_ss=new TH2D("Num_PrimCosFSI_ss","",20,0,1, 40,0,0.2);
  TH2D *Num_PrimCosFSI_os=new TH2D("Num_PrimCosFSI_os","",20,0,1, 40,0,0.2);
  //
  TH2D *Den_ss=new TH2D("Den_ss","",20,0,1, 40,0,0.2);
  TH2D *Den_os=new TH2D("Den_os","",20,0,1, 40,0,0.2);
  TH2D *DenEM_ss=new TH2D("DenEM_ss","",20,0,1, 40,0,0.2);
  TH2D *DenEM_os=new TH2D("DenEM_os","",20,0,1, 40,0,0.2);
  
  //
  TH3D *NumLambda_ss=new TH3D("NumLambda_ss","",20,0,1, 40,0,0.2, 100,0,100);
  TH3D *NumLambda_os=new TH3D("NumLambda_os","",20,0,1, 40,0,0.2, 100,0,100);
  TH2D *DenLambda_ss=new TH2D("DenLambda_ss","",20,0,1, 40,0,0.2);
  TH2D *DenLambda_os=new TH2D("DenLambda_os","",20,0,1, 40,0,0.2);
  //
  
 
  ////////////////////////
  // LPV studies
  TProfile *ThreePoint_ss = new TProfile("ThreePoint_ss","",1,0.5,1.5, -10.,10.,"");
  TProfile *ThreePoint_os = new TProfile("ThreePoint_os","",1,0.5,1.5, -10.,10.,"");
  TProfile *aEta_coscos_ss = new TProfile("aEta_coscos_ss","",10,0,1, -10,10,"");
  TProfile *aEta_sinsin_ss = new TProfile("aEta_sinsin_ss","",10,0,1, -10,10,"");
  TProfile *aEta_coscos_os = new TProfile("aEta_coscos_os","",10,0,1, -10,10,"");
  TProfile *aEta_sinsin_os = new TProfile("aEta_sinsin_os","",10,0,1, -10,10,"");
  //
  TProfile *dEta_coscos_ss = new TProfile("dEta_coscos_ss","",20,0,2, -10,10,"");
  TProfile *dEta_sinsin_ss = new TProfile("dEta_sinsin_ss","",20,0,2, -10,10,"");
  TProfile *dEta_coscos_os = new TProfile("dEta_coscos_os","",20,0,2, -10,10,"");
  TProfile *dEta_sinsin_os = new TProfile("dEta_sinsin_os","",20,0,2, -10,10,"");
  //
  TProfile *aPt_coscos_ss = new TProfile("aPt_coscos_ss","",20,0.15,2, -10,10,"");
  TProfile *aPt_sinsin_ss = new TProfile("aPt_sinsin_ss","",20,0.15,2, -10,10,"");
  TProfile *aPt_coscos_os = new TProfile("aPt_coscos_os","",20,0.15,2, -10,10,"");
  TProfile *aPt_sinsin_os = new TProfile("aPt_sinsin_os","",20,0.15,2, -10,10,"");
  //
  TProfile *dPt_coscos_ss = new TProfile("dPt_coscos_ss","",20,0,2, -10,10,"");
  TProfile *dPt_sinsin_ss = new TProfile("dPt_sinsin_ss","",20,0,2, -10,10,"");
  TProfile *dPt_coscos_os = new TProfile("dPt_coscos_os","",20,0,2, -10,10,"");
  TProfile *dPt_sinsin_os = new TProfile("dPt_sinsin_os","",20,0,2, -10,10,"");
  //
  aEta_coscos_ss->SetMarkerStyle(21); aEta_coscos_ss->SetMarkerColor(2); aEta_coscos_ss->SetLineColor(2); 
  aEta_sinsin_ss->SetMarkerStyle(25); aEta_sinsin_ss->SetMarkerColor(2); aEta_sinsin_ss->SetLineColor(2); 
  aEta_coscos_os->SetMarkerStyle(20); aEta_coscos_os->SetMarkerColor(4); aEta_coscos_os->SetLineColor(4); 
  aEta_sinsin_os->SetMarkerStyle(24); aEta_sinsin_os->SetMarkerColor(4); aEta_sinsin_os->SetLineColor(4); 
  aEta_coscos_ss->GetXaxis()->SetTitle("|#LT#eta#GT|"); aEta_coscos_ss->GetYaxis()->SetTitle("#LTcos()cos()#GT or #LTsin()sin()#GT"); 
  //
  dEta_coscos_ss->SetMarkerStyle(21); dEta_coscos_ss->SetMarkerColor(2); dEta_coscos_ss->SetLineColor(2); 
  dEta_sinsin_ss->SetMarkerStyle(25); dEta_sinsin_ss->SetMarkerColor(2); dEta_sinsin_ss->SetLineColor(2); 
  dEta_coscos_os->SetMarkerStyle(20); dEta_coscos_os->SetMarkerColor(4); dEta_coscos_os->SetLineColor(4); 
  dEta_sinsin_os->SetMarkerStyle(24); dEta_sinsin_os->SetMarkerColor(4); dEta_sinsin_os->SetLineColor(4); 
  dEta_coscos_ss->GetXaxis()->SetTitle("|#Delta#eta|"); dEta_coscos_ss->GetYaxis()->SetTitle("#LTcos()cos()#GT or #LTsin()sin()#GT");
  //
  aPt_coscos_ss->SetMarkerStyle(21); aPt_coscos_ss->SetMarkerColor(2); aPt_coscos_ss->SetLineColor(2); 
  aPt_sinsin_ss->SetMarkerStyle(25); aPt_sinsin_ss->SetMarkerColor(2); aPt_sinsin_ss->SetLineColor(2); 
  aPt_coscos_os->SetMarkerStyle(20); aPt_coscos_os->SetMarkerColor(4); aPt_coscos_os->SetLineColor(4); 
  aPt_sinsin_os->SetMarkerStyle(24); aPt_sinsin_os->SetMarkerColor(4); aPt_sinsin_os->SetLineColor(4);
  aPt_coscos_ss->GetXaxis()->SetTitle("#LTp_{T}#GT"); aPt_coscos_ss->GetYaxis()->SetTitle("#LTcos()cos()#GT or #LTsin()sin()#GT");
  //
  dPt_coscos_ss->SetMarkerStyle(21); dPt_coscos_ss->SetMarkerColor(2); dPt_coscos_ss->SetLineColor(2); 
  dPt_sinsin_ss->SetMarkerStyle(25); dPt_sinsin_ss->SetMarkerColor(2); dPt_sinsin_ss->SetLineColor(2); 
  dPt_coscos_os->SetMarkerStyle(20); dPt_coscos_os->SetMarkerColor(4); dPt_coscos_os->SetLineColor(4); 
  dPt_sinsin_os->SetMarkerStyle(24); dPt_sinsin_os->SetMarkerColor(4); dPt_sinsin_os->SetLineColor(4); 
  dPt_coscos_ss->GetXaxis()->SetTitle("|#Deltap_{T}|"); dPt_coscos_ss->GetYaxis()->SetTitle("#LTcos()cos()#GT or #LTsin()sin()#GT");
  
  
  complex<double> binput(1.0,0.0);
  
  const int MAXPIONS = 2000;
  const int NevtBuff = 1;
  float P[NevtBuff][MAXPIONS][4]={{{0}}};
  float X[NevtBuff][MAXPIONS][4]={{{0}}};
  int PID[NevtBuff][MAXPIONS]={{0}};
  bool Primary[NevtBuff][MAXPIONS]={{0}};
  int Nparticles[NevtBuff]={0};
  float P1_PRF[3]={0};
  float P2_PRF[3]={0};
  float X1_PRF[3]={0};
  float X2_PRF[3]={0};
  
  
  
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
  
  
  int seconds = time(0);
  
  for(int FileN=0; FileN<Nfiles; FileN++){
    cout<<"File #"<<FileN+1<<endl;
    
    //TString *fname=new TString("b2/event");
    //*fname += FileN+1;
    TString *fname=new TString("../ThermGenerator/events/lhyquid2dbi-RHICAuAu200c6070Ti279ti025Tf145/event");
    //TString *fname=new TString("../ThermGenerator/events/lhyquid2dbi-RHICAuAu200c5060Ti354ti025Tf145/event");
    //TString *fname=new TString("../ThermGenerator/events/lhyquid2dbi-RHICAuAu200c4050Ti398ti025Tf145/event");
    //TString *fname=new TString("../ThermGenerator/events/lhyquid2dbi-RHICAuAu200c5060Ti354ti025Tf145/event");
    //TString *fname=new TString("../ThermGenerator/events/lhyquid2dbi-RHICAuAu200c3040Ti429ti025Tf145/event");
    //TString *fname=new TString("../ThermGenerator/events/lhyquid2dbi-RHICAuAu200c2030Ti455ti025Tf145/event");
    //TString *fname=new TString("../ThermGenerator/events/lhyquid2dbi-RHICAuAu200c1020Ti476ti025Tf145/event");
    //TString *fname=new TString("../ThermGenerator/events/lhyquid2dbi-RHICAuAu200c0510Ti491ti025Tf145/event");
    //TString *fname=new TString("../ThermGenerator/events/lhyquid2dbi-RHICAuAu200c0005Ti500ti025Tf145/event");

    if(FileN<100 && FileN>=10) fname->Append("0");
    else if(FileN<10) fname->Append("00");
    else {}
    *fname += FileN;
    fname->Append(".root");
    TFile *infile=new TFile(fname->Data(),"READ");
    TTree *event_tree=(TTree*)infile->Get("events");
    TTree *particles_tree=(TTree*)infile->Get("particles");
    

    int Nevents = event_tree->GetEntries();
    int NpartList=0;
    int NpartList_past = 0;
    /////////////////////////////////////////
    // Create Pion list first
    for(int Evt=0; Evt<Nevents; Evt++){
      
      //cout<<"Event # "<<Evt+1<<endl;
      event_tree->GetEntry(Evt);
      
      TBranch *br=(TBranch*)event_tree->GetBranch("event");
      TLeaf *lf=(TLeaf*)br->GetLeaf("entries");
      NpartList_past += NpartList;
      NpartList = lf->GetValue();
      
      /////////////////////////////////////
      // past event buffer
      for(int pastEvt=NevtBuff-1; pastEvt>0; pastEvt--){
	
	for(int particle=0; particle<Nparticles[pastEvt-1]; particle++){
	  P[pastEvt][particle][0] = P[pastEvt-1][particle][0];
	  P[pastEvt][particle][1] = P[pastEvt-1][particle][1];
	  P[pastEvt][particle][2] = P[pastEvt-1][particle][2];
	  P[pastEvt][particle][3] = P[pastEvt-1][particle][3];
	  //
	  X[pastEvt][particle][0] = X[pastEvt-1][particle][0];
	  X[pastEvt][particle][1] = X[pastEvt-1][particle][1];
	  X[pastEvt][particle][2] = X[pastEvt-1][particle][2];
	  X[pastEvt][particle][3] = X[pastEvt-1][particle][3];
	  //
	  PID[pastEvt][particle] = PID[pastEvt-1][particle];
	  Primary[pastEvt][particle] = Primary[pastEvt-1][particle];
	}
	
	Nparticles[pastEvt]=Nparticles[pastEvt-1];
      }
      Nparticles[0]=0;
      /////////////////////////////////////
      
     
      // Create Pion list first
      //Int_t randomIndex[NpartList];
      //for (Int_t ri = 0; ri < NpartList; ri++) randomIndex[ri]=NpartList_past+ri;
      //Shuffle(randomIndex, 0, NpartList-1);
      //
      int Npions=0;
      
      for(int index1=NpartList_past; index1<NpartList_past + NpartList; index1++){
	
	if(Npions >= MAXPIONS) {cout<<"Too Many Pions!!!"<<endl; break;}
	particles_tree->GetEntry(index1);
	//particles_tree->GetEntry(randomIndex[index1]);
	//
	
	TLeaf *lf_pid=particles_tree->GetLeaf("pid");
	int pid = int(lf_pid->GetValue());
	if(abs(pid) != 211) continue;// pions only
	TLeaf *lf_decayed=(TLeaf*)particles_tree->GetLeaf("decayed");
	if(lf_decayed->GetValue() != 0) continue;
	
	
	TLeaf *lf_px=(TLeaf*)particles_tree->GetLeaf("px");
	TLeaf *lf_py=(TLeaf*)particles_tree->GetLeaf("py");
	TLeaf *lf_pz=(TLeaf*)particles_tree->GetLeaf("pz");
	TLeaf *lf_x=(TLeaf*)particles_tree->GetLeaf("x");
	TLeaf *lf_y=(TLeaf*)particles_tree->GetLeaf("y");
	TLeaf *lf_z=(TLeaf*)particles_tree->GetLeaf("z");
	TLeaf *lf_t=(TLeaf*)particles_tree->GetLeaf("t");
	
	float px = lf_px->GetValue(); float py = lf_py->GetValue(); float pz = lf_pz->GetValue();
	float x = lf_x->GetValue(); float y = lf_y->GetValue(); float z = lf_z->GetValue();
	float t = lf_t->GetValue();
	
	float pt = sqrt(pow(px,2)+pow(py,2));
	float eta = asinh(pz/pt);
	float E = sqrt(pow(px,2)+pow(py,2)+pow(pz,2)+pow(0.13957,2));
	if(pt<0.15 || pt>2.0) continue;// 0.16 for ALICE. 0.15 for STAR 
	if(fabs(eta)>1.0) continue;// 0.8 for ALICE. 1.0 for STAR
	//if(sqrt(pow(px,2)+pow(py,2)+pow(pz,2)) > 1.0) continue;
	//
	P[0][Npions][0] = E; P[0][Npions][1] = px; P[0][Npions][2] = py; P[0][Npions][3] = pz;
	X[0][Npions][0] = t; X[0][Npions][1] = x; X[0][Npions][2] = y; X[0][Npions][3] = z;
	PID[0][Npions] = pid;
	if(sqrt(pow(x,2)+pow(y,2)+pow(z,2)) < 1e10) Primary[0][Npions]=kTRUE; 
	else {Primary[0][Npions]=kFALSE;}
	
	
	
	if(Primary[0][Npions]){
	  float r = sqrt(pow(X[0][Npions][1],2) + pow(X[0][Npions][2],2) + pow(X[0][Npions][3],2));
	  r_dist->Fill(r);
	  x_dist->Fill(fabs(X[0][Npions][1]));
	  Pt_dist->Fill(pt);
	}
	//
	
	
	Npions++;
      }
      
      Nparticles[0]=Npions;
      
      // Shuffle
      for (Int_t i = 0; i < Npions; i++) {
	// generate random index call
	int j = (Int_t) (gRandom->Rndm() * Npions);
	// temporarily store values from random index
	float t = X[0][j][0], x = X[0][j][1], y = X[0][j][2], z = X[0][j][3];
	float E = P[0][j][0], px = P[0][j][1], py = P[0][j][2], pz = P[0][j][3];
	int pid = PID[0][j];
	bool prim = Primary[0][j];
	//
	// Swap value locations
	X[0][j][0] = X[0][i][0]; X[0][j][1] = X[0][i][1]; X[0][j][2] = X[0][i][2]; X[0][j][3] = X[0][i][3];
	P[0][j][0] = P[0][i][0]; P[0][j][1] = P[0][i][1]; P[0][j][2] = P[0][i][2]; P[0][j][3] = P[0][i][3];
	PID[0][j] = PID[0][i];
	Primary[0][j] = Primary[0][i];
	//
	X[0][i][0] = t; X[0][i][1] = x; X[0][i][2] = y; X[0][i][3] = z;
	P[0][i][0] = E; P[0][i][1] = px; P[0][i][2] = py; P[0][i][3] = pz;
	PID[0][i] = pid;
	Primary[0][i] = prim;
      }
      

      //////////////////////////////////////////////////
      // Start Correlation Analysis

      for(int EN=0; EN<NevtBuff; EN++){

	for(int index1=0; index1<Nparticles[EN]; index1++){// current or past event pion loop
	  
	  int start=index1+1;
	  if(EN>0) start=0;
	  for(int index2=start; index2<Npions; index2++){// current event pion loop (EN=0)
	    
	    float kt = sqrt(pow(P[EN][index1][1]+P[0][index2][1],2)+pow(P[EN][index1][2]+P[0][index2][2],2))/2.;
	    float qinv12 = sqrt(pow(P[EN][index1][1]-P[0][index2][1],2)+pow(P[EN][index1][2]-P[0][index2][2],2)+pow(P[EN][index1][3]-P[0][index2][3],2)-pow(P[EN][index1][0]-P[0][index2][0],2));
	    
	    if(qinv12 < 0.001) continue;
	   
	  
	    if(EN>0){
	      if(PID[EN][index1] == PID[0][index2]){
		DenEM_ss->Fill(kt, qinv12);
	      }else {
		DenEM_os->Fill(kt, qinv12);
	      }
	      continue;
	    }
	    
	    if(PID[EN][index1] == PID[0][index2]) DenLambda_ss->Fill(kt, qinv12);
	    else DenLambda_os->Fill(kt, qinv12);
	    

	  // boost into PRF
	  BoostPRF(P[EN][index1], P[0][index2], P[EN][index1], P1_PRF);
	  BoostPRF(P[EN][index1], P[0][index2], P[0][index2], P2_PRF);
	  BoostPRF(P[EN][index1], P[0][index2], X[EN][index1], X1_PRF);
	  BoostPRF(P[EN][index1], P[0][index2], X[0][index2], X2_PRF);
	 
	  //cout<<"out: "<<P1_PRF[0]<<"  "<<P2_PRF[0]<<endl;
	  //cout<<"side: "<<P1_PRF[1]<<"  "<<P2_PRF[1]<<endl;
	  //cout<<"long: "<<P1_PRF[2]<<"  "<<P2_PRF[2]<<" ++++++"<<endl;
	  //cout<<qinv12<<"  "<<sqrt(pow(P1_PRF[0]-P2_PRF[0],2)+pow(P1_PRF[1]-P2_PRF[1],2)+pow(P1_PRF[2]-P2_PRF[2],2))<<endl;
	  //
	  float rstar12_PRF = sqrt(pow(X1_PRF[0]-X2_PRF[0],2) + pow(X1_PRF[1]-X2_PRF[1],2) + pow(X1_PRF[2]-X2_PRF[2],2));
	  if(rstar12_PRF<100){
	    if(PID[EN][index1] == PID[0][index2]) NumLambda_ss->Fill(kt, qinv12, rstar12_PRF);
	    else NumLambda_os->Fill(kt, qinv12, rstar12_PRF);
	  }

	  if(rstar12_PRF < 0.1) continue;// removes pairs from the same decay
	  
	  //
	  float rho = rstar12_PRF/FmToGeV * qinv12/2.;
	  float phase = ((P1_PRF[0]-P2_PRF[0])*(X1_PRF[0]-X2_PRF[0]) + (P1_PRF[1]-P2_PRF[1])*(X1_PRF[1]-X2_PRF[1]) + (P1_PRF[2]-P2_PRF[2])*(X1_PRF[2]-X2_PRF[2]))/2.;
	  phase /= FmToGeV;
	  float Cosphase = phase/(rstar12_PRF/FmToGeV)/(qinv12/2.);
	  //float Cosphase = 2*RandomNumber->Rndm()-1.0;
	  
	  float eta = 1/(qinv12/2.*BohrR);
	  if(PID[EN][index1] != PID[0][index2]) eta = -eta;
	  complex<double> ainput(0.0, -eta);
	  complex<double> zinput12(0.0, rho + rho*Cosphase);
	  complex<double> zinput21(0.0, rho - rho*Cosphase);
	  complex<double> CHG12(1.0, 0);
	  complex<double> CHG21(1.0, 0);
	  float GamovFactor=1.0;
	  float BEE = 0.;
	  float WeightWF = 1.0;
	  // FSI
	  complex<double> Psi_alpha(0., 0.);
	  if(IncludeFSI && qinv12 < 0.2 && rstar12_PRF<100.){
	    if(Primary[EN][index1]==kTRUE && Primary[0][index2]==kTRUE){
	      if(PID[EN][index1]!=PID[0][index2]){
		double kb = sqrt(pow(masspiC,2)+pow(qinv12/2.,2) - pow(masspi0,2));
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
		Psi_alpha = fc_aa*G/(rstar12_PRF/FmToGeV);
	      }
	  
	      CHG12= conf_hyp(ainput,binput,zinput12);
	      CHG21= conf_hyp(ainput,binput,zinput21);
	      GamovFactor = Gamov(eta);
	      BEE = cos(phase);
	      
	      complex<double> planewave12(double(TMath::Cos(-rho*Cosphase)), double(TMath::Sin(-rho*Cosphase)));
	      complex<double> planewave21(planewave12.real(), -planewave12.imag());
	      //
	      complex<double> Full12(0.0,0.0);
	      Full12 += planewave12;
	      Full12 *= CHG12;
	      if(PID[EN][index1]!=PID[0][index2]) Full12 += Psi_alpha;
	      //
	      complex<double> Full21(0.0,0.0);
	      Full21 += planewave21;
	      Full21 *= CHG21;
	      //
	      complex<double> FullWF = Full12;
	      if(PID[EN][index1] == PID[0][index2]) {
		FullWF += Full21;
		FullWF /= sqrt(2.);
	      }	  

	      WeightWF = GamovFactor*pow(abs(FullWF),2);
	  
	      // Fill undiluted histos
	      if(PID[EN][index1] == PID[0][index2]){
		Num_PrimCosFSI_ss->Fill(kt, qinv12, WeightWF);
	      }else {
		Num_PrimCosFSI_os->Fill(kt, qinv12, WeightWF);
	      }
	      
	    }// primary
	  }// IncludeFSI and rstar and qinv cut
	  
	  // LPV studies
	  double phi1=atan2(P[EN][index1][2],P[EN][index1][1]);
	  double phi2=atan2(P[0][index2][2],P[0][index2][1]);
	  double pt1=sqrt(pow(P[EN][index1][1],2)+pow(P[EN][index1][2],2));
	  double pt2=sqrt(pow(P[0][index2][1],2)+pow(P[0][index2][2],2));
	  double eta1=asinh(P[EN][index1][3]/pt1);
	  double eta2=asinh(P[0][index2][3]/pt2);
	  double ccW = cos(phi1)*cos(phi2)*WeightWF;
	  double ssW = sin(phi1)*sin(phi2)*WeightWF;

	  if(PID[EN][index1] == PID[0][index2]){
	    rstarPRF_dist_ss->Fill(kt, rstar12_PRF);
	    tstar_dist_ss->Fill(kt, fabs(X[EN][index1][0]-X[0][index2][0]));
	    Num_Cos_ss->Fill(kt, qinv12, BEE);
	    Num_CosFSI_ss->Fill(kt, qinv12, WeightWF);
	    Den_ss->Fill(kt, qinv12);

	    //
	    ThreePoint_ss->Fill(1, ccW-ssW);
	    aEta_coscos_ss->Fill(fabs(eta1+eta2)/2., ccW);
	    aEta_sinsin_ss->Fill(fabs(eta1+eta2)/2., ssW);
	    dEta_coscos_ss->Fill(fabs(eta1-eta2), ccW);
	    dEta_sinsin_ss->Fill(fabs(eta1-eta2), ssW);
	    aPt_coscos_ss->Fill(fabs(pt1+pt2)/2., ccW);
	    aPt_sinsin_ss->Fill(fabs(pt1+pt2)/2., ssW);
	    dPt_coscos_ss->Fill(fabs(pt1-pt2), ccW);
	    dPt_sinsin_ss->Fill(fabs(pt1-pt2), ssW);
	  }else {
	    rstarPRF_dist_os->Fill(kt, rstar12_PRF);
	    tstar_dist_os->Fill(kt, fabs(X[EN][index1][0]-X[0][index2][0]));
	    Num_Cos_os->Fill(kt, qinv12, BEE);
	    Num_CosFSI_os->Fill(kt, qinv12, WeightWF);
	    Den_os->Fill(kt, qinv12);
	    //
	    ThreePoint_os->Fill(1, ccW-ssW);
	    aEta_coscos_os->Fill(fabs(eta1+eta2)/2., ccW);
	    aEta_sinsin_os->Fill(fabs(eta1+eta2)/2., ssW);
	    dEta_coscos_os->Fill(fabs(eta1-eta2), ccW);
	    dEta_sinsin_os->Fill(fabs(eta1-eta2), ssW);
	    aPt_coscos_os->Fill(fabs(pt1+pt2)/2., ccW);
	    aPt_sinsin_os->Fill(fabs(pt1+pt2)/2., ssW);
	    dPt_coscos_os->Fill(fabs(pt1-pt2), ccW);
	    dPt_sinsin_os->Fill(fabs(pt1-pt2), ssW);
	  }
	  
	
	}// index2
	
      }// index1
      
      }// EN buffer      
    }// Nevents
    
    infile->Close();
  }// Nfiles
  
  outfile->cd();
  //
  r_dist->Write();
  x_dist->Write();
  Pt_dist->Write();
  rstarPRF_dist_ss->Write();
  tstar_dist_ss->Write();
  rstarPRF_dist_os->Write();
  tstar_dist_os->Write();
  NumLambda_ss->Write();
  NumLambda_os->Write();
  DenLambda_ss->Write();
  DenLambda_os->Write();
  Den_ss->Write();
  Den_os->Write();
  DenEM_ss->Write();
  DenEM_os->Write();
  Num_Cos_ss->Write();
  Num_Cos_os->Write();
  Num_CosFSI_ss->Write();
  Num_CosFSI_os->Write();
  //
  Num_PrimCosFSI_ss->Write();
  Num_PrimCosFSI_os->Write();
  //
  //
  ThreePoint_ss->Write();
  ThreePoint_os->Write();
  aEta_coscos_ss->Write();
  aEta_sinsin_ss->Write();
  aEta_coscos_os->Write();
  aEta_sinsin_os->Write();
  dEta_coscos_ss->Write();
  dEta_sinsin_ss->Write();
  dEta_coscos_os->Write();
  dEta_sinsin_os->Write();
  aPt_coscos_ss->Write();
  aPt_sinsin_ss->Write();
  aPt_coscos_os->Write();
  aPt_sinsin_os->Write();
  dPt_coscos_ss->Write();
  dPt_sinsin_ss->Write();
  dPt_coscos_os->Write();
  dPt_sinsin_os->Write();
  //
  outfile->Close();

  seconds = time(0) - seconds;
  cout<<"Minutes = "<<seconds/60.<<endl;


}
  float Gamov(float eta){
  return (2*PI*eta/(exp(2*PI*eta)-1));
}
float fact(float n){
  return (n < 1.00001) ? 1 : fact(n - 1) * n;
}
void Shuffle(Int_t *iarr, Int_t i1, Int_t i2)
{
  Int_t j, k;
  Int_t a = i2 - i1;
  for (Int_t i = i1; i < i2+1; i++) {
    j = (Int_t) (gRandom->Rndm() * a);
    k = iarr[j];
    iarr[j] = iarr[i];
    iarr[i] = k;
  }
}
void BoostPRF(float P1[4], float P2[4], float V[4], float *T){
  // P1 is particle 1 momentum in reaction CMS (0=E,1=px,2=py,3=pz)
  // P2 is particle 2 momentum in reaction CMS (0=E,1=px,2=py,3=pz)
  // V is the vector to be transformed (V is in CMS)
  // T is the PRF vector

  float MT=sqrt(pow(P1[0]+P2[0],2)-pow(P1[3]+P2[3],2));
  float PT=sqrt(pow(P1[1]+P2[1],2)+pow(P1[2]+P2[2],2));
  float MINV=sqrt(pow(P1[0]+P2[0],2) - pow(P1[1]+P2[1],2) - pow(P1[2]+P2[2],2) - pow(P1[3]+P2[3],2));
  // LCMS
  T[0] = ((P1[1]+P2[1])*V[1] + (P1[2]+P2[2])*V[2])/PT;
  T[1] = ((P1[1]+P2[1])*V[2] - (P1[2]+P2[2])*V[1])/PT;
  T[2] = ((P1[0]+P2[0])*V[3] - (P1[3]+P2[3])*V[0])/MT;
  // PRF
  //float temp = T[0];
  T[0] = T[0]*MINV/MT - PT/(MT*MINV)*((P1[0]+P2[0])*V[0] - (P1[1]+P2[1])*V[1] - (P1[2]+P2[2])*V[2] - (P1[3]+P2[3])*V[3]);
  //T[0] = T[0]*MT/MINV;
  //cout<<T[0]<<"  "<<temp*MINV/MT - PT/(MT*MINV)*((P1[0]+P2[0])*V[0] - (P1[1]+P2[1])*V[1] - (P1[2]+P2[2])*V[2] - (P1[3]+P2[3])*V[3])<<endl;
  //T[0] = T[0]*masspiC/(2.*sqrt(pow(masspiC,2)+pow((P1[1]+P2[1]),2)+pow((P1[2]+P2[2]),2)));
}
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
double KinverseFunction(double kstarSq, int J){
  double ESq = kstarSq + pow(masspiC,2);
  double E = sqrt(ESq);
  double xSq = kstarSq/pow(masspiC,2);
  double Kinverse = E*(4*ESq-fzero[J][4])/(4*pow(masspiC,2)-fzero[J][4]);
  Kinverse /= fzero[J][0] + fzero[J][1]*xSq + fzero[J][2]*pow(xSq,2) + fzero[J][3]*pow(xSq,3);
  return Kinverse;
}

