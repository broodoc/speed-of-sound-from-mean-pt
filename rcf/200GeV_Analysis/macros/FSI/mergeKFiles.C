#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <Riostream.h>
#include <complex>

#include "TObject.h"
#include "TFile.h"
#include "TString.h"
#include "TF1.h"
#include "TH1.h"

using namespace std;

void mergeKFiles( vector<TString> filenames ) {

  //---------------------------------------
  // R VALUES NEED TO BE MANUALLY SET HERE!!!!!!
  int Rbins = filenames.size();
  double Rmin = 3 - 0.25;
  double Rmax = 10 + 0.25;

  TFile *fprobe = new TFile( filenames[0], "READ" );
  TH1D *probe = (TH1D*)fprobe->Get( "K2_ss" );
  int Qbins = probe->GetNbinsX();
  double Qmin = probe->GetXaxis()->GetBinLowEdge( 1 );
  double Qmax = probe->GetXaxis()->GetBinUpEdge( Qbins );
  fprobe->Close();

  TFile *outfile = new TFile("KFile.root","RECREATE");
  
  TH2D *K2_ss = new TH2D("K2_ss","K2ss Gaussian source;q_{inv} (GeV/c);R_{inv} (fm)", Qbins,Qmin,Qmax, Rbins,Rmin,Rmax);
  TH2D *K2_os = new TH2D("K2_os","K2os Gaussian source;q_{inv} (GeV/c);R_{inv} (fm)", Qbins,Qmin,Qmax, Rbins,Rmin,Rmax);

  for( int i=0; i<filenames.size(); i++ ) {

    TFile *f = new TFile(filenames[i], "READ");
    if( ! f ) { 
      cout<<"File not found"<<endl;
      break;
    }
    TH1D *ss = (TH1D*)f->Get("K2_ss");
    TH1D *os = (TH1D*)f->Get("K2_os");

    for( int q = 1; q <= ss->GetNbinsX(); q++ ) {
      K2_ss->SetBinContent( q, i + 1, ss->GetBinContent( q ) );
      K2_os->SetBinContent( q, i + 1, os->GetBinContent( q ) );
    }
  }

  outfile->cd();
  K2_ss->Write();
  K2_os->Write();

  outfile->Close();

}
