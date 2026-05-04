#include "MomentHistogram.h"
#include "TString.h"


ClassImp(MomentHistogram)

//-----------------------------------------

MomentHistogram::MomentHistogram(){
    InitHistograms();
}

//-----------------------------------------

MomentHistogram::~MomentHistogram(){}

//-----------------------------------------
void MomentHistogram::InitHistograms() {
    
    // -----**-----**-----**-----**-----**--p_{T} Graphs--**-----**-----**-----**-----**
    
    // Cumulants (Absolute)
    gpTC1 = new TGraphErrors();
    gpTC2 = new TGraphErrors();
    gpTC3 = new TGraphErrors();
    gpTC4 = new TGraphErrors();
    gpTC1->SetName("[p_{T}] C1");
    gpTC1->SetMarkerStyle(20);
    gpTC1->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gpTC1->GetYaxis()->SetTitle("C_{1}");
    gpTC1->GetXaxis()->SetRangeUser(0.8, 1.25);
    gpTC2->SetName("[p_{T}] C2");
    gpTC2->SetMarkerStyle(21);
    gpTC2->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gpTC2->GetYaxis()->SetTitle("C_{2}");
    gpTC2->GetXaxis()->SetRangeUser(0.8, 1.25);
    gpTC3->SetName("[p_{T}] C3");
    gpTC3->SetMarkerStyle(22);
    gpTC3->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gpTC3->GetYaxis()->SetTitle("C_{3}");
    gpTC3->GetXaxis()->SetRangeUser(0.8, 1.25);
    gpTC4->SetName("[p_{T}] C4");
    gpTC4->SetMarkerStyle(33);
    gpTC4->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gpTC4->GetYaxis()->SetTitle("C_{4}");
    gpTC4->GetXaxis()->SetRangeUser(0.8, 1.25);
    
    // Cumulants (Normalized to 0-5%)
    gpTC1N = new TGraphErrors();
    gpTC2N = new TGraphErrors();
    gpTC3N = new TGraphErrors();
    gpTC4N = new TGraphErrors();
    gpTC1N->SetName("[p_{T}] C1 normalized");
    gpTC1N->SetMarkerStyle(20);
    gpTC1N->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gpTC1N->GetYaxis()->SetTitle("C_{1}/C_{1}^{0-5%}");
    gpTC1N->GetXaxis()->SetRangeUser(0.8, 1.25);
    gpTC2N->SetName("[p_{T}] C2 normalized");
    gpTC2N->SetMarkerStyle(21);
    gpTC2N->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gpTC2N->GetYaxis()->SetTitle("C_{2}/C_{2}^{0-5%}");
    gpTC2N->GetXaxis()->SetRangeUser(0.8, 1.25);
    gpTC3N->SetName("[p_{T}] C3 normalized");
    gpTC3N->SetMarkerStyle(22);
    gpTC3N->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gpTC3N->GetYaxis()->SetTitle("C_{3}/C_{3}^{0-5%}");
    gpTC3N->GetXaxis()->SetRangeUser(0.8, 1.25);
    gpTC4N->SetName("[p_{T}] C4 normalized");
    gpTC4N->SetMarkerStyle(33);
    gpTC4N->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gpTC4N->GetYaxis()->SetTitle("C_{4}/C_{4}^{0-5%}");
    gpTC4N->GetXaxis()->SetRangeUser(0.8, 1.25);
    
    // Moments
    gpTM = new TGraphErrors();
    gpTVar = new TGraphErrors();
    gpTS = new TGraphErrors();
    gpTkap = new TGraphErrors();
    gpTM->SetName("[p_{T}] Mean");
    gpTM->SetMarkerStyle(24);
    gpTM->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gpTM->GetYaxis()->SetTitle("#LT[p_{T}]#GT");
    gpTM->GetXaxis()->SetRangeUser(0.8, 1.25);
    gpTVar->SetName("[p_{T}] Variance");
    gpTVar->SetMarkerStyle(25);
    gpTVar->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gpTVar->GetYaxis()->SetTitle("Var([p_{T}])");
    gpTVar->GetXaxis()->SetRangeUser(0.8, 1.25);
    gpTS->SetName("[p_{T}] Skewness");
    gpTS->SetMarkerStyle(26);
    gpTS->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gpTS->GetYaxis()->SetTitle("S([p_{T}])");
    gpTS->GetXaxis()->SetRangeUser(0.8, 1.25);
    gpTkap->SetName("[p_{T}] Kurtosis");
    gpTkap->SetMarkerStyle(27);
    gpTkap->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gpTkap->GetYaxis()->SetTitle("#kappa([p_{T}])");
    gpTkap->GetXaxis()->SetRangeUser(0.8, 1.25);
    
    // Moments (Normalized to 0-5%)
    gpTMN = new TGraphErrors();
    gpTVarN = new TGraphErrors();
    gpTSN = new TGraphErrors();
    gpTkapN = new TGraphErrors();
    gpTMN->SetName("[p_{T}] Mean normalized");
    gpTMN->SetMarkerStyle(24);
    gpTMN->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gpTMN->GetYaxis()->SetTitle("#LT[p_{T}]#GT/#LT[p_{T}]#GT^{0-5%}");
    gpTMN->GetXaxis()->SetRangeUser(0.8, 1.25);
    gpTVarN->SetName("[p_{T}] Variance normalized");
    gpTVarN->SetMarkerStyle(25);
    gpTVarN->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gpTVarN->GetYaxis()->SetTitle("Var([p_{T}])/Var([p_{T}])^{0-5%}");
    gpTVarN->GetXaxis()->SetRangeUser(0.8, 1.25);
    gpTSN->SetName("[p_{T}] Skewness normalized");
    gpTSN->SetMarkerStyle(26);
    gpTSN->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gpTSN->GetYaxis()->SetTitle("S([p_{T}])/S([p_{T}])^{0-5%}");
    gpTSN->GetXaxis()->SetRangeUser(0.8, 1.25);
    gpTkapN->SetName("[p_{T}] Kurtosis normalized");
    gpTkapN->SetMarkerStyle(27);
    gpTkapN->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gpTkapN->GetYaxis()->SetTitle("#kappa([p_{T}])/#kappa([p_{T}])^{0-5%}");
    gpTkapN->GetXaxis()->SetRangeUser(0.8, 1.25);
    
    // -----**-----**-----**-----**-----**--N_{ch} Graphs--**-----**-----**-----**-----**
    
    // Cumulants (Absolute)
    gNchC1 = new TGraphErrors();
    gNchC2 = new TGraphErrors();
    gNchC3 = new TGraphErrors();
    gNchC4 = new TGraphErrors();
    gNchC1->SetName("N_{ch} C1");
    gNchC1->SetMarkerStyle(20);
    gNchC1->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gNchC1->GetYaxis()->SetTitle("C_{1}");
    gNchC2->SetName("N_{ch} C2");
    gNchC2->SetMarkerStyle(21);
    gNchC2->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gNchC2->GetYaxis()->SetTitle("C_{2}");
    gNchC3->SetName("N_{ch} C3");
    gNchC3->SetMarkerStyle(22);
    gNchC3->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gNchC3->GetYaxis()->SetTitle("C_{3}");
    gNchC4->SetName("N_{ch} C4");
    gNchC4->SetMarkerStyle(33);
    gNchC4->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gNchC4->GetYaxis()->SetTitle("C_{4}");

    // Cumulants (Normalized to 0-5%)
    gNchC1N = new TGraphErrors();
    gNchC2N = new TGraphErrors();
    gNchC3N = new TGraphErrors();
    gNchC4N = new TGraphErrors();
    gNchC1N->SetName("N_{ch} C1 normalized");
    gNchC1N->SetMarkerStyle(20);
    gNchC1N->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gNchC1N->GetYaxis()->SetTitle("C_{1}/C_{1}^{0-5%}");
    gNchC2N->SetName("N_{ch} C2 normalized");
    gNchC2N->SetMarkerStyle(21);
    gNchC2N->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gNchC2N->GetYaxis()->SetTitle("C_{2}/C_{2}^{0-5%}");
    gNchC3N->SetName("N_{ch} C3 normalized");
    gNchC3N->SetMarkerStyle(22);
    gNchC3N->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gNchC3N->GetYaxis()->SetTitle("C_{3}/C_{3}^{0-5%}");
    gNchC4N->SetName("N_{ch} C4 normalized");
    gNchC4N->SetMarkerStyle(33);
    gNchC4N->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gNchC4N->GetYaxis()->SetTitle("C_{4}/C_{4}^{0-5%}");

    // Moments
    gNchM = new TGraphErrors();
    gNchVar = new TGraphErrors();
    gNchS = new TGraphErrors();
    gNchkap = new TGraphErrors();
    gNchM->SetName("N_{ch} Mean");
    gNchM->SetMarkerStyle(24);
    gNchM->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gNchM->GetYaxis()->SetTitle("#LTN_{ch}#GT");
    gNchVar->SetName("N_{ch} Variance");
    gNchVar->SetMarkerStyle(25);
    gNchVar->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gNchVar->GetYaxis()->SetTitle("Var(N_{ch})");
    gNchS->SetName("N_{ch} Skewness");
    gNchS->SetMarkerStyle(26);
    gNchS->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gNchS->GetYaxis()->SetTitle("S(N_{ch})");
    gNchkap->SetName("N_{ch} Kurtosis");
    gNchkap->SetMarkerStyle(27);
    gNchkap->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gNchkap->GetYaxis()->SetTitle("#kappa(N_{ch})");

    // Moments (Normalized to 0-5%)
    gNchMN = new TGraphErrors();
    gNchVarN = new TGraphErrors();
    gNchSN = new TGraphErrors();
    gNchkapN = new TGraphErrors();
    gNchMN->SetName("N_{ch} Mean normalized");
    gNchMN->SetMarkerStyle(24);
    gNchMN->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gNchMN->GetYaxis()->SetTitle("#LTN_{ch}#GT/#LTN_{ch}#GT^{0-5%}");
    gNchVarN->SetName("N_{ch} Variance normalized");
    gNchVarN->SetMarkerStyle(25);
    gNchVarN->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gNchVarN->GetYaxis()->SetTitle("Var(N_{ch})/Var(N_{ch})^{0-5%}");
    gNchSN->SetName("N_{ch} Skewness normalized");
    gNchSN->SetMarkerStyle(26);
    gNchSN->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gNchSN->GetYaxis()->SetTitle("S(N_{ch})/S(N_{ch})^{0-5%}");
    gNchkapN->SetName("N_{ch} Kurtosis normalized");
    gNchkapN->SetMarkerStyle(27);
    gNchkapN->GetXaxis()->SetTitle("#LTdN_{ch}/d#eta#GT/#LTdN_{ch}/d#eta#GT_{0-5%}");
    gNchkapN->GetYaxis()->SetTitle("#kappa(N_{ch})/#kappa(N_{ch})^{0-5%}");
}
