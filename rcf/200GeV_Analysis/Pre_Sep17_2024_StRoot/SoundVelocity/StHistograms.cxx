#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"

#include "StHistograms.h"
#include "StParameters.h"
#include "TString.h" /// <-- Necessary?

ClassImp(StHistograms)

  //------------------------------------------------------------------------
StHistograms::StHistograms( StParameters *params )
{

  InitHistograms( params );

}

//------------------------------------------------------------------------
StHistograms::~StHistograms()
{
}

//------------------------------------------------------------------------
void StHistograms::InitHistograms( StParameters *params ) 
{
 /// pt_plot
 /// N_{ch}
 /// Missing definitions in StHistograms.h and appropriate references in StSoundVelocity and other files
    TString eta_label; TString temp_title; TString temp_plot_title; TString pt_label;
    for (int i = 0; i < 13; i+=4) {
        if (i == 0) {eta_label = "Eta0p5_";} /// abs(#eta) < 0.5
        else if (i == 4) {eta_label = "Eta0p8_";} /// abs(#eta) < 0.8
        else if (i == 8) {eta_label = "Eta1p0_";} /// abs(#eta) < 1.0
        else if (i == 12) {eta_label = "Full_TPC_Eta_1p5_";} /// abs(#eta) < 1.5
        else {continue;}
        for (int j = 0; j < 4; j++) {
            if (j == 0) {pt_label = "_pTlow0p0_";} /// abs(pT) > 0
            else if (j == 1) {pt_label = "_pTlow0p15_";} /// abs(pT) > 0.15
            else if (j == 2) {pt_label = "_pTlow0p4_";} /// abs(pT) > 0.40
            else if (j == 3) {pt_label = "_pTlow0p8_";} /// abs(pT) > 0.80
            else {continue;}
            ///--------------------------ETA DIST---------------------------------------
            temp_title = eta_label + pt_label + "hEta";
        temp_plot_title = "#eta Distribution (Acceptance" + eta_label + pt_label + "); #eta; Tracks";
            hEta[i + j] = new TH1D(temp_title.Data(), temp_plot_title.Data(), 320, -1.6,1.6 );
            ///--------------------------ETA DIST---------------------------------------
            ///--------------------------ETA/PHI DIST---------------------------------------
            temp_title = eta_label + pt_label + "hEtaPhi";
        temp_plot_title = "#eta  vs  #phi Distribution (Acceptance_" + eta_label + pt_label + ");#eta;#phi";
            hEtaPhi[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),60,-3,3, 100,-TMath::Pi(),TMath::Pi());
            ///--------------------------ETA/PHI DIST---------------------------------------
            ///--------------------------Min Bias p_{T} DIST---------------------------------------
//            temp_title = eta_label + pt_label + "hTransvMomentum";
//        temp_plot_title = "p_{T}  Distribution (Acceptance" + eta_label + pt_label + ");p_{T} (GeV/c); Tracks";
//            hTransvMomentum[i + j] = new TH1D(temp_title.Data(), temp_plot_title.Data(), 1100, 0, 11);
            ///--------------------------Min Bias p_{T} DIST---------------------------------------
            ///--------------------------2D pT_Nch/RefMult DISTS---------------------------------------
            temp_title = eta_label + pt_label + "hPt_Nch";
        temp_plot_title = "2D p_{T}  vs N_{ch} (Acceptance" + eta_label + pt_label + "); N_{ch}; p_{T} (GeV/c)";
            hPt_Nch[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),1501,-0.5,1500.5,1100,0,11);
	temp_title = eta_label + pt_label + "hPt_Eta";
        temp_plot_title = "p_{T} vs Eta (Acceptance_" + eta_label+ pt_label + ");p_{T} (GeV/c);#eta";
             hPtEta[i + j] = new TH2D(temp_title.Data(), "p_{T} vs Eta;p_{T} (GeV);#eta", 500,0,5, 320,-1.6,1.6);
//            temp_title = eta_label + pt_label + "hPt_RefMult";
//        temp_plot_title = "2D p_{T}  vs RefMult (Acceptance" + eta_label + pt_label + "); RefMult; p_{T} (GeV/c)";
//            hPt_RefMult[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),1501,-0.5,1500.5,1100,0,11);
//            temp_title = eta_label + pt_label + "hPt_dNdEta";
//        temp_plot_title = "2D p_{T}  vs dN_{ch}/d{#eta} (approx.) (Acceptance" + eta_label + pt_label + "); dN_{ch}/d#eta; p_{T} (GeV/c)";
//            hPt_dNdEta[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),1501,-0.5,1500.5,1100,0,11);
            ///--------------------------2D pT_Nch/RefMult DISTS---------------------------------------
            ///--------------------------1D Multiplicity DISTS---------------------------------------
            temp_title = eta_label + pt_label + "hNch";
        temp_plot_title = "N_{ch} Distribution (Acceptance" + eta_label + pt_label + ");N_{ch}; Events";
            hMult[i + j] = new TH1D(temp_title.Data(),temp_plot_title.Data(),1501,-0.5,1500.5); /// <- filled w/ mytracks.size() w/ all necessary cuts
//            temp_title = "hChargedTracks" + eta_label;
//            temp_plot_title = "N_{ch} Distribution (p_{T} #geq 0.15 GeV/c & Acceptance" + eta_label + "); N_{ch}; Events";
//            hChargedTracks[i + j] = new TH1D(temp_title.Data(),temp_plot_title.Data(),1501,-0.5,1500.5); /// <- filled w/ mytracks.size() but only: | eta | & p_{T} >= 0.15 GeV/c cuts
            ///--------------------------1D Multiplicity DISTS---------------------------------------
            ///--------------------------Charge_p_{T} DIST-----------------------------------------
//            temp_title = eta_label + pt_label + "hChargePt";
//        temp_plot_title = "Charge vs p_{T} (Acceptance" + eta_label + pt_label + ");Charge;p_{T} (GeV/c)";
//            hChargePt[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),2,-1.5,1.5, 1100,0,11);
            ///--------------------------Charge_p_{T} DIST-----------------------------------------
            ///--------------------------DCA_p_{T} DIST-----------------------------------------
            temp_title = eta_label + pt_label + "hPt_DCA";
        temp_plot_title = "p_{T} vs DCA (Acceptance" + eta_label+ pt_label + ");p_{T} (GeV/c);DCA (cm)";
            hPtDCA[i + j] = new TH2D(temp_title.Data(), "p_{T} vs DCA;p_{T} (GeV);DCA (cm)", 500,0,5, 1100,0,11);
            ///--------------------------DCA_p_{T} DIST-----------------------------------------
        }
    }
///--------------------------Vertex DIST-----------------------------------------------
		hVertex = new TH3D("hVertex","Primary Vertex Distribution;Vx (cm);Vy (cm);Vz (cm)",40,-2.,2., 80,-4.,4., 240,-30.,210.);
///--------------------------Vertex DIST-----------------------------------------------
///--------------------------hMultRaw and hMultCent (no eta cut)
		hMultRaw = new TH1D("hMultRaw","Raw N_{ch} Distribution;Nch",2001,-0.5,2000.5);
		hMultCent = new TH2D("hMultCent",";Refmult;centrality1to16", 1501,-0.5,1500.5, 16,-0.5,15.5);
		hRefMult = new TH1D("RefMult","RefMult Distribution;RefMult",1501,-0.5,1500.5); /// <- filled w/ StCorr->RefMultCorr Class
}
