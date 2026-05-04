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
    TString eta_label; TString temp_title; TString temp_plot_title; TString pt_label;
     int j{0};
    for (int i = 0; i < 1; i++) {
	if (i == 0) {eta_label = "dEta1p0_";} /// dEta is 1p0
        else if ( i == 1 ) {eta_label = "Centrality_";} // (Centrality determination QA) 
            if (j == 0) {pt_label = "_pTlow0p20_";} /// abs(pT) > 0.20
	    ///--------------------------Vz DIST---------------------------------------
	    temp_title = eta_label + pt_label + "hVz";
	    temp_plot_title = "Vz Distribution (Acceptance_" + eta_label + pt_label + "); Vz (cm); Events";
	    hVz[i + j] = new TH1D(temp_title.Data(), temp_plot_title.Data(), 6, -6,6 );
	    ///--------------------------Vz DIST--------------------------------------
            ///--------------------------ETA DIST---------------------------------------
        temp_title = eta_label + pt_label + "hEta";
        temp_plot_title = "#eta Distribution (Acceptance_" + eta_label + pt_label + "); #eta; Tracks";
            hEta[i + j] = new TH1D(temp_title.Data(), temp_plot_title.Data(), 31, -1.55,1.55 );
            ///--------------------------ETA DIST---------------------------------------
            ///--------------------------ETA/PHI DIST---------------------------------------
        temp_title = eta_label + pt_label + "hEtaPhi";
        temp_plot_title = "#eta  vs  #phi Distribution (Acceptance_" + eta_label + pt_label + ");#eta;#phi";
            hEtaPhi[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),31,-1.55,1.55, 100,-TMath::Pi(),TMath::Pi());
            ///--------------------------ETA/PHI DIST---------------------------------------
            ///--------------------------Eta vs Nch DIST---------------------------------------
        temp_title = eta_label + pt_label + "hEta_Nch";
        temp_plot_title = "#eta  vs  N_{ch} Distribution (Acceptance_" + eta_label + pt_label + ");#eta;N_{ch}";
            hEta_Nch[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),31,-1.55,1.55, 801,-0.5,800.5);
            ///--------------------------Eta vs Nch DIST---------------------------------------
	    ///--------------------------phi vs Nch DIST---------------------------------------
        temp_title = eta_label + pt_label + "hPhi_Nch";
        temp_plot_title = "#phi  vs  N_{ch} Distribution (Acceptance_" + eta_label + pt_label + ");#phi;N_{ch}";
            hPhi_Nch[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),100,-TMath::Pi(),TMath::Pi(), 801,-0.5,800.5);
            ///--------------------------phi vs Nch DIST---------------------------------------
            ///--------------------------Eta vs Nhits DIST---------------------------------------
       temp_title = eta_label + pt_label + "hEta_Nhits";
       temp_plot_title = "#eta  vs  Nhits Distribution (Acceptance_" + eta_label + pt_label + ");#eta_{lab};Nhits";
            hEta_Nhits[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),31,-1.5,1.5, 81,-0.05,80.55);
            ///--------------------------Eta vs Nhits DIST---------------------------------------     
	    ///--------------------------Nch vs NhitsFit DIST---------------------------------------
	    temp_title = eta_label + pt_label + "hNch_NhitsFit";
       temp_plot_title = "NhitsFit  vs  N_{ch} Distribution (Acceptance_" + eta_label + pt_label + ");N_{ch}; NhitsFit";
            hNch_NhitsFit[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),801,-0.5,800.5, 81,-0.5,80.5);
	    ///--------------------------Nch vs NhitsFit DIST---------------------------------------
	    ///--------------------------Nch vs NhitsdEdx DIST---------------------------------------
            temp_title = eta_label + pt_label + "hNch_NhitsdEdx";
       temp_plot_title = "NhitsdEdx  vs  N_{ch} Distribution (Acceptance_" + eta_label + pt_label + ");N_{ch}; NhitsdEdx";
            hNch_NhitsdEdx[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),801,-0.5,800.5, 81,-0.5,80.5);
            ///--------------------------Nch vs Nhits DIST---------------------------------------
	    ///--------------------------Nch vs DCA DIST---------------------------------------
            temp_title = eta_label + pt_label + "hNch_DCA";
       temp_plot_title = "DCA  vs  N_{ch} Distribution (Acceptance_" + eta_label + pt_label + ");N_{ch}; DCA (cm)";
            hNch_DCA[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),801,-0.5,800.5, 21,-0.005,2.505);
            ///--------------------------Nch vs DCA DIST---------------------------------------
	temp_title = eta_label + pt_label + "hPt_dNdEta";
        temp_plot_title = "2D p_{T}  vs dN_{ch}/d{#eta} (approx.) (Acceptance_" + eta_label + pt_label + "); dN_{ch}/d#eta; p_{T} (GeV/c)";
            //hPt_dNdEta[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),601,-0.5,600.5,1100,0,11);
            ///--------------------------2D pT_Nch/RefMult DISTS---------------------------------------
            ///--------------------------1D Multiplicity DISTS-----------------------------------------
        temp_title = eta_label + pt_label + "hNch";
        temp_plot_title = "N_{ch} Distribution (Acceptance_" + eta_label + pt_label + ");N_{ch}; Events";
            hMult[i + j] = new TH1D(temp_title.Data(),temp_plot_title.Data(),801,-0.5,800.5); /// <- filled w/ mytracks.size() w/ 
            ///--------------------------1D Multiplicity DISTS---------------------------------------
            ///--------------------------1D Multiplicity DISTS-----------------------------------------
        temp_title = eta_label + pt_label + "hRefMultCorr";
        temp_plot_title = "RefMultCorr Distribution (Acceptance_" + eta_label + pt_label + ");N_{ch}; Events";
            hRefMultCorr[i + j] = new TH1D(temp_title.Data(),temp_plot_title.Data(),801,-0.5,800.5); /// <- filled w/ CorrectedRefMult + weight 
            ///--------------------------1D Multiplicity DISTS---------------------------------------
	    ///--------------------------1D Multiplicity DISTS---------------------------------------
	    temp_title = eta_label + pt_label + "hRefMult";
        temp_plot_title = "RefMult Distribution (Acceptance_" + eta_label + pt_label + ");N_{ch}; Events";
            hRefMult[i + j] = new TH1D(temp_title.Data(),temp_plot_title.Data(),801,-0.5,800.5); /// <- filled w/ refMult (no weights or corrections)
	    ///--------------------------1D Multiplicity DISTS---------------------------------------
	    ///--------------------------Charge_p_{T} DIST-----------------------------------------
        temp_title = eta_label + pt_label + "hChargePt";
        temp_plot_title = "Charge vs p_{T} (Acceptance_" + eta_label + pt_label + ");Charge;p_{T} (GeV/c)";
        //    hChargePt[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),2,-1.5,1.5, 110,0,11);
            ///--------------------------Charge_p_{T} DIST-----------------------------------------
            ///--------------------------DCA_p_{T} DIST-----------------------------------------
        temp_title = eta_label + pt_label + "hDCA_pT";
        temp_plot_title = "DCA vs p_{T} (Acceptance_" + eta_label+ pt_label + ");DCA (cm); p_{T} (GeV/c)";
             hPtDCA[i + j] = new TH2D(temp_title.Data(), "DCA vs p_{T};DCA (cm);p_{T} (GeV/c)",325,0,3.25,500,0,5);
            ///--------------------------DCA_p_{T} DIST-----------------------------------------
	    ///--------------------------Eta_p_{T} DIST-----------------------------------------
	 temp_title = eta_label + pt_label + "hEta_pT";
        temp_plot_title = "#eta vs p_{T} (Acceptance_" + eta_label+ pt_label + ");#eta; p_{T} (GeV/c)";
             hPtEta[i + j] = new TH2D(temp_title.Data(), "Eta vs p_{T}; #eta; p_{T} (GeV/c)",31,-1.55,1.55,500,0,5);
	    ///--------------------------Eta_p_{T} DIST-----------------------------------------
            ///--------------------------2D FXTMult_nBTOFmatch DIST-----------------------------------------
        temp_title = eta_label + pt_label + "hRefMult_nBTOFmatch";
        temp_plot_title = "2D FXTMult vs nBTOFmatch (Acceptance_" + eta_label + pt_label + "); nBTofMatch; RefMult";
            hFXTMult_nBTOFmatch[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),601,-0.5,600.5,1001,-0.5,1000.5);
            ///--------------------------2D FXTMult_nBTOFmatch DIST-----------------------------------------
            ///--------------------------2D FXTMult_RunID DIST-----------------------------------------
	temp_title = eta_label + pt_label + "hRefMult_RunID";
        temp_plot_title = "2D RefMult vs Run ID (Acceptance_" + eta_label + pt_label + "); Run ID; refMult";
            hFXTMult_RunID[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),6000,17.038e6,17.179012e6,801,-0.5,800.5);
	    ///--------------------------2D FXTMult_RunID DIST-----------------------------------------
	    ///--------------------------2D nBTOFmatch_RunID DIST-----------------------------------------
	//temp_title = eta_label + pt_label + "hnBTofMatch_RunID";
        //temp_plot_title = "2D nBTofMatch vs Run ID (Acceptance_" + eta_label + pt_label + "); Run ID; nBTOFmatch";
        //    hnBTOFmatch_RunID[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),6000,19.151e6,19.157e6,401,-0.5,400.5);
	    ///--------------------------2D nBTOFmatch_RunID DIST-----------------------------------------
	    ///--------------------------2D BBCx_RunID DIST-----------------------------------------
	//temp_title = eta_label + pt_label + "hBBCx_RunID";
        //temp_plot_title = "2D BBCx vs Run ID (Acceptance " + eta_label + pt_label + "); Run ID; BBCx";
            //hBBCx_RunID[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),6000,19.151e6,19.157e6,451,-0.5,4500.5);
	    ///--------------------------2D BBCx_RunID DIST-----------------------------------------
	    //cout << "\n Initializing BBCx_nBTOFmatch Histogram[ " << i+j << " ]" << endl;
	    ///--------------------------2D BBCx_nBTofMatch DIST-----------------------------------------
	temp_title = eta_label + pt_label + "hBBCx_nBTOFmatch";
        temp_plot_title = "2D BBCx vs nBTOFmatch (Acceptance_" + eta_label + pt_label + "); nBTOFmatch; BBCx";
            hBBCx_nBTOFmatch[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),601,-0.5,600.5,1101,2e4,10e4);	    
	    ///--------------------------2D BBCx_nBTofMatch DIST-----------------------------------------
	    //cout << "\n Initializing BBCx_FXTMult Histogram[ " << i+j << " ]\n" << endl;
	    ///--------------------------2D BBCx_FXTMult DIST-----------------------------------------
	temp_title = eta_label + pt_label + "hBBCx_FXTMult";
        temp_plot_title = "2D BBCx vs FXTMult (Acceptance_" + eta_label + pt_label + "); RefMult; BBCx";
            hBBCx_FXTMult[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),701,-0.5,700.5,1101,2e4,10e4);
	    ///--------------------------2D BBCx_FXTMult DIST-----------------------------------------
	    ///--------------------------2D ZDCx_FXTMult DIST-----------------------------------------
	temp_title = eta_label + pt_label + "hZDCx_FXTMult";
	temp_plot_title = "2D ZDCx vs FXTMult (Acceptance_" + eta_label + pt_label + "); RefMult; ZDCx";
		hZDCx_FXTMult[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),701,-0.5,700.5,100,2e4,10e4);
	    ///--------------------------2D ZDCx_FXTMult DIST-----------------------------------------
	    ///--------------------------2D ZDCx_nBTofmatch DIST-----------------------------------------
	temp_title = eta_label + pt_label + "hZDCx_nBTOFmatch";
	temp_plot_title = "2D ZDCx vs nBTOFmatch (Acceptance_" + eta_label + pt_label + "); nBTOFmatch; ZDCx";
		hZDCx_nBTOFmatch[i + j] = new TH2D(temp_title.Data(),temp_plot_title.Data(),601,-0.5,600.5,100,2e4,10e4);
	    ///--------------------------2D ZDCx_nBTofmatch DIST-----------------------------------------
	//}
    }
	///--------------------------Nch vs StRefMultCorr------------------------------------------
	temp_title = "hNch_RefMultCorr";
        temp_plot_title = "StRefMultCorr - N_{ch} Correlation; StRefMultCorr; N_{ch}";
        hNch_RefMultCorr = new TH2D(temp_title.Data(),temp_plot_title.Data(), 800, 0, 800, 800, 0, 800);
	///--------------------------Nch vs StRefMultCorr------------------------------------------
	///--------------------------ZDC East vs West----------------------------------------------
	temp_title = "hZdcEastWest";
	temp_plot_title = "ZDC East vs West; ZDC East; ZDC West";
	//hZdcEastWest = new TH2D(temp_title.Data(),temp_plot_title.Data(), 4500, 0, 4500, 4500, 0, 4500);
	///--------------------------ZDC East vs West----------------------------------------------
	///--------------------------RefMult vs ZDC----------------------------------------------
	temp_title = "hRefMult_Zdc";
        temp_plot_title = "RefMult vs ZDC; gRefMult; ZDC [East + West]";
        //hRefMult_Zdc = new TH2D(temp_title.Data(),temp_plot_title.Data(), 701, -0.5, 700.5, 4000, 0, 8000);
	///--------------------------RefMult vs ZDC  ---------------------------------------------
	///--------------------------ZDC East + West ----------------------------------------------
	temp_title = "hZdcEastpWest";
        temp_plot_title = "ZDC; ZDC [East + West]; Counts (a.u.)";
        //hZdcEastPWest = new TH1D(temp_title.Data(),temp_plot_title.Data(), 4000, 0, 8000);
	///--------------------------ZDC East + West ----------------------------------------------
	///--------------------------2D pT (analysis criteria) vs Centrality Bin DIST--------------
 	temp_title = "pT_Centrality";
	temp_plot_title = "Centrality vs pT;; p_{T} (GeV/c) [Analysis criteria]";
	//hpTCent = new TH2D(temp_title.Data(),temp_plot_title.Data(),17, 1, 18, 1001, -0.05, 10.05);
	//std::vector<std::string> binLabels = {
    //"75-80%", "70-75%", "65-70%", "60-65%", "55-60%", "50-55%", 
    //"45-50%", "40-45%", "35-40%", "30-35%", "25-30%", "20-25%", 
    //"15-20%", "10-15%", "5-10%", "0-5%"
	//};
	//for (int i = 1; i < static_cast<int>( binLabels.size() ); i++) {
    	//hpTCent->GetXaxis()->SetBinLabel(i, binLabels[i].c_str()); 
	//}
	///--------------------------------------------------------------------------------------
	///--------------------------2D Nch (centrality criteria) vs Centrality Bin DIST---------
	temp_title = "Nch_Centrality";
        temp_plot_title = "Centrality vs N_{ch};; N_{ch} [Centrality criteria]";
        //hNchCent = new TH2D(temp_title.Data(),temp_plot_title.Data(),17, 1, 18, 801,-0.5,800.5);
        //for (int i = 1; i < static_cast<int>( binLabels.size() ); i++) {
        //hNchCent->GetXaxis()->SetBinLabel(i, binLabels[i].c_str());
        //}	 
	///--------------------------------------------------------------------------------------
	
	///--------------------------2D pT (analysis criteria) vs 9 Centrality Bin DIST--------------
        //temp_title = "pT_Centrality9";
        //temp_plot_title = "Centrality vs pT;; p_{T} (GeV/c) [Analysis criteria]";
        //hpTCent9 = new TH2D(temp_title.Data(),temp_plot_title.Data(),17, 1, 18, 1001, -0.05, 10.05);
        //std::vector<std::string> binLabels9 = {
    //"70-80%","60-70%", "50-60%", "40-50%", "30-40%", "20-30%", "10-20%",
    //"5-10%", "0-5%"};
        //for (int i = 1; i < static_cast<int>( binLabels9.size() ); i++) {
        //hpTCent9->GetXaxis()->SetBinLabel(i, binLabels9[i].c_str());
        //}
        ///--------------------------------------------------------------------------------------
        ///--------------------------2D Nch (centrality criteria) vs 9 Centrality Bin DIST---------
        //temp_title = "Nch_Centrality9";
        //temp_plot_title = "Centrality vs N_{ch};; N_{ch} [Centrality criteria]";
        //hNchCent9 = new TH2D(temp_title.Data(),temp_plot_title.Data(),17, 1, 18, 801,-0.5,800.5);
        //for (int i = 1; i < static_cast<int>( binLabels9.size() ); i++) {
        //hNchCent9->GetXaxis()->SetBinLabel(i, binLabels9[i].c_str());
        //}
        ///--------------------------------------------------------------------------------------	
	
	///--------------------------3D pT_Nch_eta DISTS---------------------------------------
        temp_title = eta_label + pt_label + "hPt_Nch_eta";
        temp_plot_title = "3D p_{T} [Ana]  vs N_{ch} [Cent] vs eta [Ana]; N_{ch}; p_{T} (GeV/c); #eta";
            hPt_Nch_eta = new TH3D(temp_title.Data(),temp_plot_title.Data(),801,-0.5,800.5,50,0,5,300, -1.5, 1.5);
        ///--------------------------3D pT_Nch_eta DISTS---------------------------------------
	///--------------------------2D pT_Nch DIST---------------------------------------
        temp_title = eta_label + pt_label + "hPt_Nch";
        temp_plot_title = "2D p_{T} [Analysis] vs N_{ch} [Charged Tracks]; N_{ch}; p_{T} (GeV/c)";
            hPt_Nch = new TH2D(temp_title.Data(),temp_plot_title.Data(),801,-0.5,800.5,1000,0,10);
        ///--------------------------2D pT_Nch/FXTMult DISTS---------------------------------------
	///--------------------------2D pT_Nch DIST---------------------------------------
	temp_title = eta_label + pt_label + "hPt_RefMultCorr";
        temp_plot_title = "2D p_{T} [Analysis] vs N_{ch} [RefMultCorr]; N_{ch}; p_{T} (GeV/c)";
            hPt_RefMultCorr = new TH2D(temp_title.Data(),temp_plot_title.Data(),801,-0.5,800.5,1000,0,10);
	///--------------------------2D pT_Nch DIST---------------------------------------	
	///--------------------------2D <pT> vs Nch DIST-----------------------------------------
	temp_title = "<pT>_Nch";
	temp_plot_title = "<pT> vs Nch; N_{ch}; <p_{T}> (GeV/c)";
	    hMeanpTNch = new TH2D(temp_title.Data(),temp_plot_title.Data(),801, -0.5, 800.5, 801, -0.05, 8.05); 
	///--------------------------------------------------------------------------------------
	///--------------------------Event Cut Statistics-----------------------------------------
	temp_title = "EventCuts";
	temp_plot_title = "Event Cut Statistics;;Events";
	hEventCuts = new TH1D(temp_title.Data(), temp_plot_title.Data(), 10, 0, 10);
	hEventCuts->GetXaxis()->SetBinLabel(1, "Before Cuts");
	hEventCuts->GetXaxis()->SetBinLabel(2, "Triggers");
	hEventCuts->GetXaxis()->SetBinLabel(3, "|V_{z}| < 6cm");
	hEventCuts->GetXaxis()->SetBinLabel(4, "|V_{z} - Vz_{VPD}|<3cm");
	hEventCuts->GetXaxis()->SetBinLabel(5, "!(V_{x}=0&V_{y}=0&V_{z}=0)");
	hEventCuts->GetXaxis()->SetBinLabel(6, "|V_{r}| < 2cm");
	hEventCuts->GetXaxis()->SetBinLabel(7, "refMultCorr run rej. ");
	hEventCuts->GetXaxis()->SetBinLabel(8, "manual run rej. ");
	hEventCuts->GetXaxis()->SetBinLabel(9, "refMultCorr pileup rej. ");
	hEventCuts->GetXaxis()->SetBinLabel(10, "manual pileup rej.");
	///--------------------------------------------------------------------------------------	
	///--------------------------Track Cut Statistics-----------------------------------------
        temp_title = "TrackCuts";
        temp_plot_title = "Track Cut Statistics;;Tracks";
        hTrackCuts = new TH1D(temp_title.Data(), temp_plot_title.Data(), 10, 0, 10);
	hTrackCuts->GetXaxis()->SetBinLabel(1, "Before Cuts");
	hTrackCuts->GetXaxis()->SetBinLabel(2, "!pTrack");
	hTrackCuts->GetXaxis()->SetBinLabel(3, "!isPrimary");
	//hTrackCuts->GetXaxis()->SetBinLabel(4, "HFT cut");
	hTrackCuts->GetXaxis()->SetBinLabel(4, "DCA cut");
	//hTrackCuts->GetXaxis()->SetBinLabel(6, "nHits cut");
	hTrackCuts->GetXaxis()->SetBinLabel(5, "nHitsFit cut");
	hTrackCuts->GetXaxis()->SetBinLabel(6, "nHitsDdDx cut");
	hTrackCuts->GetXaxis()->SetBinLabel(7, "nHitsFit/nHitsMax cut");
	hTrackCuts->GetXaxis()->SetBinLabel(8, "phi cut (TPC)");
	hTrackCuts->GetXaxis()->SetBinLabel(9, "#eta cut");
	hTrackCuts->GetXaxis()->SetBinLabel(10, "p_{T} cuts");	
	///--------------------------2D Multiplicity DISTS-----------------------------------------
        temp_title = "hNch2D";
        temp_plot_title = "N_{ch} Correlation; Refmult [|#eta| < 0.5]; RefMult2 [0.5 < |#eta| < 1.0]";
            hMult2D = new TH2D(temp_title.Data(),temp_plot_title.Data(),801, -0.5, 800.5, 801, -0.5, 800.5);// <- filled w/ mytracks.size() w/
        ///--------------------------2D Multiplicity DISTS---------------------------------------
//		temp_title = "hPt_FXTMult";
    //    	temp_plot_title = "2D p_{T}  vs FXTMult; FXTMult; p_{T} (GeV/c)";
  //          hPt_FXTMult = new TH2D(temp_title.Data(),temp_plot_title.Data(),601,-0.5,600.5,1100,0,11);
///--------------------------Vertex DIST-----------------------------------------------
		hVertex = new TH3D("hVertex","Primary Vertex Distribution;Vx (cm);Vy (cm);Vz (cm)",50,-2.5,2.5, 80,-4.,4., 80, -40,40.);
///--------------------------Vertex DIST-----------------------------------------------
///--------------------------hMultRaw and hMultCent (no eta cut)
	//	hMultRaw = new TH1D("hMultRaw","Raw N_{ch} Distribution;Nch",2001,-0.5,2000.5);
	//	hMultCent = new TH2D("hMultCent",";Refmult;centrality1to16", 1501,-0.5,1500.5, 16,-0.5,15.5);
	//	hRefMult = new TH1D("FXTMult","FXTMult Distribution;FXTMult",1501,-0.5,1500.5); /// <- filled w/ StCorr->RefMultCorr Class
}
