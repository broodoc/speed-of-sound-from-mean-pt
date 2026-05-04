#ifndef STUTILITIES_CXX
#define STUTILITIES_CXX

#include "TH1F.h"
#include "TH2F.h"

#include "StUtilities.h"
#include "StParameters.h"
#include "StHistograms.h"

ClassImp(StUtilities)

  //------------------------------------------------------------------------
StUtilities::StUtilities( TString setting )
{
  mSetting = setting;
  mFXT = (mSetting.Contains("FXT", TString::ECaseCompare::kIgnoreCase)) ? true : false;
//  LoadFSIfile( Kfilename );
//  LoadMRCfile( MRCfilename );
//  LoadWeightfile( Weightfilename );
}

//------------------------------------------------------------------------
void StUtilities::FillPicoRawTracks( StPicoDst *picoDst, StParameters *params, vector<RawTrackStruct> *tracks, StHistograms *hm ) {

  TVector3 priVtx = picoDst->event()->primaryVertex();

  SetiTPCCase( picoDst->event()->runId() );

  for( UInt_t i = 0; i < picoDst->numberOfTracks(); i++ ) {
    RawTrackStruct myRawTrack;
    hm->hTrackCuts->Fill(0.5);
    StPicoTrack *pTrack = picoDst->track( i );
    if( ! pTrack ) continue;
    hm->hTrackCuts->Fill(1.5);
    Bool_t isPrimary = pTrack->isPrimary();
    if( !isPrimary ) continue;
    hm->hTrackCuts->Fill(2.5);
    TVector3 momentum= pTrack->pMom();
    
    //myRawTrack.mTrack.mFlag = pTrack->flag(); // ??
    
    myRawTrack.mTrack.mCharge = pTrack->charge();
    myRawTrack.mTrack.mFirstWordTopologyMap = pTrack->topologyMap(0);
    myRawTrack.mTrack.mSecondWordTopologyMap = pTrack->topologyMap(1);
    myRawTrack.mTrack.mThirdWordTopologyMap = pTrack->iTpcTopologyMap();
    myRawTrack.mHFThitsMap = pTrack->hftHitsMap();
    myRawTrack.mTrack.mP = momentum;
    myRawTrack.mTrack.mPt = momentum.Perp();
    myRawTrack.mTrack.mPhi = momentum.Phi();
    myRawTrack.mTrack.mEta = momentum.PseudoRapidity();
    myRawTrack.mTrack.mMass = constants::mass_pion; // pion mass assumption
    myRawTrack.mTrack.mE = sqrt( momentum.Mag2() + pow(myRawTrack.mTrack.mMass,2) );

    myRawTrack.mRawTrackID = pTrack->id();
    myRawTrack.mDCA = pTrack->gDCA(priVtx.x(), priVtx.y(), priVtx.z());
    myRawTrack.mNhits = pTrack->nHits();
    myRawTrack.mNhitsFit = pTrack->nHitsFit();
    myRawTrack.mNhitsMax = pTrack->nHitsMax();
    myRawTrack.mNhitsDedx = pTrack->nHitsDedx();
    myRawTrack.mNsigmaPion = pTrack->nSigmaPion();
    myRawTrack.mNsigmaElectron = pTrack->nSigmaElectron();
    myRawTrack.mNsigmaKaon = pTrack->nSigmaKaon();
    myRawTrack.mNsigmaProton = pTrack->nSigmaProton();

    myRawTrack.mBTOFindex = pTrack->bTofPidTraitsIndex();
    if( myRawTrack.mBTOFindex >= 0 ) {
      StPicoBTofPidTraits *TOFtraits = (StPicoBTofPidTraits*)picoDst->btofPidTraits( myRawTrack.mBTOFindex );
      if( ! TOFtraits ) { myRawTrack.mBTOFmatchFlag = 0; }
      else { 
        myRawTrack.mBTOFmatchFlag = TOFtraits->btofMatchFlag();
        if( myRawTrack.mBTOFmatchFlag != 0 ) {
          myRawTrack.mNsigmaPionTOF = TOFtraits->nSigmaPion();
          myRawTrack.mBTOFylocal = TOFtraits->btofYLocal();
          myRawTrack.mBTOFzlocal = TOFtraits->btofZLocal();
          double beta = TOFtraits->btofBeta(); 
          double betaPion = momentum.Mag() / (sqrt( momentum.Mag2() + pow(constants::mass_pion,2) ));
          double m2 = momentum.Mag2() * ( 1./(beta*beta) - 1);
          myRawTrack.mM2 = m2;
          myRawTrack.mBeta = beta;
          myRawTrack.mBetaPion = betaPion;
        }
      }
    }
    else {
      myRawTrack.mBTOFmatchFlag = -1;
      myRawTrack.mNsigmaPionTOF = -1;
      myRawTrack.mBTOFylocal = -1;
      myRawTrack.mBTOFzlocal = -1;
      myRawTrack.mM2 = -1;
      myRawTrack.mBeta = -1;
      myRawTrack.mBetaPion = -1;
    }

    tracks->push_back( myRawTrack );
  }

}

//------------------------------------------------------------------------
void StUtilities::FillMuRawTracks( StMuDst *muDst, StParameters *params, vector<RawTrackStruct> *tracks ) {

  StThreeVectorF Stvtx = muDst->primaryVertex()->position(); 
  TVector3 priVtx( Stvtx.x(), Stvtx.y(), Stvtx.z() );

  SetiTPCCase( muDst->event()->runId() );

  // numberOfGlobalTracks() or numberOfPrimaryTracks()??
  for( UInt_t i = 0; i < muDst->numberOfPrimaryTracks(); i++ ) {
    RawTrackStruct myRawTrack;

    StMuTrack *muTrack = muDst->primaryTracks( i );
    if( ! muTrack ) continue;
    short isPrimary = muTrack->type(); //returns: 0=global; 1=primary
    if( isPrimary == 0 ) continue;

    StThreeVectorF Stmom = muTrack->momentum();
    TVector3 momentum( Stmom.x(), Stmom.y(), Stmom.z() );

    myRawTrack.mTrack.mCharge = muTrack->charge();
    myRawTrack.mTrack.mFirstWordTopologyMap = muTrack->topologyMap().data(0);
    myRawTrack.mTrack.mSecondWordTopologyMap = muTrack->topologyMap().data(1);
    myRawTrack.mTrack.mThirdWordTopologyMap = muTrack->topologyMap().data(2);
    myRawTrack.mTrack.mP = momentum;
    myRawTrack.mTrack.mPt = momentum.Perp();
    myRawTrack.mTrack.mPhi = momentum.Phi();
    myRawTrack.mTrack.mEta = momentum.PseudoRapidity();
    myRawTrack.mTrack.mMass = constants::mass_pion; // pion mass assumption
    myRawTrack.mTrack.mE = sqrt( momentum.Mag2() + pow(myRawTrack.mTrack.mMass,2) );

    myRawTrack.mRawTrackID = muTrack->id();
    myRawTrack.mDCA = muTrack->dcaGlobal().mag();
    myRawTrack.mNhits = muTrack->nHits();
    myRawTrack.mNhitsDedx = muTrack->nHitsDedx();
    myRawTrack.mNsigmaPion = muTrack->nSigmaPion();
    myRawTrack.mNsigmaElectron = muTrack->nSigmaElectron();
    myRawTrack.mNsigmaKaon = muTrack->nSigmaKaon();
    myRawTrack.mNsigmaProton = muTrack->nSigmaProton();

    myRawTrack.mBTOFindex = muTrack->index2BTofHit();
    StMuBTofPidTraits TOFtraits = muTrack->btofPidTraits();
    myRawTrack.mBTOFmatchFlag = TOFtraits.matchFlag();
    if( myRawTrack.mBTOFindex >= 0 ) {
      if( myRawTrack.mBTOFmatchFlag != 0 ) {
        myRawTrack.mNsigmaPionTOF = TOFtraits.sigmaPion();
        myRawTrack.mBTOFylocal = TOFtraits.yLocal();
        myRawTrack.mBTOFzlocal = TOFtraits.zLocal();
        double beta = TOFtraits.beta(); 
        double betaPion = momentum.Mag() / (sqrt( momentum.Mag2() + pow(constants::mass_pion,2) ));
        double m2 = momentum.Mag2() * ( 1./(beta*beta) - 1);
        myRawTrack.mM2 = m2;
        myRawTrack.mBeta = beta;
        myRawTrack.mBetaPion = betaPion;
      }
    }
    else {
      myRawTrack.mNsigmaPionTOF = -1;
      myRawTrack.mBTOFylocal = -1;
      myRawTrack.mBTOFzlocal = -1;
      myRawTrack.mM2 = -1;
      myRawTrack.mBeta = -1;
      myRawTrack.mBetaPion = -1;
    }

//    double QB = myRawTrack.mTrack.mCharge * muDst->event()->magneticField(); // B needs to be in kiloGauss
//    CalculateCrossings( params, momentum, QB, myRawTrack.mTrack.mCrossings );
//
//    // MC matching
//    myRawTrack.mTrack.mID_MC = muTrack->idTruth();
//    GetMCtruthQuantities( muDst, &myRawTrack );

    tracks->push_back( myRawTrack );
  }

}


//------------------------------------------------------------------------
bool StUtilities::GoodPicoEvent( StParameters *params, StPicoEvent *event, StHistograms *hm ) {
  //return true; //Debugging: Remove this or comment for official submission!
  TVector3 priVtx = event->primaryVertex();
  float VpdVz = event->vzVpd();
  /// Trigger
  bool goodTrigger = false;
  for( uint i = 0; i < params->paramMap["Trigger"].size(); i++ ) {
    if( event->isTrigger( round(params->paramMap["Trigger"][i]) ) ) { goodTrigger = true; }
  }
  if( ! goodTrigger ) { return false; }
  hm->hEventCuts->Fill(1.5);
  /// Vz
  if( priVtx.z() < params->paramMap["VzMin"][0] ) { return false; }
  if( priVtx.z() > params->paramMap["VzMax"][0] ) { return false; }
  if( priVtx.z() == 0 ) { return false; } 
  hm->hEventCuts->Fill(2.5);
  /// Vz - Vpd
  /// |Vz_{TPC} - Vz_{VPD}| < cut
 if (!mFXT){
 if (priVtx.z() - VpdVz > params->paramMap["VzVpdMax"][0] || priVtx.z() - VpdVz < params->paramMap["VzVpdMin"][0]) {return false;}
 hm->hEventCuts->Fill(3.5);
 if( std::abs(priVtx.x()) <  10e-6 && std::abs(priVtx.y()) < 10e-6 && std::abs(priVtx.z()) < 10e-6 ) { return false; }
 hm->hEventCuts->Fill(4.5); 
 }
//cout<<"past Vz"<<endl;
  /// Vr
  double VyShift = 0;
  if( mFXT ) {
    VyShift = -2; // -2 cm Y offset of fixed target
  }
  double Vr = sqrt( pow( priVtx.x(), 2) + pow( priVtx.y() - VyShift, 2) );
  if( Vr > params->paramMap["VrCut"][0] ) { return false; }
  hm->hEventCuts->Fill(5.5);
//Ignore this:  //if( BadRun( params, event->runId() ) ) { return false; }
	//cout<<"past Vr"<<endl;
  
  return true;
}

//------------------------------------------------------------------------
bool StUtilities::GoodMuEvent( StParameters *params, StMuEvent *event ) {

  StThreeVectorF priVtx = event->primaryVertexPosition();
//  float VpdVz = event->vpdVz();
  // Trigger
  if( ! event->triggerIdCollection().nominal().isTrigger( round(params->paramMap["Trigger"][0]) ) ) { return false; }

  /// Vz
  //if( priVtx.z() < params->paramMap["VzMin"][0] ) { return false; }
  //if( priVtx.z() > params->paramMap["VzMax"][0] ) { return false; }
  //if (priVtx.z() - VpdVz > params->paramMap["VzVpdMax"][0] || priVtx.z() - VpdVz < params->paramMap["VzVpdMin"][0]) {return false;}   
  /// Vr
  //double VyShift = 0;
  //if( mFXT ) {
  //  VyShift = -2; // -2 cm Y offset of fixed target
  //}
  //double Vr = sqrt( pow( priVtx.x(), 2) + pow( priVtx.y() - VyShift, 2) );
  //if( Vr > params->paramMap["VrCut"][0] ) { return false; }
  
  //Bad runs
  if( BadRun( params, event->runId() ) ) { return false; }

  return true;
}

//------------------------------------------------------------------------
bool StUtilities::GoodTrack( StParameters *params, RawTrackStruct rawTrack, StHistograms *hm ) {
// Apply track cuts for pT vs Nch histogram: nHitsMin, DCA, nHitsMin is 15, pT > 0.15
//return true; //Debugging: Remove this or comment for official submission!
bool firstCall = true;
  //mHFThitsMap HFT hits (PXL1, PXL2, IST, SSD)
//  bool isPXL1 = (rawTrack.mHFThitsMap & PXL1_MASK) != 0;  // PXL1 has been hit
//  bool isPXL2 = (rawTrack.mHFThitsMap & PXL2_MASK) != 0;  // PXL2 has been hit
//  bool isIST  = (rawTrack.mHFThitsMap & IST_MASK) != 0;   // IST has been hit
//  bool isSSD  = (rawTrack.mHFThitsMap & SSD_MASK) != 0;   // SSD has been hit
//  if(!(isPXL1 && isPXL2 && (isIST || isSSD))) {return false;} 
//  hm->hTrackCuts->Fill(3.5);
  // DCA
  if( rawTrack.mDCA >= params->paramMap["dcaMax"][0] ) {return false;}
  if( firstCall ) {hm->hTrackCuts->Fill(3.5);}
  // nHits
  //if( rawTrack.mNhits < round(params->paramMap["nHitsMin"][0]) ) { return false; }
  //hm->hTrackCuts->Fill(4.5);
  // nHitsFit
  if( rawTrack.mNhitsFit < round(params->paramMap["nHitsFitMin"][0]) ) { return false; }
  if( firstCall ) {hm->hTrackCuts->Fill(4.5);}
  // nHitsDedx
  if( rawTrack.mNhitsDedx < round(params->paramMap["nHitsDedxMin"][0]) ) { return false; }
  if( firstCall ) {hm->hTrackCuts->Fill(5.5);}
  // nHitsFit/nHitsMax < 0.52
  if( (static_cast<double>(rawTrack.mNhitsFit) / static_cast<double>(rawTrack.mNhitsMax)) < 0.52 ) { return false; } //EVENTUALLY PUT THIS INTO CUT-TABLE.txt
  if( firstCall ) {hm->hTrackCuts->Fill(6.5);}
  
  //Phi cut (poor TPC sectors from Run16)
  //TrackStruct myTrack;
  //myTrack = rawTrack.mTrack;
  //double phi = myTrack.mPhi;
  //if( (phi >= -2.760 && phi <= -1.38) || (phi >= -1.26 && phi <= 0.13)
  //	|| (phi >= 2.375 && phi <= 3.02) || (phi >= 2.13 && phi <= 2.34) ) {return false; }
  if( firstCall ) {hm->hTrackCuts->Fill(7.5);}
  // p
//  if( rawTrack.mTrack.mP.Mag() < params->paramMap["pMin"][0] || rawTrack.mTrack.mP.Mag() > params->paramMap["pMax"][0] ) {
//    return false;
//  }
//
//  // pT
//  if( rawTrack.mTrack.mP.Perp() < params->paramMap["ptMin"][0] || rawTrack.mTrack.mP.Perp() > params->paramMap["ptMax"][0] ) {
//    return false;
//  }
//
//  // eta
//  if( rawTrack.mTrack.mEta < params->paramMap["etaMin"][0] || rawTrack.mTrack.mEta > params->paramMap["etaMax"][0] ) { return false; }
firstCall = false;
  return true;
}
//------------------------------------------------------------------------
void StUtilities::FillAnalyzerTracks(StHistograms *hm, double refMult, RawTrackStruct rawTrack, TH2D* hWeighting, double refMultWeight){
	TrackStruct myTrack;
        myTrack = rawTrack.mTrack;
	double eta = myTrack.mEta; double phi = myTrack.mPhi, weight{0.0};
        int bin = hWeighting->FindBin(eta, phi);
        //Note: if the weighting is still not working, consider the default binning being assigned to hWeighting inside StUtils (is it the same as the TH2D in StSoundVelocity?)
        if ( hWeighting->GetBinContent(bin) > 0 ) weight = 1.0 / ( hWeighting->GetBinContent(bin) ) ; //Note: trying 1/weight to enfore upweighting
 //cout << "Utils weight: " << weight << endl;
        if (weight > 0) weight = weight * refMultWeight;
		//weight = 1; //Temporary (for Run11)
		hm->hPtDCA[ 0 ]->Fill( rawTrack.mDCA, myTrack.mPt, weight ); /// <-- DCA acceptance here 
                hm->hPtEta[ 0 ]->Fill( myTrack.mEta, myTrack.mPt, weight );
                hm->hEta[ 0 ]->Fill(myTrack.mEta, weight );
                hm->hEtaPhi[ 0 ]->Fill( myTrack.mEta, myTrack.mPhi, weight ); /// <-- Eta acceptance filling here
                hm->hEta_Nhits[ 0 ]->Fill( myTrack.mEta, rawTrack.mNhits, weight );

                hm->hEta_Nch[ 0 ]->Fill( myTrack.mEta, refMult, weight );
		hm->hPhi_Nch[ 0 ]->Fill( myTrack.mPhi, refMult, weight );
                //hm->hPt_Nch_eta->Fill(refMult, myTrack.mPt, myTrack.mEta, weight );
                hm->hPt_RefMultCorr->Fill( refMult,myTrack.mPt, weight ); //For <pT> and dNdEta computation
		hm->hNch_NhitsFit[ 0 ]->Fill( refMult, rawTrack.mNhitsFit, weight );
                hm->hNch_NhitsdEdx[ 0 ]->Fill( refMult, rawTrack.mNhitsDedx, weight );
                hm->hNch_DCA[ 0 ]->Fill( refMult, rawTrack.mDCA, weight );
              //hm->hPt_dNdEta[ k+j ]->Fill( myTracks[k].size() / abs(eta_high_cut[eta_index] - eta_low_cut[eta_index]),myTrack.mPt );
              //hm->hChargePt[k+j]->Fill( myTrack.mCharge, myTrack.mPt ); /// <-- Pt of charged particles filling here 
}

//------------------------------------------------------------------------
int StUtilities::GetMultBin( int centBin1to16 ) {

  int mbin = -1;

  if( centBin1to16 == 16 ) { mbin = 0; }
  else if( centBin1to16 == 15 ) { mbin = 1; }
  else if( centBin1to16 == 14 ) { mbin = 2; }
  else if( centBin1to16 == 13 ) { mbin = 3; }
  else if( centBin1to16 == 12 ) { mbin = 4; }
  else if( centBin1to16 == 11 ) { mbin = 5; }
  else if( centBin1to16 == 10 ) { mbin = 6; }
  else if( centBin1to16 == 9 ) { mbin = 7; }
  else if( centBin1to16 == 8 ) { mbin = 8; }
  else if( centBin1to16 == 7 ) { mbin = 9; }
  else if( centBin1to16 == 6 ) { mbin = 10; }
  else if( centBin1to16 == 5 ) { mbin = 11; }
  else if( centBin1to16 == 4 ) { mbin = 12; }
  else { mbin = 13; }

  return mbin;
}

//------------------------------------------------------------------------
int StUtilities::GetVzBin( StParameters *params, double vz ) {

    //One Vz Bin for FXT
  int mbin = int( (vz - params->paramMap["VzMin"][0]) / ParamConsts::ZvtxWidthMix);

  return mbin;
}

//------------------------------------------------------------------------
int StUtilities::CalcRefMult( StPicoDst *pDst ) {

  int refMult = 0;

  if( mFXT ) {
    for(UInt_t i = 0; i < pDst->numberOfTracks(); i++ ) {

      StPicoTrack *picoTrack = pDst->track( i );
      if( ! picoTrack ) { continue; }
      if( ! picoTrack->isPrimary() ) { continue; }

      refMult++;
    }
  }
  else {

    refMult = pDst->event()->refMult();

  }

  return refMult;
}
//------------------------------------------------------------------------
int StUtilities::CalcRefMult2( StPicoDst *pDst ){
  int refMult2 = 0;

  if( mFXT ) { // This does not specify the outer region
    for(UInt_t i = 0; i < pDst->numberOfTracks(); i++ ) {

      StPicoTrack *picoTrack = pDst->track( i );
      if( ! picoTrack ) { continue; }
      if( ! picoTrack->isPrimary() ) { continue; }
    
      refMult2++;
    }
  } 
  else {

    refMult2 = pDst->event()->refMult2();

  }

  return refMult2;
}
//------------------------------------------------------------------------
int StUtilities::CalcRefMult( StMuDst *muDst ) {

  int refMult = 0;

  if( mFXT ) {
    refMult = muDst->numberOfPrimaryTracks();
  }
  else {
    refMult = muDst->event()->refMult();
  }

  return refMult;
}

//------------------------------------------------------------------------
bool StUtilities::BadRun( StParameters *params, int run ) {

  for( UInt_t i = 0; i < params->mBadRunList.size(); i++ ) {
    if( run == params->mBadRunList[ i ] ) { return true; }
  }

  return false;
}

//------------------------------------------------------------------------
void StUtilities::Shuffle(Int_t *iarr, Int_t i1, Int_t i2){
  Int_t j, k;
  Int_t a = (i2 - i1) - 1;
  for (Int_t i = i1; i < i2; i++) {
    j = (Int_t) (gRandom->Rndm() * a);
    k = iarr[j];
    iarr[j] = iarr[i];
    iarr[i] = k;
  }
}

//------------------------------------------------------------------------
void StUtilities::SetiTPCCase( int run ){
  miTPC_In = (run >= ParamConsts::FirstRunWithiTPC) ? true : false;
}

//------------------------------------------------------------------------
//Unncessary function but runAnalysis breaks if not included:
void StUtilities::LoadFSIfile( TString filename ) {

  TFile *f = new TFile( filename, "READ" );

  mFSI_ss = (TH2D*)f->Get("K2_ss");
  mFSI_os = (TH2D*)f->Get("K2_os");
  if( mFSI_ss ) { mFSI_ss->SetDirectory(0); }
  if( mFSI_os ) { mFSI_os->SetDirectory(0); }

  f->Close();
}

//------------------------------------------------------------------------
void StUtilities::PrintTrackInfo( RawTrackStruct rawtrack ){
  cout<<"DCA:"<<rawtrack.mDCA<<"  nHits:"<<rawtrack.mNhits<<"  nHitsDedx:"<<rawtrack.mNhitsDedx<<"  p:"<<rawtrack.mTrack.mP.Mag()<<"  pT:"<<rawtrack.mTrack.mP.Perp()<<"  eta:"<<rawtrack.mTrack.mEta<<"  nSigPi:"<<rawtrack.mNsigmaPion<<"  nSigEl:"<<rawtrack.mNsigmaElectron<<"  nSigK:"<<rawtrack.mNsigmaKaon<<"  nSigPr:"<<rawtrack.mNsigmaProton<<"  TOFidx:"<<rawtrack.mBTOFindex<<"  matchFlag:"<<rawtrack.mBTOFmatchFlag<<"  ylocal:"<<rawtrack.mBTOFylocal<<"  zlocal:"<<rawtrack.mBTOFzlocal<<"  m2:"<<rawtrack.mM2<<"  beta:"<<rawtrack.mBeta<<"  betaPion:"<<rawtrack.mBetaPion<<endl;
}
//------------------------------------------------------------------------
bool StUtilities::passnTofMatchRefmultCut(Double_t refmult, Double_t ntofmatch, bool window1, int etaIndex){
	// window1 is centrality flag
	double x_cutoff{0};
	bool notPileUp{false};
	//Double_t a0{}, a1{}, a2{}, a3{}, a4{};
        Double_t b0{}, b1{}, b2{}, b3{}, b4{};
        Double_t c0{}, c1{}, c2{}, c3{}, c4{};
        Double_t refmultcutmax{};
        Double_t refmultcutmin{};

	// 200 GeV Run 16 (2016) dataset (pTlow0p15, Eta0p5, dca3, Nhits > 15):
	//b0 = 9.87558635118469; 
	//b1 = 2.50070701452315;
	//b2 = -0.00780922650423221;
	//b3 = 3.76985911388424e-05;
	//b4 = -5.60999307825242e-08;
	//c0 = -13.760653149845;
	//c1 = 1.3339416056808;
	//c2 =  -0.00257814142758717;
	//c3 = 2.2003477167024e-05;
	//c4 = -3.80153554499135e-08;	

	    // 200 GeV Run 16 (2016) dataset (pTlow0p20, Eta0p5, dca2,|Vz| < 6, Vr < 2.0,Nhits >= 15, NhitsFit >= 15, NhitsdEdx > 10, Self-correlation):
        //if (!Outer){
        // b0 = 5.85177120016888;
        // b1 = 2.66485814250352;
        // b2 = -0.0117415885664514;
        // b3 = 6.09079135569088e-05;
        // b4 = -9.93241426465257e-08;
        // c0 = -13.8316939526706;
        // c1 = 1.22093472348819;
        // c2 = -0.0028917554569357;
        // c3 = 2.76182456085501e-05;
        // c4 = -5.3285351307685e-08;
        //}

     // 200 GeV Run 16 (2016) dataset (pTlow0p20, Eta0p5, dca2,|Vz| < 6, Vr < 2.0,Nhits >= 15, NhitsFit >= 15, NhitsdEdx > 10, Self-correlation):
        //if (Outer){
        // b0 = 10.0114288670218;
        // b1 = 2.84977730988112;
        // b2 = -0.00924854195875469;
        // b3 = 5.25377391875824e-05;
        // b4 = -1.012450389272e-07;
        // c0 = -12.116034361624;
        // c1 = 1.30111712686582;
        // c2 = 3.35516852019208e-05;
        // c3 = 1.68922512168319e-05;
        // c4 = -4.84700279878035e-08;
        //} 	
	
	// 200 GeV Run 16 (2016) dataset (pTlow0p20, Eta0p5, dca2, Nhits >= 15, NhitsFit >= 15, NhitsdEdx > 10, Self-correlation):
	//if (!Outer){
	// b0 = 6.58034859679174;
	// b1 = 2.39814922547886;
	// b2 = -0.00828628943723686;
	// b3 = 4.36530058077724e-05;
	// b4 = -7.19278168932131e-08;
	// c0 = -11.3561841461277;
	// c1 = 1.21234684586826;
	// c2 = -0.00153758226430252 ;
	// c3 = 1.93971329302196e-05;
	// c4 = -3.93007209768407e-08;
	//}

	// 200 GeV Run 16 (2016) dataset (pTlow0p20, EtaOuter0p5, dca2, Nhits >= 15, NhitsFit >= 15, NhitsdEdx > 10, Seperate):
	//if (Outer){
        // b0 = 9.21407373818314;
        // b1 = 2.65507533480053;
        // b2 = -0.00718215859573945;
        // b3 = 4.12673635671513e-05;
        // b4 = -8.04404009597605e-08;
        // c0 = -10.0534490323685;
        // c1 = 1.33475494406829;
        // c2 = 0.000577782322061653;
        // c3 = 1.23163901100102e-05;
        // c4 = -3.88511289180943e-08;
	//}
	
	// 200 GeV Run 16 (2016) dataset (StRefMultCorr, gRefMult):
         //b0 = 6.3403945998022;
         //b1 = 1.77996736031512;
         //b2 = -0.0011631293389005;
         //b3 = 1.87446630846821e-05;
         //b4 = -3.84341923743788e-08;
         //c0 = -10.3361319471878;
         //c1 = 1.44381447703204;
         //c2 = -0.00357054330548962;
         //c3 = 2.92642368635652e-05;
         //c4 = -5.01265033040289e-08;
  
	 //200 GeV Run 16 (2016) dataset (StRefMultCorr, StRefMultCorr weights, Trigger IDs: 520001,520011,520021,520031,520041,520051)
	 b0 = 6.68427162711842;
	 b1 = 2.96924103162064;
	 b2 = -0.012065602533994;
	 b3 = 5.57565034339225e-05;
	 b4 = -8.08207787078094e-08;
	 c0 = -12.1560731653035;
         c1 = 1.41187314508252 ;
         c2 = -0.0023446615654021;
         c3 = 2.01269384281336e-05;
         c4 = -3.43722371433138e-08;

        refmultcutmax = calcPileUpRefMult(ntofmatch, b0, b1, b2, b3, b4); // --- Check these definitions if the cuts aren't working
        refmultcutmin = calcPileUpRefMult(ntofmatch, c0, c1, c2, c3, c4); // ---
            if ( isInPileUpRefMultLimits( refmult, refmultcutmin, refmultcutmax ) ) {notPileUp = true;}
	    if (x_cutoff != 0) { if (ntofmatch > x_cutoff) {notPileUp = false;} } 
	//if (!notPileUp) {cout << "\n refmult: " << refmult << " | ntofmatch: " << ntofmatch << " | refmultcutmax: " << refmultcutmax << " | refmultcutmin: " << refmultcutmin << " | notPileUp: " << notPileUp; }    
	return notPileUp;

} 
//------------------------------------------------------------------------
Double_t StUtilities::calcPileUpRefMult(Double_t ntofmatch, Double_t x0, Double_t x1, Double_t x2, Double_t x3, Double_t x4) {
    return ( x0 + x1*(ntofmatch) + x2*pow(ntofmatch,2) + x3*pow(ntofmatch,3) + x4*pow(ntofmatch,4) );
}
//------------------------------------------------------------------------
bool StUtilities::isInPileUpRefMultLimits(Double_t refmult,Double_t refmultcutmin,Double_t refmultcutmax ){
    bool inLimits{false};
    if (refmult > refmultcutmin && refmult < refmultcutmax) {inLimits = true;}
    return inLimits;
}
StUtilities::~StUtilities()
{
}

#endif

