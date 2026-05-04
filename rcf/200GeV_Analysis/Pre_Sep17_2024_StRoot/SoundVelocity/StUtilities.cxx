#ifndef STUTILITIES_CXX
#define STUTILITIES_CXX

#include "TH1F.h"
#include "TH2F.h"

#include "StUtilities.h"
#include "StParameters.h"

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
void StUtilities::FillPicoRawTracks( StPicoDst *picoDst, StParameters *params, vector<RawTrackStruct> *tracks ) {

  TVector3 priVtx = picoDst->event()->primaryVertex();

  SetiTPCCase( picoDst->event()->runId() );

  for( UInt_t i = 0; i < picoDst->numberOfTracks(); i++ ) {
    RawTrackStruct myRawTrack;

    StPicoTrack *pTrack = picoDst->track( i );
    if( ! pTrack ) continue;

    TVector3 momentum= pTrack->pMom();
    
    //myRawTrack.mTrack.mFlag = pTrack->flag(); // ??
    
    myRawTrack.mTrack.mCharge = pTrack->charge();
    myRawTrack.mTrack.mFirstWordTopologyMap = pTrack->topologyMap(0);
    myRawTrack.mTrack.mSecondWordTopologyMap = pTrack->topologyMap(1);
    myRawTrack.mTrack.mThirdWordTopologyMap = pTrack->iTpcTopologyMap();
    myRawTrack.mTrack.mP = momentum;
    myRawTrack.mTrack.mPt = momentum.Perp();
    myRawTrack.mTrack.mPhi = momentum.Phi();
    myRawTrack.mTrack.mEta = momentum.PseudoRapidity();
    myRawTrack.mTrack.mMass = constants::mass_pion; // pion mass assumption
    myRawTrack.mTrack.mE = sqrt( momentum.Mag2() + pow(myRawTrack.mTrack.mMass,2) );

    myRawTrack.mRawTrackID = pTrack->id();
    myRawTrack.mDCA = pTrack->gDCA(priVtx.x(), priVtx.y(), priVtx.z());
    myRawTrack.mNhits = pTrack->nHits();
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

//    double QB = myRawTrack.mTrack.mCharge * picoDst->event()->bField(); // B needs to be in kiloGauss
//    CalculateCrossings( params, momentum, QB, myRawTrack.mTrack.mCrossings );
//
//    // MC matching (only fully implemented for MuDsts)
//    myRawTrack.mTrack.mID_MC = pTrack->idTruth();
//    myRawTrack.mTrack.mP_MC = TVector3(0,0,0);

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
bool StUtilities::GoodPicoEvent( StParameters *params, StPicoEvent *event ) {

  TVector3 priVtx = event->primaryVertex();
  float VpdVz = event->vzVpd();
  /// Trigger
  bool goodTrigger = false;
  for( uint i = 0; i < params->paramMap["Trigger"].size(); i++ ) {
    if( event->isTrigger( round(params->paramMap["Trigger"][i]) ) ) { goodTrigger = true; }
  }
  if( ! goodTrigger ) { return false; }

  /// Vz
  if( priVtx.z() < params->paramMap["VzMin"][0] ) { return false; }
  if( priVtx.z() > params->paramMap["VzMax"][0] ) { return false; }
  /// Vz - Vpd
  /// |Vz_{TPC} - Vz_{VPD}| < cut
  if (priVtx.z() - VpdVz > params->paramMap["VzVpdMax"][0] || priVtx.z() - VpdVz < params->paramMap["VzVpdMin"][0]) {return false;}  
//cout<<"past Vz"<<endl;
  /// Vr
  double VyShift = 0;
  //if( mFXT ) {
  //  VyShift = -2; // -2 cm Y offset of fixed target
  //}
  double Vr = sqrt( pow( priVtx.x(), 2) + pow( priVtx.y() - VyShift, 2) );
  if( Vr > params->paramMap["VrCut"][0] ) { return false; }
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
bool StUtilities::GoodTrack( StParameters *params, RawTrackStruct rawTrack ) {

  // DCA
  if( rawTrack.mDCA > params->paramMap["dcaMax"][0] ) {return false;}

  // nHits
  if( rawTrack.mNhits < round(params->paramMap["nHitsMin"][0]) ) { return false; }

  // nHitsDedx
  //if( rawTrack.mNhitsDedx <= round(params->paramMap["nHitsDedxMin"][0]) ) { return false; }

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

  return true;
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
StUtilities::~StUtilities()
{
}

#endif

