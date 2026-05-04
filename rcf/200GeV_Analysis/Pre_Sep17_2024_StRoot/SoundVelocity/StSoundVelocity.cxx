/// Notes: Removing eta and pT cuts from cutTable and enforcing them here
#ifndef STMULTIPIONFEMTO_CXX
#define STmULTIPIONFEMTO_CXX
/// NOTES: Need to add p_{T} distribution, eta distribution
///                   ^ for the purpose of creating a 2D histogram (Nch, pt) for the purpose of doing Omar's analysis
#include "TH1F.h"
#include "TH2F.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"

#include "StRoot/StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StRoot/StMuDSTMaker/COMMON/StMuDst.h"
#include "StRoot/StMuDSTMaker/COMMON/StMuEvent.h"
#include "StRoot/StMuDSTMaker/COMMON/StMuTrack.h"
#include "StRoot/StMuDSTMaker/COMMON/StMuMcTrack.h"
#include "StRoot/StMuDSTMaker/COMMON/StMuMcVertex.h"
#include "StRoot/StMuDSTMaker/COMMON/StMuPrimaryVertex.h"

#include "StRefMultCorr/CentralityMaker.h"
#include "StRefMultCorr/StRefMultCorr.h"

#include "StSoundVelocity.h"
#include "StHistograms.h"
#include "StParameters.h"
#include "StEventCollection.h"
#include "StUtilities.h"

ClassImp(StSoundVelocity)

  //------------------------------------------------------------------------
StSoundVelocity::StSoundVelocity( TString setting, bool picoCase ) : StMaker()
{
  /// setting at the moment is "_200p0_COL"
  mSetting = setting;
  mPicoCase = picoCase;
}

//------------------------------------------------------------------------
StSoundVelocity::~StSoundVelocity()
{
}

//------------------------------------------------------------------------
Int_t StSoundVelocity::Init() //Initialization
{
  mDebug = 0;

  //cout<<"StSoundVelocity Init()"<<endl;
//  cout<<"2-pion only? "<<mTwoPionOnly<<endl;
  // All the Utility functions
//  utils = new StUtilities( mSetting, mKfile, mMRCfile, mWeightfile );
  utils = new StUtilities( mSetting );
  // reference multiplicity class
  if( utils->mFXT ) { 
    refMultCorr = CentralityMaker::instance()->getRefMultCorrFxt(); /// Fixed target mode
      /// Function description: https://www.star.bnl.gov/webdata/dox/html/classStRefMultCorr.html
      /// refMultCorr is initialized in StSoundVelocity.h as an StRefMultCorr object
  }
  else {
    refMultCorr = CentralityMaker::instance()->getRefMultCorr();    /// Collider mode (Throws error(?))
      /// Funciton description: https://www.star.bnl.gov/webdata/dox/html/classStRefMultCorr.html
      /// refMultCorr is initialized in StSoundVelocity.h as an StRefMultCorr object
  }

  // Parameters
  params = new StParameters( mSetting );
  if( ! params->LoadParamTable() ) { return kStFatal; } /// <- Error status (couldn't LoadParamTable())

  int ZvtxMixBins = round( (params->paramMap["VzMax"][0] - params->paramMap["VzMin"][0]) / ParamConsts::ZvtxWidthMix );
  if( ZvtxMixBins == 0 ) { ZvtxMixBins++; } // FXT will have only 1 bin

  // Histogram manager
  hm = new StHistograms( params );

  //cout<<"Done StSoundVelocity Init()"<<endl;
  return kStOK; //<- STATUS: OK
}

//------------------------------------------------------------------------
Int_t StSoundVelocity::Make()
{
    static const int kNumEtaRanges = 4;
    static const int kNumpTRanges = 4;
    double eta_cut[kNumEtaRanges] = {0.5, 0.8, 1.0, 0.65};
    double pt_low_cut[kNumpTRanges] = {0.0, 0.15, 0.4, 0.8};
    double eta; // Pulling eta from mytrack.Eta
    double pt;// double p;
    //double pt_low = 1.e-10; //Same/similar cuts as RefMult from StRefMultCorr (StuRefMult.hh)
    //double p_low = 1.e-10; //Same/similar cuts as RefMult from StRefMultCorr (StuRefMult.hh)
	// After Apr 3^
   //Before Apr 3:  double pt_low = 0.15; double pt_high = 10; /// 0.15 and 10 GeV/c resp
   //Before Apr 3:   double p_low = 0.15; double p_high = 10; /// 0.15 and 10 GeV/c resp
    
  StPicoDstMaker *picoDstMaker;
  StPicoDst *picoDst;
  StPicoEvent *picoEvent;

  StMuDstMaker *muDstMaker;
  StMuDst *muDst;
  StMuEvent *muEvent;

  TVector3 priVtx;
  int refMult = 0; 
  int runID = 0;
  //int eventID = 0;
  float ZDCx = 0;
  int nBTOFmatch = 0;
    
  // Select good events

  if( mPicoCase ) { // PicoDst 

    picoDstMaker = (StPicoDstMaker*) GetMaker("picoDst");
    if( ! picoDstMaker ) { return kStWarn; }
    picoDst = picoDstMaker->picoDst();
    if( ! picoDst ) { return kStWarn; }
    picoEvent = picoDst->event();
    if( ! picoEvent ) { return kStWarn; }

    if( ! utils->GoodPicoEvent( params, picoEvent ) ) { return kStOK; }

    priVtx = picoEvent->primaryVertex();
    refMult = utils->CalcRefMult( picoDst );
    //eventID = picoEvent->eventId();
    runID = picoEvent->runId();
    ZDCx = picoEvent->ZDCx();
    nBTOFmatch = picoEvent->nBTOFMatch();
    //  cout << "PicoDST Done" << endl;
  }
  else { // MuDst

    muDstMaker = (StMuDstMaker*) GetMaker("MuDst");
    if( ! muDstMaker ) { return kStWarn; }
    muDst = muDstMaker->muDst();
    if( ! muDst ) { return kStWarn; }
    muEvent = muDst->event();
    if( ! muEvent ) { return kStWarn; }

    if( ! utils->GoodMuEvent( params, muEvent ) ) { return kStOK; }

    StThreeVectorF Stvtx = muDst->primaryVertex()->position(); 
    priVtx.SetXYZ( Stvtx.x(), Stvtx.y(), Stvtx.z() );
    refMult = utils->CalcRefMult( muDst );
    //eventID = muEvent->eventId();
    runID = muEvent->runId();
    ZDCx = 0;
    nBTOFmatch = 0;
  }

//  if( mDebug == 1 ) { LOG_INFO << "PASSED MY EVENT SELECTION" <<endm; }
  
  // Vz bin (ETA CUT HERE?)
  //int VzIdx = utils->GetVzBin( params, priVtx.z() );
  hm->hVertex->Fill( priVtx.x(), priVtx.y(), priVtx.z() );
  //cout << "Vz bin filled" << endl;

  // StRefMultCorr initialize
  if( ! refMultCorr ) { return kStWarn; }
  refMultCorr->init( runID );
  //cout << "RefMultCorr initialized" << endl;
    
  // Bad run cut (Removes bad runs)
  if( refMultCorr->isBadRun( runID ) ) { return kStOK; }  
    refMultCorr->initEvent( refMult, priVtx.z(), ZDCx );
    // Pile-up cut (Pileup rejection here)
  if( ! refMultCorr->passnTofMatchRefmultCut( double(refMult), double(nBTOFmatch) ) ) { return kStOK; }
  
  hm->hMultCent->Fill( refMult, 1 + refMultCorr->getCentralityBin16() ); /// use getCentralityBin16 instead?
  hm->hRefMult->Fill ( refMult );  
//cout << '\n' << "refMult: " << refMult << endl;
// Multiplicity bin
//commented out b/c I don't want to discriminate centralities
//  int multIdx = utils->GetMultBin( 1 + refMultCorr->getCentralityBin9() ); /// use getCentralityBin16 instead?
//  if( multIdx < 0 ) {
//    return kStOK;
//  }
//    cout << "GetRefMult bin = " << multIdx << endl;
  vector<RawTrackStruct> myRawTracks;

  if( mPicoCase ) {
    utils->FillPicoRawTracks( picoDst, params, &myRawTracks );
  }
  else {
    utils->FillMuRawTracks( muDst, params, &myRawTracks );
  }
  // ETA CUT HERE?
  hm->hMultRaw->Fill( myRawTracks.size() );
  //cout << '\n' << "MultRaw Filled" << endl;
  ////////////////////////////////////////
  // Select good tracks
  //
    std::vector<TrackStruct> myTracks[16]; //Array of vectors, see table at bottom for index renferences
    
//    vector<TrackStruct> myTracks;
//    vector<TrackStruct> myTracks0; /// eta0p5
//    vector<TrackStruct> myTracks1; /// eta0p8
//    vector<TrackStruct> myTracks2; /// eta1p0
//    vector<TrackStruct> myTracks3; /// et
    

  // Shuffle indices to shuffle tracks
  const int nRawTracks = myRawTracks.size();
  Int_t randomIndices[ nRawTracks ];
  for( Int_t i = 0; i < nRawTracks; i++ ) { randomIndices[ i ] = i; }
  utils->Shuffle( randomIndices, 0, nRawTracks );
	//cout << '\n' << "Selecting good tracks";
  // store good tracks
  for( Int_t i = 0; i < nRawTracks; i++ ) {
    TrackStruct myTrack;
    RawTrackStruct rawTrack = myRawTracks[ randomIndices[ i ] ];
    if( ! utils->GoodTrack( params, rawTrack ) ) { continue; }

    myTrack = rawTrack.mTrack;

    //ETA CUT APPLIED
      eta = myTrack.mEta;
      pt = myTrack.mPt;
  //    p = myTrack.mP.Mag();
      for (int k = 0; k < 13; k+=4) {
	int eta_index = 0;
	if (k == 4){eta_index = 1;}
	else if (k == 8) {eta_index = 2;}
	else if (k == 12) {eta_index = 3;}
 // Event selection based on centrality criteria (eta and pT)       
  if (std::abs(eta) <= eta_cut[eta_index]){
	//cout << '\n' << "INSIDE ETA CUT LOOP eta: " << eta_cut[k];
              for (int j = 0; j < 4; j++) {
                if (pt >= pt_low_cut[j]) {
		//	cout << '\n' << "============================";
		//	cout << '\n' << "low pT cut: " << pt_low_cut[j];
 		//	cout << '\n' << "k = " << k << ", j = " << j;
		//	cout << '\n' << "pT loop k+j = " << ( k + j );
  //                    hm->hChargePt[k+j]->Fill( myTrack.mCharge, myTrack.mPt ); /// <-- Pt of charged particles filling here (Seg. Violation here)
  //                    hm->hTransvMomentum[k+j]->Fill(myTrack.mPt); /// <-- Mean p_{T} filling here
                      hm->hEta[k+j]->Fill(myTrack.mEta); /// <-- Exclusive eta acceptance filling here
                      hm->hEtaPhi[k+j]->Fill( myTrack.mEta, myTrack.mPhi ); /// <-- Eta acceptance filling here
                      hm->hPtDCA[k+j]->Fill( myTrack.mPt, rawTrack.mDCA ); /// <-- DCA acceptance here (Seg. Violation here)
                      myTracks[k+j].push_back( myTrack ); //This is my Nch++
		      hm->hPtEta[k]->Fill( myTrack.mPt, myTrack.mEta );
                  }
              }
          }
      }
  }
//Filling Multiplicity Histogram (events)
    for (int i = 0; i < 16; i++) {
        if ( myTracks[i].size() > ParamConsts::MaxMult ){
            return kStOK;
        }
        ///ETA CUT characteristic saved inside myTracks[X] (X = 0,1,2,3,...)
	//cout << '\n' << "myTracks[" << i << "].size(): "<< myTracks[i].size();
        hm->hMult[i]->Fill( myTracks[i].size() ); /// <--Nch filling her
}
  
     for( Int_t i = 0; i < nRawTracks; i++ ) { randomIndices[ i ] = i; }
     utils->Shuffle( randomIndices, 0, nRawTracks );
    ///cout << '\n' << "Filling 2D Histograms";
    ///Track selection for generic criteria (pT selection for eta and pT)
    for (int k = 0; k < 16; k++) {
        for( Int_t i = 0; i < nRawTracks; i++ ){ //CHANGES MADE HERE AND EVERYTHING WITH "2" AT THE END IS FOR THE SECOND TRACK-BASED FOR LOOP THAT FULLFILLS THE TRACK FILLING CONDITION (AS OPPOSED TO THE CENTRALITY QUALIFIER CONDITIONS IN THE FOR LOOP ABOVE)
	    TrackStruct myTrack;
            RawTrackStruct rawTrack = myRawTracks[ randomIndices[ i ] ];
            if( ! utils->GoodTrack( params, rawTrack ) ) { continue; }

            myTrack = rawTrack.mTrack;
            
//            TrackStruct myTrack =  myTracks[k][i]; //<==== Check this (mistake here) THIS IS COMMENTED OUT?
            if (myTrack.mEta <= 0.5){
		if ( myTrack.mPt <= 10 && myTrack.mPt >= 0.15) {
                hm->hPt_Nch[k]->Fill( myTracks[k].size(),myTrack.mPt );
                //usage of refMult correct?
  //              hm->hPt_RefMult[k]->Fill( refMult, myTrack.mPt );
//                hm->hPt_dNdEta[k]->Fill( myTracks[k].size(),myTrack.mPt ); // eta gap is 1 so no scaling for per unit eta
            }
	  }
        }
    }
   //cout << '\n' << "End of StSoundVelocity.cxx Init()";
  return kStOK;
}
#endif
//index | Eta | pT_low (centrality qualifiers)
//0     | 0.5 | 0
//1     | 0.5 | 0.15
//2     | 0.5 | 0.4
//3     | 0.5 | 0.8
//4     | 0.8 | 0
//5     | 0.8 | 0.15
//6     | 0.8 | 0.4
//7     | 0.8 | 0.8
//8     | 1.0 | 0
//9     | 1.0 | 0.15
//10    | 1.0 | 0.4
//11    | 1.0 | 0.8
//12    | 0.65| 0
//13    | 0.65| 0.15
//14    | 0.65| 0.4
//15    | 0.65| 0.8
