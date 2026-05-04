// 200 GeV SoundVelocity Analysis Maker written by Caleb Broodo
/// Notes: Removing eta and pT cuts from cutTable and enforcing them here
#ifndef STMULTIPIONFEMTO_CXX
#define STmULTIPIONFEMTO_CXX

#include "TH1F.h"
#include "TH2F.h"
#include "TSystem.h"
#include "TChain.h"
#include "TChainElement.h"

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
  /// setting at the moment is "_3p0_FXT"
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
  //  refMultCorr = CentralityMaker::instance()->getRefMultCorrFxt(); /// Fixed target mode
      /// Function description: https://www.star.bnl.gov/webdata/dox/html/classStRefMultCorr.html
      /// refMultCorr is initialized in StSoundVelocity.h as an StRefMultCorr object
  }
  else {
        //refMultCorr = CentralityMaker::instance()->getRefMultCorr();    /// Collider mode (Throws error(?))
	refMultCorr = CentralityMaker::instance()->getgRefMultCorr_P16id();
      /// Funciton description: https://www.star.bnl.gov/webdata/dox/html/classStRefMultCorr.html
      /// refMultCorr is initialized in StSoundVelocity.h as an StRefMultCorr object
  }
	//std::cout << "GetChain() result: " << GetChain() << std::endl;

  //cout << "Inside Init() for StSoundVelocity" << endl;
  StChain* chain = dynamic_cast<StChain*>(GetChain());
  StPicoDstMaker* picoDstMaker = dynamic_cast<StPicoDstMaker*>(chain->GetMaker("picoDst"));
  
  if (picoDstMaker) {
  TChain* picoChain = picoDstMaker->chain();
    if (picoChain) {
        // Get the list of input files
        TObjArray* fileList = picoChain->GetListOfFiles();
        if (fileList && fileList->GetEntries() > 0) {
            for (int i = 0; i < fileList->GetEntries(); ++i) {
		TChainElement* fileEntry = static_cast<TChainElement*>(fileList->At(i));
		if (!fileEntry) continue;
		std::string fileName = fileEntry->GetTitle();
                //std::cout << "Input file #" << i + 1 << ": " << fileName << std::endl;
            if (fileName.find("st_sst") != std::string::npos) {
                mStream = 1;
                //std::cout << "Stream: st_sst" << std::endl;
            } else if (fileName.find("st_nosst") != std::string::npos) {
                mStream = 1;
		//std::cout << "Stream: st_nosst" << std::endl;
            } else if (fileName.find("st_physics") != std::string::npos) {
                mStream = 0;
		//std::cout << "Stream: st_physics" << std::endl;
            } else {
                return kStWarn;
		//std::cout << "Unknown stream in file name: " << fileName << std::endl;
            } 
	}
    } 
  }
}  

    refMultCorr->setVzForWeight(6, -6.0, 6.0); //Vz cut    
    if (mStream == 0) refMultCorr->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_VpdnoVtx_Vpd5_Run16.txt"); 
    //for new StRefMultCorr, Run16, SL16j, production 1, physics stream
    if (mStream == 1) refMultCorr->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_VpdnoVtx_Vpd5_Run16_prod2.txt"); 
    //for new StRefMultCorr, Run16, SL16j, production 2, st_nosst and st_sst stream


    std::ifstream inFile("./StRoot/SoundVelocity/weighting_txt.txt");
    std::string line;
    while (std::getline(inFile, line)) {
        std::istringstream iss(line);
        int binX, binY; double xCoord, yCoord, content;
        iss >> xCoord >> yCoord >> content;
        //cout << "  " << xCoord << " | " << yCoord << " | " << content << endl;
        // Fill the histogram
        binX = hWeighting->GetXaxis()->FindBin(xCoord);
        binY = hWeighting->GetYaxis()->FindBin(yCoord);
	//content = 1; //Removing weighting to see effect in final SOS plot on dN/deta
        hWeighting->SetBinContent(binX, binY, content);
    }
 

  //int bin = hWeighting->FindBin(0.1, 1.5);
  //double weight = hWeighting->GetBinContent(bin);
  //cout << "eta: 0.1 phi: 1.5 and bin: " << bin << endl;
  //cout << "TESTING Weighting histo reading -  weight: " << weight << endl;
  

  params = new StParameters( mSetting );
  if( ! params->LoadParamTable() ) { return kStFatal; } /// <- Error status (couldn't LoadParamTable())
  params->SetBadRuns(); //Sets manually inserted Bad Runs (Make sure the energy/Run config is uncommented in StParameters.cxx)
  int ZvtxMixBins = round( (params->paramMap["VzMax"][0] - params->paramMap["VzMin"][0]) / ParamConsts::ZvtxWidthMix );
  if( ZvtxMixBins == 0 ) { ZvtxMixBins++; } // FXT will have only 1 bin
  //cout<<"Initializing Histogram Object"<<endl;
  // Histogram manager
  hm = new StHistograms( params );

  //cout<<"Done StSoundVelocity Init()"<<endl;
  return kStOK; //<- STATUS: OK
}

//------------------------------------------------------------------------
Int_t StSoundVelocity::Make()
{
  //cout << "Inside SoundVel Make()" << endl;
    double eta; // Pulling eta from mytrack.Eta
    double pt;// double p;
    
  StPicoDstMaker *picoDstMaker;
  StPicoDst *picoDst;
  StPicoEvent *picoEvent;

  StMuDstMaker *muDstMaker;
  StMuDst *muDst;
  StMuEvent *muEvent;

  TVector3 priVtx;
  double refMult = 0;
  double refMult2 = 0; 
//  Int_t Nch_NoCuts{0};
//  Int_t nBTofMatch_NoCuts{0};
  int runID = 0;
  double ZdcEast{0.0};
  double ZdcWest{0.0};
  //int eventID = 0;
  float ZDCx = 0;
  float BBCx{0};
  double Vz{0};
  //int nBTOFmatch = 0;
  hm->hEventCuts->Fill(0.5);  
  // Select good events
  // cout << "Filled initial event cuts histo" << endl;
  if( mPicoCase ) { // PicoDst 
  //   cout << "creating PicoDST maker" << endl;
    picoDstMaker = (StPicoDstMaker*) GetMaker("picoDst");
    if( ! picoDstMaker ) { return kStWarn; }
    picoDst = picoDstMaker->picoDst();
    if( ! picoDst ) { return kStWarn; }
    picoEvent = picoDst->event();
    if( ! picoEvent ) { return kStWarn; }
    if( ! utils->GoodPicoEvent( params, picoEvent, hm ) ) { return kStOK; }
    priVtx = picoEvent->primaryVertex();
    refMult = utils->CalcRefMult( picoDst );
    //refMult = picoEvent->grefMult(); //Using grefMult
    refMult2 = utils->CalcRefMult2( picoDst ); //Need to find a way to get refMult2 via StRefMultCorr (?) so you can compare apples to apples
    //eventID = picoEvent->eventId();
    runID = picoEvent->runId();
  //  cout << "RunID: " << runID << endl;
    ZDCx = picoEvent->ZDCx();
    BBCx = picoEvent->BBCx();
    Vz = priVtx.z();
    ZdcEast = picoEvent->ZdcSumAdcEast();
    ZdcWest = picoEvent->ZdcSumAdcWest();
    //nBTOFmatch = picoEvent->nBTOFMatch();
    // cout << "PicoDST Done" << endl;
  }
  else { // MuDst
//	cout << "creating MuDST maker" << endl;
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
    BBCx = 0;
    ZDCx = 0;
    //nBTOFmatch = 0;
//	cout << "MuDST Done" << endl;
  }

  if( ! refMultCorr ) { return kStWarn; }
  refMultCorr->init( runID ); // StRefMultCorr initialize
  refMultCorr->initEvent( refMult, priVtx.z(), ZDCx );
  if( refMultCorr->isBadRun( runID ) ) { return kStOK; } // Bad run cut (Removes bad runs)
  hm->hEventCuts->Fill(6.5);
  //refMult = refMultCorr->getRefMultCorr(); //refMult definition here
  //refMultCorr->setVerbose(true);
  //refMultCorr->print();


 //Run16 Run rejection (filter via StUtils eventually). Link: https://drupal.star.bnl.gov/STAR/system/files/D0AnaNote_Feb21_2019.pdf#page=62&zoom=100,113,700
//  if ( ( runID < 17062047 ) || ( runID >=17065002 && runID <=17068053 ) || ( runID >= 17063041 && runID <= 17063043 ) ) { return kStOK; }
  if ( runID > 17062000 ) { return kStOK; }
  
  hm->hEventCuts->Fill(7.5);

  if( params->isBadRun( runID ) ) { return kStOK; } //manual bad run cut (see StParameters.cxx)  
  double Correctedrefmult = refMultCorr->getRefMultCorr();   
  double refMultWeight = refMultCorr->getWeight();
  //if( ! refMultCorr->passnTofMatchRefmultCut( double(refMult), double(nBTOFmatch) ) ) { return kStOK; } // Pile-up cut (Pileup rejection here)
  //if( refMultCorr->isPileUpEvent( double(refMult), double(nBTOFmatch) ) ) { return kStOK; } // Pile-up cut (Pileup rejection here)
  hm->hEventCuts->Fill(8.5); 
  hm->hVertex->Fill( priVtx.x(), priVtx.y(), priVtx.z() );
  vector<RawTrackStruct> myRawTracks;
  if( mPicoCase ) {
    utils->FillPicoRawTracks( picoDst, params, &myRawTracks, hm );
  }
  else {
    utils->FillMuRawTracks( muDst, params, &myRawTracks );
  }

  Int_t centrality{refMultCorr->getCentralityBin16()};
  //bool bin_0_to_10{centrality == 15 || centrality = 14}; //Flag to QA plots in 0-10% 
  Int_t centrality9{refMultCorr->getCentralityBin9()};

  std::vector<TrackStruct> myTracks[2]; //Array of vectors, see table at bottom for index renferences
//ETA INDEX SET AND PILEUP REMOVAL BOOLEAN HERE (Don't forget to adjust etaIndex in StUtils for pileup!):       
    int etaIndex{0}; // Eta Index is SET. specifies the eta config for centrality det. & analysis det. (see table below)
    bool RemovePileup{true}; // To remove pileup manually using my own criteria approved by Centrality PWG 
    int nBTofmatch[2] = {0};
//cout << "Accepted event: ";
//cout << "RunID: " << runID << " | ";
//cout << "Centrality index (Bin16): " << centrality << " | Centrality index (Bin9): " << centrality9 << endl; 
//cout << "ZDCx: " << ZDCx << " | BBCx: " << BBCx << endl;
//------------------------------------------------ACCEPTANCE TABLE----------------------------------------------------------
//myTracks index |  Centrality/Analysis   | eta_index(SelfCorr) | eta_index(InsideOut)    | eta_index(HalfHalf) | eta_index(Complete) 
//0              |  Analysis              | |eta| < 0.5         | |eta| < 0.5             | 0 < eta < 1.0       | |eta| < 0.5
//1              |  Centrality            | |eta| < 0.5         |-1<eta<-0.5 & 0.5<eta<1  | -1.0 < eta < 0      | |eta| < 1.0
//
//etaIndex | Configuration
//0        | SelfCorr
//1        | InsideOut 
//2        | Complete
//3        | HalfHalf
//--------------------------------------------------------------------------------------------------------------------------

//-------------------------------------------------------

  double meanpt{0.0}; // Mean pt for event-by-event analysis
  const int nRawTracks = myRawTracks.size();  // Shuffle indices to shuffle tracks
  Int_t randomIndices[ nRawTracks ];
  for( Int_t i = 0; i < nRawTracks; i++ ) { randomIndices[ i ] = i; }
  utils->Shuffle( randomIndices, 0, nRawTracks );
  
  for( Int_t i = 0; i < nRawTracks; i++ ) { //1rst Track Loop: Adding charged tracks and nTofMatch hits
    int k = 0; //-->Switching to Analysis
    TrackStruct myTrack;
    RawTrackStruct rawTrack = myRawTracks[ randomIndices[ i ] ];
    myTrack = rawTrack.mTrack;
      eta = myTrack.mEta;
      pt = myTrack.mPt;
    if ( ( eta > -0.5 && eta < 0.5 ) && ( rawTrack.mBTOFindex >= 0 ) && ( rawTrack.mBTOFmatchFlag != 0 ) ){nBTofmatch[ k ]++;} 
         //^<-- nBTOFmatch filling here (Using base centrality criteria to compare w/ RefMult)
    if( ! utils->GoodTrack( params, rawTrack, hm ) ) { continue; }

      k = 0; //-->Switching to Analysis
 // Track selection based on analysis criteria (eta)
  bool AddToAnalysis{false};
  if (etaIndex == 0 || etaIndex == 1 || etaIndex == 2) { if ( eta > -0.5 && eta < 0.5 ) {AddToAnalysis = true;}} // SelfCorr, InsideOut, Complete
  if (etaIndex == 3) {if (eta > 0 && eta < 1.0){AddToAnalysis = true;}} // HalfHalf
  hm->hTrackCuts->Fill(8.5);
     if ( AddToAnalysis ) { // Conditional for adding Ncharged tracks for Analysis (this is QA for analysis + comparing tracks)	
		if (pt > 0.20 && pt <= 10) {
		      hm->hTrackCuts->Fill(9.5);
		      myTracks[k].push_back( myTrack ); //<-- Nch++ here
		      //if ( ( rawTrack.mBTOFindex >= 0 ) && ( rawTrack.mBTOFmatchFlag != 0 ) ){nBTofmatch[ k ]++;} //<-- nBTOFmatch filling here
			 meanpt += pt;
		}
	}
} // <--- End of 1rst Track loop 

//-------------------------------------------------------

//For pileup, true applies to centrality window1, false applies to analysis window2 (always cross-check in StUtilities)
	int i{0}; //--> Analysis 
	//Analysis
	//QA BBCx and ZDCx Plots
        hm->hBBCx_nBTOFmatch[ i ]->Fill( nBTofmatch[ i ] , BBCx );
        hm->hBBCx_FXTMult[ i ]->Fill( Correctedrefmult, BBCx, refMultWeight );
        hm->hZDCx_nBTOFmatch[ i ]->Fill( nBTofmatch[ i ], ZDCx );
        hm->hZDCx_FXTMult[ i ]->Fill( Correctedrefmult, ZDCx,refMultWeight );	
	//Uncomment to include pileup removal:
	//if (!RemovePileup || ( RemovePileup && utils->passnTofMatchRefmultCut((double)myTracks[ i ].size(), (double)nBTofmatch[ i ], false, etaIndex ) ) ){
	if (!RemovePileup || ( RemovePileup && utils->passnTofMatchRefmultCut( Correctedrefmult, (double)nBTofmatch[ i ], false, etaIndex ) ) ){
		hm->hRefMult[ i ]->Fill( refMult ); /// <--Nch filling here (trying w/o refMultWeight + Vz/luminosity corrections)
		hm->hMult[ i ]->Fill( myTracks[ i ].size() ); //<-- charged tracks filling here
		hm->hRefMultCorr[ i ]->Fill( Correctedrefmult, refMultWeight ); // (trying w/ refMultWeight + Vz/luminosity corrections)
                hm->hNch_RefMultCorr->Fill( Correctedrefmult, myTracks[ i ].size(), refMultWeight);
		hm->hFXTMult_nBTOFmatch[ i ]->Fill( nBTofmatch[ i ], Correctedrefmult, refMultWeight ); // <-- Nch vs nBTOFmatch
                hm->hVz[ i ]->Fill( Vz, refMultWeight );
		meanpt = meanpt/myTracks[ 0 ].size();
		hm->hMeanpTNch->Fill( myTracks[ 0 ].size(), meanpt );
		hm->hEventCuts->Fill(9.5);
		hm->hFXTMult_RunID[ 0 ]->Fill( runID, Correctedrefmult, refMultWeight );
   		//if ( Correctedrefmult > 100){ //Avoiding UPC events
		//hm->hZdcEastWest->Fill( ZdcEast, ZdcWest );
                //hm->hRefMult_Zdc->Fill( Correctedrefmult, ZdcEast + ZdcWest, refMultWeight );
                //hm->hZdcEastPWest->Fill( ZdcEast + ZdcWest );
		//}		
		if (!RemovePileup || ( RemovePileup && utils->passnTofMatchRefmultCut( Correctedrefmult, (double)nBTofmatch[ i ], true, etaIndex ) ) ){
		hm->hMult2D->Fill ( Correctedrefmult, refMult2, refMultWeight ); 
		// ^ Nch Correlation filling (Inner [|eta|<0.5] vs Outer [0.5<|eta|<1.0] window) 
		}	
	}

//-------------------------------------------------------


     //2nd Track Loop: More QA plots + pT vs Nch (Unless Centrality config: Nch(refMult) from Centrality, pT from Analysis)
     for( Int_t i = 0; i < nRawTracks; i++ ) { randomIndices[ i ] = i; }
     utils->Shuffle( randomIndices, 0, nRawTracks );
	for( Int_t i = 0; i < nRawTracks; i++ ) {	
	    TrackStruct myTrack;
            RawTrackStruct rawTrack = myRawTracks[ randomIndices[ i ] ];
            if( ! utils->GoodTrack( params, rawTrack, hm ) ) { continue; }

            myTrack = rawTrack.mTrack;
            eta = myTrack.mEta;
	    double phi = myTrack.mPhi, weight{0.0};
	    int bin = hWeighting->FindBin(eta, phi);
	    if ( hWeighting->GetBinContent(bin) > 0 ) weight = 1.0/ ( hWeighting->GetBinContent(bin) ); //For reference: https://root.cern.ch/doc/master/classTH1.html#filling-histograms
  		// Note: trying 1/weight to enfore upweighting
	//cout << "Weight " <<	weight  << endl; // temporary
      	    if (weight > 0) weight = weight * refMultWeight; // This combines real data weighting (Track level weighting for detector effect, inefficiency) and RefMultCorr weighting (Event weighting)
	    //weight = 1;
	    pt = myTrack.mPt;
	    int k{0}; //--> Analysis 
	    bool AddToAnalysis{false};
          if (etaIndex == 0 || etaIndex == 1 || etaIndex == 2) { if ( eta > -0.5 && eta < 0.5 ) {AddToAnalysis = true;}} // SelfCorr, InsideOut, Complete
          if (etaIndex == 3) {if (eta > 0 && eta < 1.0){AddToAnalysis = true;}} // HalfHalf
		if ( AddToAnalysis ) { // Conditional for adding Ncharged tracks for Analysis (this is QA for analysis + comparing tracks)
		if (pt > 0.20 && pt <= 10) {
		if ( !RemovePileup || ( RemovePileup && utils->passnTofMatchRefmultCut( Correctedrefmult, (double)nBTofmatch[ k ], false, etaIndex ) ) ){
			utils->FillAnalyzerTracks(hm, Correctedrefmult, rawTrack, hWeighting, refMultWeight);
			//hm->hPt_Nch->Fill( myTracks[ 0 ].size(), pt, weight );  
			hm->hPt_Nch->Fill( refMult, pt, weight ); 
		        //hm->hpTCent->Fill( centrality, pt, weight);
			//hm->hNchCent->Fill( centrality, Correctedrefmult, weight);
			//hm->hpTCent9->Fill( centrality9, pt, weight);
			//hm->hNchCent9->Fill( centrality9, Correctedrefmult, weight);
			}
		  }
	  }
} //End of 2nd track loop 
  //cout << "\n Event processed with : " << myTracks[0].size() << endl;           
  return kStOK;
}
#endif
