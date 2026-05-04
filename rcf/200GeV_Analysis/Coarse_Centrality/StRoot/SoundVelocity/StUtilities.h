#ifndef STUTILITIES_HH
#define STUTILITIES_HH

#include "TVector3.h"
#include "TRandom3.h"
#include "TH2D.h"

#include "StMaker.h"
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

#include "StParameters.h"
#include "StHistograms.h"

struct TrackStruct {
  Int_t mIdx;
  Float_t mCharge;
  UInt_t mFirstWordTopologyMap;
  UInt_t mSecondWordTopologyMap;
  ULong64_t mThirdWordTopologyMap;
  TVector3 mP;
  Float_t mPt;
  Float_t mPhi;
  Float_t mEta;
  Float_t mMass;
  Float_t mE;
  Short_t mFlag;
  Int_t mBTOFindex;
  Int_t mNhits;
  Int_t mBTOFmatchFlag;
  TVector3 mCrossings[ ParamConsts::NTPCPadRowsforCrossings ];
  Int_t mID_MC;
  TVector3 mP_MC;
};

struct RawTrackStruct {
  TrackStruct mTrack;
  // other params which are not of interest after cuts
  Int_t mRawTrackID;
  Float_t mDCA;
  Int_t mNhits;
  Int_t mNhitsFit;
  Int_t mNhitsMax;
  Int_t mNhitsDedx;
  UInt_t mHFThitsMap;
  Float_t mNsigmaPion;
  Float_t mNsigmaElectron;
  Float_t mNsigmaKaon;
  Float_t mNsigmaProton;
  Int_t mBTOFindex;
  Int_t mBTOFmatchFlag;
  Short_t mFlag;
  Float_t mNsigmaPionTOF;
  Float_t mBTOFylocal;
  Float_t mBTOFzlocal;
  Float_t mM2;
  Float_t mBeta;
  Float_t mBetaPion;
};


struct PairStruct {
  Double_t qInv;
  Double_t kT;
  Double_t qOut;
  Double_t qSide;
  Double_t qLong;
  TVector3 mom1;
  TVector3 mom2;
  Double_t idx1;
  Double_t idx2;
  Double_t charge1;
  Double_t charge2;
  Double_t SL;
  Double_t FMH;
  TVector3 mom1_MC;
  TVector3 mom2_MC;
  Double_t qInv_MC;
};

struct EventStruct {
  vector<TrackStruct> tracks;
  //vector<MCpartStruct> mcParticles;
};

enum SqrtS {
  _3p0_FXT,
  _3p2_FXT,
  _3p5_FXT,
  _3p9_FXT,
  _4p5_FXT,
  _7p7_FXT,
  _7p7_COL,
  _9p0_COL,
  _11p5_COL,
  _14p5_COL,
  _19p6_COL,
  _27p0_COL,
  _39p0_COL,
  _54p0_COL,
  _62p0_COL,
  _200p0_COL,
  _NumEnergies
};

class StUtilities : public TObject {
  public:
//    StUtilities( TString setting, TString Kfilename, TString MRCfilename, TString Weightfilename );
    StUtilities( TString setting );
    ~StUtilities();

    bool GoodPicoEvent( StParameters *params, StPicoEvent *event, StHistograms *hm );
    bool GoodMuEvent( StParameters *params, StMuEvent *event );

    void FillPicoRawTracks( StPicoDst *picoDst, StParameters *params, vector<RawTrackStruct> *tracks, StHistograms *hm );
    void FillMuRawTracks( StMuDst *muDst, StParameters *params, vector<RawTrackStruct> *tracks );

    void FillAnalyzerTracks(StHistograms *hm, double refMult, RawTrackStruct rawTrack, TH2D* hWeighting, double refMultWeight);

    bool GoodTrack( StParameters *params, RawTrackStruct rawTrack, StHistograms *hm );
    int GetMultBin( int centBin1to16 );
    int GetVzBin( StParameters *params, double vz );
    int CalcRefMult( StPicoDst *pDst );
    int CalcRefMult2( StPicoDst *pDst );
    int CalcRefMult( StMuDst *muDst );
    bool BadRun( StParameters *params, int run );
    void Shuffle( Int_t *iarr, Int_t i1, Int_t i2 );
    void SetiTPCCase( int run );
    void LoadFSIfile( TString filename );
    void PrintTrackInfo( RawTrackStruct track );    
    bool passnTofMatchRefmultCut(Double_t refmult, Double_t ntofmatch, bool window1, int etaIndex);
    Double_t calcPileUpRefMult(Double_t ntofmatch, Double_t x0, Double_t x1, Double_t x2, Double_t x3, Double_t x4);
    bool isInPileUpRefMultLimits( Double_t refmult, Double_t refmultcutmin, Double_t refmultcutmax );

    // Constant Bit Masks for HFT layers
    const static UInt_t PXL1_MASK = 1 << 0;  // Bit 0 for PXL1
    const static UInt_t PXL2_MASK = 1 << 1;  // Bit 1 for PXL2
    const static UInt_t IST_MASK  = 1 << 2;  // Bit 2 for IST
    const static UInt_t SSD_MASK  = 1 << 3;  // Bit 3 for SSD
 
    TString mSetting;
    bool miTPC_In;
    bool mFXT;
    TH2D *mFSI_ss;
    TH2D *mFSI_os;

    ClassDef(StUtilities,1);

};


#endif
