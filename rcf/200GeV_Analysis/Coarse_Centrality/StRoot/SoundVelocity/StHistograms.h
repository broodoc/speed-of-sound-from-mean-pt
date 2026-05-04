#ifndef STHISTOGRAMS_HH
#define STHISTOGRAMS_HH

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TObject.h"

#include "StMaker.h"

#include "StParameters.h"

class StHistograms : public TObject {

  public:
    StHistograms( StParameters *params );
    ~StHistograms();

    void InitHistograms( StParameters *params );
   
    static const int kNumEtaRanges = 2;
    static const int kNumPtRanges = 1; ///
    /// pT_low = 0.20,
    TH1D *hEta[kNumEtaRanges*kNumPtRanges];
    TH1D *hTransvMomentum[kNumEtaRanges*kNumPtRanges];
    TH1D *hMult[kNumEtaRanges*kNumPtRanges];
    TH1D *hRefMultCorr[ kNumEtaRanges*kNumPtRanges ];
    TH1D *hRefMult[ kNumEtaRanges*kNumPtRanges ];
    TH1D *hVz[kNumEtaRanges*kNumPtRanges];
//    TH1D *hChargedTracks[kNumEtaRanges*kNumPtRanges];

    TH2D *hEtaPhi[kNumEtaRanges*kNumPtRanges];
    TH2D *hPhi_Nch[ kNumEtaRanges*kNumPtRanges ];
    TH2D *hPt_RefMultCorr;
    TH2D *hPt_Nch;
    TH2D *hPt_dNdEta[kNumEtaRanges*kNumPtRanges];
    TH2D *hChargePt[kNumEtaRanges*kNumPtRanges];
    TH2D *hPtDCA[kNumEtaRanges*kNumPtRanges];
    TH2D *hPtEta[kNumEtaRanges*kNumPtRanges];
    TH2D *hEta_Nch[ kNumEtaRanges*kNumPtRanges ];    
    TH2D *hNch_NhitsFit[ kNumEtaRanges*kNumPtRanges ];
    TH2D *hEta_Nhits[ kNumEtaRanges*kNumPtRanges ];
    TH2D *hNch_NhitsdEdx[ kNumEtaRanges*kNumPtRanges ];
    TH2D *hNch_DCA[ kNumEtaRanges*kNumPtRanges ];
    TH2D *hMeanpTNch;
    TH2D *hMult2D;
    TH1D *hEventCuts;
    TH1D *hTrackCuts;     
    TH2D *hpTCent;
    TH2D *hNchCent;
    TH2D *hNchCent9;
    TH2D *hpTCent9;
 
    TH2D *hFXTMult_nBTOFmatch[kNumEtaRanges*kNumPtRanges];
    TH2D *hFXTMult_RunID[kNumEtaRanges*kNumPtRanges];
    TH2D *hnBTOFmatch_RunID[kNumEtaRanges*kNumPtRanges];
    TH2D *hFXTMult_nBTOFmatch_NoCuts;

    TH2D *hBBCx_RunID[kNumEtaRanges*kNumPtRanges];
    TH2D *hBBCx_nBTOFmatch[kNumEtaRanges*kNumPtRanges];
    TH2D *hBBCx_FXTMult[kNumEtaRanges*kNumPtRanges];
    TH2D *hZDCx_FXTMult[ kNumEtaRanges*kNumPtRanges ];
    TH2D *hZDCx_nBTOFmatch[ kNumEtaRanges*kNumPtRanges ];

    TH3D *hPt_Nch_eta;
// Event histograms
    //TH1D *hRefMult;
    TH3D *hVertex;
    TH1D *hMultRaw;
    TH2D *hMultCent;
    //TH1D *hVz;
    TH2D *hPt_FXTMult;

    TH2D *hNch_RefMultCorr;
    TH2D *hZdcEastWest;
    TH2D *hRefMult_Zdc;
    TH1D *hZdcEastPWest;

    ClassDef(StHistograms,1);

};
#endif
