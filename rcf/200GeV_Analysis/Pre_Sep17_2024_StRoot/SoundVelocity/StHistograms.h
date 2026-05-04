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
   
    static const int kNumEtaRanges = 4; // Eta = +/- 0.5, Eta = +/- 0.8, Eta = +/- 1.0, Eta = +/- 1.5
    static const int kNumPtRanges = 4; // pT_low = 0, pT_low = 0.15, pT_low = 0.4, pT_low = 0.8
    TH1D *hEta[kNumEtaRanges*kNumPtRanges];
    TH1D *hTransvMomentum[kNumEtaRanges*kNumPtRanges];
    TH1D *hMult[kNumEtaRanges*kNumPtRanges];
//    TH1D *hChargedTracks[kNumEtaRanges*kNumPtRanges];
	//Track-based histograms
    TH2D *hEtaPhi[kNumEtaRanges*kNumPtRanges];
    TH2D *hPt_Nch[kNumEtaRanges*kNumPtRanges];
    TH2D *hPt_RefMult[kNumEtaRanges*kNumPtRanges];
    TH2D *hPt_dNdEta[kNumEtaRanges*kNumPtRanges];
    TH2D *hChargePt[kNumEtaRanges*kNumPtRanges];
    TH2D *hPtDCA[kNumEtaRanges*kNumPtRanges];
    TH2D *hPtEta[kNumEtaRanges*kNumPtRanges];

    // Event histogram
    TH1D *hRefMult;
    TH3D *hVertex;
    TH1D *hMultRaw;
    TH2D *hMultCent;

    ClassDef(StHistograms,1);

};
#endif
