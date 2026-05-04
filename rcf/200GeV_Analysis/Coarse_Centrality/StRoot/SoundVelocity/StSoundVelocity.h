#ifndef STMULTIPIONFEMTO_HH /// Where do these come from? Need an equivalent for StSoundVelocity.h?
#define STMULTIPIONFEMTO_HH /// Where do these come from? Need an equivalent for StSoundVelocity.h?

#include "TString.h"

#include "StMaker.h"
#include "StRefMultCorr/CentralityMaker.h"
#include "StRefMultCorr/StRefMultCorr.h"

#include "StParameters.h"
#include "StUtilities.h"
#include "StEventCollection.h"
#include "StHistograms.h"

#include <fstream>
#include <iostream>
#include <string>
#include "TFile.h"
#include "TH2D.h"

class StRefMultCorr; //Added for new StRefMultCorr in analysis

class StSoundVelocity : public StMaker {
 public:
  StSoundVelocity(TString setting, bool picoCase );
  ~StSoundVelocity();

  int Init();
  int Make();
  inline void SetDebug( int db ) 	    { mDebug = db; }
//  inline void SetKfile( TString name )      { mKfile = name; }
//  inline void SetMRCfile( TString name )    { mMRCfile = name; }
//  inline void SetWeightfile( TString name ) { mWeightfile = name; }


 protected:
  StUtilities *utils;
  StParameters *params;
  StHistograms *hm;
  EventStruct *evt;
  StEventCollection ***EC;

 private:

  virtual const char *GetCVS() const {
    static const char cvs[]="Tag $Name:  $Id: built " __DATE__ " " __TIME__ ; return cvs;
  }

  int mDebug;
//  TString mKfile;
//  TString mMRCfile;
//  TString mWeightfile;
//  bool mTwoPionOnly;
  TString mSetting;
  bool mPicoCase;
  Int_t mStream;
  //const char* weightFileName;
//  vector<vector<vector<int>>> pairFinder;
//  vector<vector<vector<int>>> pairFinderInit;
  StRefMultCorr *refMultCorr;
  TH2D* hWeighting = new TH2D("hWeighting", "hWeighting", 31, -1.55, 1.55, 100, -3.14159, 3.14159);
//Added for updated StRefMultCorr in Analysis
  void setRefMutCorr(StRefMultCorr* gRefMultCorr) { refMultCorr = gRefMultCorr; }
  StRefMultCorr* getRefMultCorr() { return refMultCorr; }

  ClassDef(StSoundVelocity, 1);
};
#endif
