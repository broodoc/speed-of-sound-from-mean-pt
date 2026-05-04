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
//  vector<vector<vector<int>>> pairFinder;
//  vector<vector<vector<int>>> pairFinderInit;
  StRefMultCorr *refMultCorr;

  ClassDef(StSoundVelocity, 1);
};
#endif
