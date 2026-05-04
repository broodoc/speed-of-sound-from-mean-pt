#ifndef STEVENTCOLLECTION_CXX
#define STEVENTCOLLECTION_CXX

#include "TH1F.h"
#include "TH2F.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"

#include "StEventCollection.h"

ClassImp(StEventCollection)

//------------------------------------------------------------------------
StEventCollection::StEventCollection()
{

  for( int i = 0; i < ParamConsts::MixingBufferSize; i++ ) {
    mEvtStr[ i ].tracks.clear();
  }
}

//------------------------------------------------------------------------
StEventCollection::~StEventCollection()
{
}

//------------------------------------------------------------------------
void StEventCollection::FIFOShift() 
{
  for( UShort_t i = ParamConsts::MixingBufferSize - 1; i > 0; i-- ) {
    
    mEvtStr[ i ].tracks.clear(); // Clear the vector before it's refilled
    
    for( UInt_t j = 0; j < mEvtStr[i-1].tracks.size(); j++) {
      mEvtStr[ i ].tracks.push_back( mEvtStr[ i-1 ].tracks[ j ] );
    }
  }
  
  mEvtStr[ 0 ].tracks.clear(); // Clear vector for "current" event
}

#endif
