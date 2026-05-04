#ifndef STEVENTCOLLECTION_HH
#define STEVENTCOLLECTION_HH

#include "StMaker.h"

#include "StUtilities.h"
#include "StParameters.h"

class StEventCollection : public TObject {
 public:
  StEventCollection();
  ~StEventCollection();

  void FIFOShift();

  EventStruct mEvtStr[ ParamConsts::MixingBufferSize ];

  ClassDef(StEventCollection,1);
};
#endif
