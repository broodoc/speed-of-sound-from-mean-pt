#ifndef STPARAMERS_HH
#define STPARAMERS_HH

#include "TString.h"

#include "StMaker.h"

namespace ParamConsts {

  // One iTPC sector in 2018, fully installed in 2019
  // First run of 2018 isobar data is ~19071038.  First of 2019 is ~20000000
  const int FirstRunWithiTPC = 20000000;

  const int NTPCPadRowsforCrossings = 45;

  const double FmToGeV = 0.197327;

  const int MaxMult = 5000;
  const int MixingBufferSize = 3;
  const double ZvtxWidthMix = 1; // cm
  const int ExDimBins = 1;
  const int MultBins = 4;
  const int KtBins = 3;

  const double Q2LowerBound = 0.005; // GeV
  const double Q2UpperBound = 0.3;
  const double Q2BinSize = 0.005;

  const double KtLowerBound = 0.15;
  const double KtUpperBound = 0.65;
  const double KtBinSize = 0.05;

  const double Q3UpperBound = 0.3;

  const int NKt3bins = 2;
  const double Kt3limits[NKt3bins + 1] = {0, 0.3, 1.0};
  
  // For track splitting cuts 
  // https://www.star.bnl.gov/webdata/dox/html/StTrackTopologyMap_8cxx_source.html
  // Inner sectors are bits 8-20 of first word
  // Outer sectors are bits 21-31 of first word and 0-20 of second word
  const unsigned long tpcMaskWord1 = 0xFFFFFF00; // selects bits 8-31 of first word
  const unsigned long tpcMaskWord2 = 0x1FFFFF; // selects bits 0-20 of second word
  const unsigned long long tpcMaskWord3 = 0x1FFFFFFFFFE; // selects bins 0-40 of third word (iTPC)

  // For track merging cuts
  // hit separation in mm for "r*phi" and z direction, for inner (1-13) and outer (14-45) padrows
  const double MaxDuInner = 1.2;
  const double MaxDuOuter = 2.1;
  const double MaxDzInner = 4.5;
  const double MaxDzOuter = 4.8;

  // For track merging cuts
  // The radius in cm of each TPC pad row
  const double TPCrowRadius[ NTPCPadRowsforCrossings ] = {
    60,        // row1
    64.8,      // row2
    69.6,      // row3
    74.4,      // row4
    79.2,      // row5
    84,        // row6
    88.8,      // row7
    93.6,      // row8
    98.8,      // row9
    104,       // row10
    109.2,     // row11
    114.4,     // row12
    119.6,     // row13
    127.195,   // row14
    129.195,   // row15
    131.195,   // row16
    133.195,   // row17
    135.195,   // row18
    137.195,   // row19
    139.195,   // row20
    141.195,   // row21
    143.195,   // row22
    145.195,   // row23
    147.195,   // row24
    149.195,   // row25
    151.195,   // row26
    153.195,   // row27
    155.195,   // row28
    157.195,   // row29
    159.195,   // row30
    161.195,   // row31
    163.195,   // row32
    165.195,   // row33
    167.195,   // row34
    169.195,   // row35
    171.195,   // row36
    173.195,   // row37
    175.195,   // row38
    177.195,   // row39
    179.195,   // row40
    181.195,   // row41
    183.195,   // row42
    185.195,   // row43
    187.195,   // row44
    189.195};  // row45

}

namespace constants {
  const double mass_pion = 0.13957;
};

class StParameters : public TObject {
  public:
    StParameters( TString setting );
    ~StParameters();

    bool LoadParamTable();
    void SetBadRuns();
    bool isBadRun( int RunID );

    // potentially several values for one paramater, ie triggers
    std::map<string, vector<double>> paramMap;
    TString mSetting;
    vector<int> mBadRunList;

    ClassDef(StParameters,1);
};
#endif
