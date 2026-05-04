#ifdef __CINT__

#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

// includes all header files
#include "Parameterization.h"
#include "Bootstrap_Yield.h"
#include "Run_Analysis.h"
#include "Shade_MultHisto.h"
#include "MomentHistogram.h"

// All classes
#pragma link C++ class Parameterization + ;
#pragma link C++ class MomentHistogram + ;

// all functions
#pragma link C++ function runMe;
#endif
