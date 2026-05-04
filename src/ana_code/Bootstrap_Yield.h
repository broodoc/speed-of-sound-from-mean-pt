#ifndef BOOTSTRAP_YIELD_H
#define BOOTSTRAP_YIELD_H

#include "TFile.h"
#include "TStyle.h"
#include "TString.h"
#include <TH1.h>
#include <TH2.h>
#include "TLatex.h"
#include <TCanvas.h>
#include <TStyle.h>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <random>
#include "TRandom.h"
#include <TError.h>
#include <TGraph.h>
#include <TMath.h>
#include <TTree.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include "TRandom3.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

//void Bootstrap_Yield(string dir = "",string prepdir = "",string refMult = "", string pTMult = "", bool debug = false);

void Bootstrap_Yield(string dir = "", string prepdir = "", string refMult = "", string pTMult = "", string InEx = "", string DebugVerbose = "");

#endif
