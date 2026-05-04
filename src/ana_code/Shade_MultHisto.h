#ifndef SHADE_MULTHISTO_H
#define SHADE_MULTHISTO_H

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

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

int customColors[] = {
    kRed, kBlue, kGreen+2, kMagenta, kCyan, kOrange+7, kViolet, kTeal, kPink+10,
    kAzure+10, kSpring+9, kYellow+2, kOrange+5, kViolet+7, kGreen+4, kCyan+3,
    kMagenta+2, kPink+5, kRed+2, kBlue+3, kTeal+2, kAzure+7, kSpring+5, kPink+7,
    kOrange+3, kYellow+3, kGreen+3, kMagenta+4, kCyan+4, kBlue+2, kRed+3,
    kViolet+6, kSpring+3, kTeal+3, kAzure+3
};

void Shade_MultHisto(string dir = "",string prepdir = "",string refMult = "", bool debug = false);

#endif

