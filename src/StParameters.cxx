#ifndef STPARAMERS_CXX
#define STPARAMERS_CXX

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "TString.h"

#include "StParameters.h"

using namespace std;

//------------------------------------------------------------------------
StParameters::StParameters() {
  paramMap = {
      {"Directory", {}},
      {"Multname", {}},
      {"MultFromData", {}},
      {"pTMultname", {}},
      {"Glaubername", {}},
      {"COMenergy", {}},
      {"FXT/COL", {}},
      {"System", {}},
      {"PrepDir", {}},
      {"AnaEta1", {}},
      {"AnaEta2", {}},
      {"AnaEta3", {}},
      {"AnaEta4", {}},
      {"CentEta1", {}},
      {"CentEta2", {}},
      {"CentEta3", {}},
      {"CentEta4", {}},
      {"CentpTlow", {}},
      {"Inclusive/Exclusive(I/E)", {}},
      {"DebugVerbose(Y/N)", {}},
      {"Overwrite(Y/N)", {}},
      {"Glauber(Y/N)", {}},
      {"DisplayFitRanges(Y/N)", {}},
      {"GenerateFigures(Y/N)", {}},
      {"MomentAna(Y/N)", {}},
  };
}

//------------------------------------------------------------------------
StParameters::~StParameters() {}

//------------------------------------------------------------------------
bool StParameters::LoadParamTable() {
  std::ifstream input("/Users/calebbroodo/Desktop/sos_analysis/input_table.txt");

  if (!input) {
    std::cout << "Input table not found in StParameters::LoadParamTable" << std::endl;
    return false;
  }

  std::string line;
  std::cout << "---------- Reading input table ----------" << std::endl;

  while (std::getline(input, line)) {
    if (line.empty()) continue;
    std::istringstream iss(line);
    std::string key, value;
    if (!(iss >> key >> value)) continue;

    if (paramMap.find(key) != paramMap.end()) {
      paramMap[key].push_back(value);
      std::cout << key << " = " << value << std::endl;
    } else {
      std::cerr << "Warning: Unknown parameter '" << key << "'" << std::endl;
    }
  }

  return true;
}


#endif
