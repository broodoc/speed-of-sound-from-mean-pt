#ifndef STPARAMERS_HH
#define STPARAMERS_HH

#include "TString.h"
#include "TObject.h"

#include <map>
#include <vector>
#include <string>

class StParameters : public TObject {
  public:
    StParameters();
    ~StParameters();

    bool LoadParamTable();

    std::map<std::string, std::vector<std::string>> paramMap;

    ClassDef(StParameters, 1);
};

#endif
