#ifndef STPARAMERS_CXX
#define STPARAMERS_CXX

#include <cstdlib>
//#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "TString.h"

#include "StParameters.h"

using namespace std;

//------------------------------------------------------------------------
StParameters::StParameters( TString setting ) {

  mSetting = setting;
/// Eta cuts removed (check )
  paramMap = { 
    {"Trigger", {}},
    {"VzMin", {}},
    {"VzMax", {}},
    {"VrCut", {}},
    {"nHitsMin", {}},
    {"nHitsDedxMin", {}},
    {"dcaMax", {}},
    {"VzVpdMin", {}},
    {"VzVpdMax", {}},
  };
}

//------------------------------------------------------------------------
StParameters::~StParameters() {}

//------------------------------------------------------------------------
bool StParameters::LoadParamTable() {
  string line;
  ifstream input( "/star/u/cbroodo/200GeV_Analysis/RCF/Dhevan_Run11_Modified_Analysis/Run11_Mar_6_2024/cutTable.txt" ); /// Check this repo address and make sure it's referenced here correctly!!
    
  if( ! input ) {
    cout<<"Cut Table not found in StParameters::LoadParamTable"<<endl;
    return false;
  }

  int columnSetting = 0;
  int lineNumber = 0;

  cout<<"----------Reading parameter file----------"<<endl;  
  cout<<"----------Setting "<<mSetting<<"----------"<<endl;

  while( std::getline(input, line).good() ) { 
    //cout<<"Inside while loop"<<endl; /// debugging
    std::istringstream iss(line);
    string segment = "";
    int column = 0;

    if( lineNumber == 0 ) { // Header line with all energy setting strings

      while( iss >> segment ) {
        if( segment.compare( mSetting.Data() ) == 0 ) {
          columnSetting = column;
        }
        column++;
      }
      if( columnSetting == 0 ) { 
        cout<<"No matching energy setting found in StParameters::LoadParamTable"<<endl;
        return false;
      }

    } // Header line
    else { // Cut values

      string cutName = "";

      while( iss >> segment ) {
        if( column == 0 ) {
          // verify that param string in txt file matches one in paramMap
          if( paramMap.find( segment ) != paramMap.end() ) { 
            cutName = segment;
          } else { 
            cout<<"No parameter found in StParameters::LoadParamTable"<<endl; 
            return false;
          }
        }
        else {
          if( column == columnSetting ) {
            // parameter could have multiple values, i.e. triggers
            stringstream ss( segment );
            while( ss.good() ) { // loop over multiple params
              string substr;
              getline( ss, substr, ',' );
              paramMap[ cutName ].push_back( std::stod( substr ) );
            }
            cout<<cutName<<"  ";
            for(uint i = 0; i < paramMap[cutName].size(); i++) { 
              cout<<paramMap[cutName][i]; 
              if(paramMap[cutName].size() > 1) {cout<<",";} 
            }
            cout<<endl;
          }
        }

        column++;
      }
    } // Lines with cut values

    lineNumber++;
  } // main while loop over each line

  cout<<"----------Done reading parameter file----------"<<endl;  

  
  return true;
}

void StParameters::SetBadRuns() {}

#endif
