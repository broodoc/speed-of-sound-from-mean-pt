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
    {"nHitsFitMin",{}},
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
  ifstream input( "/star/u/cbroodo/200GeV_Analysis/RCF/200GeV_Analysis/Coarse_Centrality/cutTable.txt" ); /// Check this repo address and make sure it's referenced here correctly!!
    
  if( ! input ) {
    cout<<"Cut Table not found in StParameters::LoadParamTable"<<endl;
    return false;
  }

  int columnSetting = 0;
  int lineNumber = 0;

//  cout<<"----------Reading parameter file----------"<<endl;  
//  cout<<"----------Setting "<<mSetting<<"----------"<<endl;

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
            //cout<<cutName<<"  ";
            //for(uint i = 0; i < paramMap[cutName].size(); i++) { 
              //cout<<paramMap[cutName][i]; 
              //if(paramMap[cutName].size() > 1) {
			//cout<<",";
		//	} 
            //}
            //cout<<endl;
          }
        }

        column++;
      }
    } // Lines with cut values

    lineNumber++;
  } // main while loop over each line

  //cout<<"----------Done reading parameter file----------"<<endl;  

  
  return true;
}

void StParameters::SetBadRuns() {
// Dongsheng Li's bad run list for Run21 3 GeV Analysis:
// https://drupal.star.bnl.gov/STAR/system/files/QA_board_Run21_3GeV_dsli_240809.pdf
   	vector<int> BadRuns_3p0_FXT = {22121045,22122003,22122020,22125011,22158032,22158033,22158036,22159001,22159008,22159009,22159011,22159013,22160002,22161015,22162009,22162023,22163002,22163003,22164033,22163004,22163005,22166024,22166029,22167003,22169008,22171022,22174003,22174004,22174005,22174006,22174007,22174010,22174011,22174012,22174013,22176006,22178013,22171028};
	mBadRunList.insert( mBadRunList.end(), BadRuns_3p0_FXT.begin(), BadRuns_3p0_FXT.end() );
}

bool StParameters::isBadRun( int RunID ) {
//this->SetBadRuns();
for (int i = 0; i < static_cast<int>(mBadRunList.size()); i++){
    if (RunID == mBadRunList[ i ]) {return true;}
} 

return false;
}
#endif
