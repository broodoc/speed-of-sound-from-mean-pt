#include <iostream>
#include "TSystem.h"
#include "TH1F.h"
#include "TTree.h"
#include "TClonesArray.h"

// settings
// "_3p0_FXT"  "_3p2_FXT"  "_3p5_FXT"  "_3p9_FXT"  "_4p5_FXT"  "_7p7_FXT"  "_7p7_COL"  ... "_200p0_COL"

/// NOTE: Don't worry about bool twoPionOnly = true. That section of code is commented out
void runAnalysis( bool twoPionOnly = true, 
		  const char* setting = "_200p0_COL", 
                  const char* fileList = "./picoLists/picoList_physics.list",
                  const char* outfile = "output.root" ){
// /star/u/cbroodo/Coarse_Centrality_Shortcut/picoLists/picoList_physics.list
// "filesPico_200p0_2011.list"
//echo "'`pwd`/input/WeightingHistogram.root'"
//"'$(pwd)'/input/WeightingHistogram.root" 
 //cout << "weightfile : " <<  weightfile << endl; 
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  // Dynamically link needed shared lib
  gSystem->Load("libStTableUtilities");
  gSystem->Load("libTable");
  gSystem->Load("libPhysics");
  gSystem->Load("St_base");
  gSystem->Load("StChain");
  gSystem->Load("St_Tables");
  gSystem->Load("StUtilities"); // StUtilities and other dependencies from /StRoot linked here?
  gSystem->Load("StTreeMaker");
  gSystem->Load("StIOMaker");
  gSystem->Load("StarClassLibrary");
  gSystem->Load("StTriggerDataMaker"); 
  gSystem->Load("StBichsel");
  gSystem->Load("StEvent");
  gSystem->Load("StEventUtilities");
  gSystem->Load("StDbLib");
  gSystem->Load("StEmcUtil");
  gSystem->Load("StPreEclMaker");
  gSystem->Load("StStrangeMuDstMaker");
  gSystem->Load("StRefMultCorr"); ///  <--- StRefMultCorr here
 
  gSystem->Load("StMuDSTMaker"); ///   <--- StMuDSTMaker here
  
  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker"); /// <--- StPicoDstMaker here
  
  cout << '\n' << "About to load in SoundVelocity" << endl;
  gSystem->Load("SoundVelocity"); /// <--- NEED TO CREATE A FILE NAMED SoundVelocity/ and put StRoot/ Inside it

  // create the chain
  StChain *chain = new StChain; 

  StPicoDstMaker* picoMaker;
  StMuDstMaker*  muDstMaker;

  // get first line of input list to checker whether PicoDst or MuDst
  ifstream input( fileList );
  string line;
  getline(input, line);
  input.close();
  TString line_root( line.data() ); 

  int nevents = 0;
  bool pico = true;
  cout<<"FileList name: "<<fileList<<endl;
  cout<<"First line of fileList: "<<line_root<<endl;
  if( line_root.Contains("pico",TString::ECaseCompare::kIgnoreCase) ) {
    cout<<"------------PicoDst mode-------------"<<endl;
    picoMaker = new StPicoDstMaker( StPicoDstMaker::IoRead, fileList, "picoDst");
    pico = true;  
  }
  else {
    cout<<"------------MuDst mode-------------"<<endl;
    muDstMaker = new StMuDstMaker( 0, 0, "", fileList, "MuDst", 100 );
    pico = false;
  }
  //cout << '\n' << "About to create StSoundVelocity object" << endl;
  // SoundVelocity Maker
  StSoundVelocity *anaMaker = new StSoundVelocity( setting, pico ); // <- // StMultiPionFemto referenced here
  // Input root files

  anaMaker->SetDebug( 0 ); //Debug setting

  chain->Init();

  if( pico ) { nevents = picoMaker->chain()->GetEntries(); }
  else { nevents = muDstMaker->chain()->GetEntries(); }
  
  cout<<"total entries: "<<nevents<<endl;

  int istat=0,iev=1;
  EventLoop: if (iev<=nevents && istat!=2) {
    if(iev%10000000==0) cout << "Working on eventNumber " << iev << endl;
    chain->Clear();
    istat = chain->Make(iev); // This should call the Make() method in ALL makers
    if (istat == 2) { cout << "Last  Event Processed. Status = " << istat << endl; }
    if (istat == 3) { cout << "Error Event Processed. Status = " << istat << endl; }
    iev++; goto EventLoop;
          
  } // Event Loop
        
  chain->Finish();

  TFile *fout = new TFile(outfile,"RECREATE");
    fout->mkdir("SoundVelocity");
  
  TIter next( anaMaker->GetHistList() );
  TObject *obj;
  while( obj = next() ) {
     TString name(obj->GetName());
     fout->cd();
     obj->Write();
  }
  
  fout->Close();
  delete fout;
}
