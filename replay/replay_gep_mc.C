#include <iostream>
#include <unordered_map> 

#include "TSystem.h"
#include "THaGlobals.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TObject.h"

#include "THaEvData.h"
#include "THaRun.h"
#include "THaAnalyzer.h"
#include "THaVarList.h"
#include "THaInterface.h"

#include "SBSBigBite.h"
#include "SBSBBShower.h"
#include "SBSBBTotalShower.h"
#include "SBSGRINCH.h"
#include "SBSEArm.h"
#include "SBSHCal.h"
#include "SBSGEMSpectrometerTracker.h"
#include "SBSTimingHodoscope.h"

#include "SBSSimDecoder.h"


TDatime get_datime(uint gepconfig)
/* Returns TDatime for a given GEn configuration */
{
  std::unordered_map<uint,TDatime> m = {{1, "2024-10-01 00:00:00"},
					{2, "2024-10-01 00:00:00"},
					{3, "2024-10-01 00:00:00"}};
  if (m.find(gepconfig)==m.end()) 
    throw std::invalid_argument("Invalid SBS config!! Valid options are: 1,2,3");
  return m[gepconfig];
}

void replay_gep_mc(const char* filebase, uint gepconfig, uint nev = -1, TString experiment="gep")
{
  SBSEArm* earm = new SBSEArm("earm", "GEP electron arm" );
  //earm->AddDetector( new SBSCalorimeter("ecal", "ECal") );
  //earm->AddDetector( new SBSCDet("cdet", "coordinate detector") );
  //gHaApps->Add(earm);
  SBSEArm *harm = new SBSEArm("sbs","Hadron Arm with HCal");
  harm->AddDetector( new SBSHCal("hcal","HCAL") );
  harm->AddDetector( new SBSGEMSpectrometerTracker("gemFT", "Front tracker") );
  harm->AddDetector( new SBSGEMSpectrometerTracker("gemFPP", "Focal Plane Polarimeter") );
  gHaApps->Add(harm);

  //bigbite->SetDebug(2);
  //harm->SetDebug(2);

  THaAnalyzer* analyzer = new THaAnalyzer;
  
  THaInterface::SetDecoder( SBSSimDecoder::Class() );
  
  TString run_file = Form("%s.root", filebase);
  if(std::getenv("DATA_DIR")){
    run_file = Form("%s/%s.root", std::string(std::getenv("DATA_DIR")).c_str(), filebase);
  }
  
  if( gSystem->AccessPathName(run_file) ) {
    Error("replay.C", "Input file does not exist: %s", run_file.Data() );
    exit(-1);
  }
  
  THaRunBase *run = new SBSSimFile(run_file.Data(), "gmn", "");
  run->SetFirstEvent(0);

  cout << "Number of events to replay (-1=all)? ";
  //if( nev > 0 )
  //run->SetFirstEvent(110);
  run->SetLastEvent(nev);
  
  run->SetDataRequired(0);
  run->SetDate(get_datime(gepconfig));
  //run->SetDate(TDatime());
  
  TString out_dir = gSystem->Getenv("OUT_DIR");
  if( out_dir.IsNull() )
    out_dir = ".";
  TString out_file = out_dir + "/" + Form("replayed_%s.root", filebase);
  
  analyzer->SetOutFile( out_file.Data() );
  cout << "output file " << out_file.Data() << " set up " << endl; 
  // File to record cuts accounting information
  analyzer->SetSummaryFile("sbs_hcal_test.log"); // optional

  // Change the cratemap to point to the sim one
  analyzer->SetCrateMapFileName("sbssim_cratemap");

  cout << "sim crate map setup " << endl;

  TString prefix = gSystem->Getenv("SBS_REPLAY");
  prefix += "/replay/";
  
  TString odef_filename = "replay_gep_mc.odef";
  odef_filename.Prepend( prefix );
  analyzer->SetOdefFile( odef_filename );
  
  //added cut list in order to have 
  TString cdef_filename = "replay_gep_mc.cdef";
  cdef_filename.Prepend( prefix );
  analyzer->SetCutFile( cdef_filename );

  //analyzer->SetCutFile( "replay_gep.cdef" );
  //analyzer->SetOdefFile( "replay_gep_mc.odef" );
  
  cout << "cut file and out file processed " << endl;
  
  analyzer->SetVerbosity(2);  // write cut summary to stdout
  analyzer->EnableBenchmarks();
  
  run->Print();
  
  cout << "about to process " << endl;
  analyzer->Process(run);

  // Clean up

  analyzer->Close();
  delete analyzer;
  //gHaCuts->Clear();
  gHaVars->Clear();
  gHaPhysics->Delete();
  gHaApps->Delete();
 
  
}

int main(int argc, char *argv[])
{
  new THaInterface( "The Hall A analyzer", &argc, argv, 0, 0, 1 );
  //const 
  string filebase; 
  uint nev = -1;
  if(argc<2 || argc>3){
    cout << "Usage: replay_gen filebase(char*) nev(uint)" << endl;
    return -1;
  }
  filebase = argv[1];
  if(argc==3) nev = atoi(argv[2]);

  replay_gep_mc(filebase.c_str(), nev);
  return 0;
}
