#include <iostream>

#include "TSystem.h"
#include "TString.h"
#include "TFile.h"

#include "THaShower.h"
#include "THaEvent.h"
#include "MultiFileRun.h"
#include "THaEvData.h"
#include "THaRun.h"
#include "THaAnalyzer.h"
#include "THaVarList.h"

#include "SBSBigBite.h"
#include "SBSBBTotalShower.h"
#include "SBSBBShower.h"

#include "SBSEArm.h"
#include "SBSECal.h"
#include "SBSGEPEArm.h"


void replay_ECal(UInt_t runnum=10491, Long_t nevents=-1, Long_t firstevent=1, const char *fname_prefix="e1209016", UInt_t firstsegment=0, UInt_t maxsegments=10, Int_t maxstream=2, Int_t pedestalmode=0, Int_t cmplots=1, Int_t usesbsgems=1)
		 //int run_number = 124, uint nev = -1, TString start_name = "e1209019", uint nseg = 0)
{

  THaAnalyzer* analyzer = new THaAnalyzer;  

  SBSGEPEArm* earm = new SBSGEPEArm("earm", "GEP electron arm" );
  SBSECal* ecal = new SBSECal("ecal", "ECal");
  ecal->SetDataOutputLevel(1);
  earm->AddDetector(ecal);
  gHaApps->Add(earm);  

  //--- Set up the run we want to replay ---

  TString prefix = gSystem->Getenv("DATA_DIR");

  vector<string> pathlist;
  if( !prefix.IsNull() )
    pathlist.push_back( prefix.Data() );

  if( prefix != "/adaqeb[1-3]/data1" )
    pathlist.push_back( "/adaqeb[1-3]/data1" );

  if( prefix != "/adaq1/data1/sbs" )
    pathlist.push_back( "/adaq1/data1/sbs" );

  // Do not use wildcard filenames with this directory
  if( prefix != "/cache/mss/halla/sbs/raw" )
    pathlist.push_back( "/cache/mss/halla/sbs/raw" );

  if( prefix != "/cache/mss/halla/sbs/GEnRP/raw" )
    pathlist.push_back( "/cache/mss/halla/sbs/GEnRP/raw" );

  for( const auto& path: pathlist ) {
    cout << "search paths = " << path << endl;
  }

  TString ftest(fname_prefix);
  bool test_data = ( ftest == "bbgem" || ftest == "e1209019_trigtest" );

  // Build list of files we'd like to analyze. These are candidates. If some of them
  // do not exist, they will be quietly skipped by the MultiFileRun class. A list of
  // the actually found files will be printed during initialization.
  vector<string> filenames;
  Int_t actual_maxstream = test_data ? 0 : maxstream;
  for( Int_t istream = 0; istream <= actual_maxstream; ++istream ) {
    for( UInt_t iseg = firstsegment; iseg < firstsegment + maxsegments; ++iseg ) {
      TString codafilename;
      if( test_data )
	codafilename.Form("%s_%u.evio.%u", fname_prefix, runnum, istream);
      else
	codafilename.Form("%s_%u.evio.%d.%u", fname_prefix, runnum, istream, iseg);
      cout << "codafilename = " << codafilename << endl;
      filenames.emplace_back(codafilename.Data());
    }
  }
  if( filenames.empty() ) {
    cout << "No files. Check your parameters." << endl;
    return;
  }

  auto* run = new Podd::MultiFileRun(pathlist, filenames);
  run->SetFirstSegment(firstsegment);  // Starting segment number
  run->SetMaxSegments(maxsegments);    // Number of segments
  run->SetMaxStreams(maxstream+1);  // Number of streams

  // Count physics events sequentially, i.e. always starting at 1, for each invocation of
  // analyzer->Process(). This count will be compared to the run's event range, so it
  // can be specified independently of absolute event numbers. This is relevant if
  // the analysis starts with a segment > 0. Actual (raw, CODA) physics event numbers
  // will then start with some largish number (the first event number in the continuation
  // segment), but the analyzer's count will still start at 1 for that analysis,
  // so you can give the run an event range of, say, 1-50000 to analyze 50k events,
  // regardless of what segments are being replayed. The event header written for each
  // entry in the output ROOT tree does contain the actual CODA event number, which is
  // what you usually want for analysis.

  analyzer->SetCountMode(THaAnalyzer::kCountPhysics);

  if( nevents > 0 ) run->SetLastEvent(firstevent+nevents-1);
  run->SetFirstEvent( firstevent );
  run->SetDataRequired(THaRunBase::kDate|THaRunBase::kRunNumber);
  Int_t st = run->Init();
  if( st != THaRunBase::READ_OK ) {
    cerr << "========= Error initializing run" << endl;
    return;
  }

  cout << "n segments to analyze = " << run->GetNFiles() << endl;

  prefix = gSystem->Getenv("OUT_DIR");

  TString outfilename;
  TString outfile_prefix = "ecal";
  UInt_t lastsegment = run->GetLastSegment();
  outfilename.Form( "%s/%s_%d_%ld.root", prefix.Data(), outfile_prefix.Data(), runnum, nevents);
  // if( nevents > 0 ){
  //   outfilename.Form( "%s/%s_replayed_%u_stream%d_%d_seg%u_%u_firstevent%ld_nevent%ld.root", prefix.Data(), fname_prefix, runnum,
  // 		      0, maxstream, firstsegment, lastsegment, firstevent, nevents );
  // } else {
  //outfilename.Form( "%s/%s_fullreplay_%u_stream%d_%d_seg%u_%u.root", prefix.Data(), outfile_prefix, runnum,
  //		    0, maxstream, firstsegment, lastsegment );
  // }


  analyzer->EnableHelicity();

  // analyzer->SetVerbosity(2);
  analyzer->SetMarkInterval(1000);

  analyzer->EnableBenchmarks();

  // Define the analysis parameters
  analyzer->SetOutFile( outfilename.Data() );

  // File to record cuts accounting information
  // prefix = gSystem->Getenv("LOG_DIR");
  // analyzer->SetSummaryFile(Form("%s/%s_%d_stream%d_%d_seg%d_%d.log", prefix.Data(), fname_prefix, runnum,
  // 				0, maxstream, firstsegment, lastsegment));

  prefix = gSystem->Getenv("SBS_REPLAY");
  prefix += "/replay/";

  TString odef_filename = "replay_ECal.odef";

  //no need to define SBS GEM output if we aren't analyzing the data
  if( usesbsgems == 0 ) odef_filename = "replay_gen_noSBSGEMs.odef";

  odef_filename.Prepend( prefix );

  analyzer->SetOdefFile( odef_filename );

  //added cut list in order to have
  //on the farm we only write out events with a track in BigBite
  //TString cdef_filename = "replay_gen_farm.cdef";
  TString cdef_filename = "replay_ECal.cdef";

  //the above cdef file includes almost no significant cuts (basically just any hit in BigBite shower or preshower
  cdef_filename.Prepend( prefix );

  analyzer->SetCutFile( cdef_filename );
  //analyzer->SetCompressionLevel(0); // turn off compression

  st = analyzer->Init(run);
  if( st != 0 ) {
    cerr << "========= Error initializing analyzer" << endl;
    return;
  }
  analyzer->Process(run);     // start the actual analysis


  
  // /////////////////////////////////////////////
  // TString exp = "ecal";
  // // Create file name patterns.
  // string firstname = "ecal_%d";
  // string endname = Form(".evio.%d",nseg);
  // //string endname = Form(".evio");
  // string combined(string(firstname)+endname);
  // const char* RunFileNamePattern = combined.c_str();
  // vector<TString> pathList;
  // pathList.push_back("/cache/halla/sbs/raw");
  // pathList.push_back("/cache/halla/sbs/GEnRP/raw");
  // pathList.push_back(Form("%s/data","/adaqfs/home/a-onl/sbs"));
  // ////////////////////////////////////////////
  // //Kip: Think this needs to be changed to search for the proper format of an ecal cosmic run?
  
  
  // THaAnalyzer* analyzer = new THaAnalyzer;
  // THaEvent* event = new THaEvent;
  // THaRun* run = 0;
  // int seg = 0;
  // bool seg_ok = true;
  // while(seg_ok) {
  //   TString data_fname;
  //   data_fname = TString::Format("%s/%s_%d.evio.0.%d",getenv("DATA_DIR"),start_name.Data(),run_number,seg);
  //   //new THaRun( pathList, Form(RunFileNamePattern, run_number) );
  //   //  for (UInt_t nn=0; nn < pathList.size(); nn++) {
  //   // data_fname = TString::Format("%s/%s_%d.evio.%d",pathList[nn].Data(),start_name.Data(),run_number,seg);
  //   std::cout << "Looking for segment " << seg << " file " << data_fname.Data() << std::endl;
  //   if( gSystem->AccessPathName(data_fname)) {
  //     seg_ok = false;
  //     std::cout << "Segment " << seg << " not found. Exiting" << std::endl;
  //     continue;
  //   }
  //   //  }
  //   run = new THaRun(data_fname);
  //   run->SetLastEvent(nev);

  //   run->SetDataRequired(0);//for the time being
  //   run->SetDate(TDatime());
  
  

  //   analyzer->SetEvent( event );
  //   TString out_dir = gSystem->Getenv("OUT_DIR");
  //   if( out_dir.IsNull() )  out_dir = ".";
  //   TString out_file = out_dir + "/" + exp + Form("_%d_%d.root", run_number,nev);

  //   analyzer->SetOutFile( out_file );
  
  //   analyzer->SetCutFile( "/adaqfs/home/a-onl/sbs/ECAL_replay/git-repo/sbs_devel/SBS-replay/replay/replay_ECal.cdef" );
  //   analyzer->SetOdefFile( "/adaqfs/home/a-onl/sbs/ECAL_replay/git-repo/sbs_devel/SBS-replay/replay/replay_ECal.odef" );

  //   analyzer->SetVerbosity(2);  // write cut summary to stdout
  //   analyzer->EnableBenchmarks();

  //   //--- Analyze the run ---
  //   // Here, one could replay more than one run with
  //   // a succession of Process calls. The results would all go into the
  //   // same ROOT output file

  //   run->Print();
  

  //   analyzer->Process(run);
  //   // Cleanup this run segment
  //   delete run;
    
  //   seg++; // Increment for next search
  // }

  // // Clean up

  // analyzer->Close();
  // delete analyzer;
  // //gHaCuts->Clear();
  // gHaVars->Clear();
  // gHaPhysics->Delete();
  // gHaApps->Delete();

  // // Open the ROOT file so that a user doing interactive analysis can 
  // // immediately look at the results. Of course, don't do this in batch jobs!
  // // To close the file later, simply type "delete rootfile" or just exit

  // //TFile* rootfile = new TFile( out_file, "READ" );
}
