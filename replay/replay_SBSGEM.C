#include "TSystem.h"
#include "TList.h"
#include "THaRun.h"
#include "THaEvent.h"
#include "THaAnalyzer.h"
#include "THaApparatus.h"
#include "TString.h"
#include "TClonesArray.h"
#include "SBSGenericDetector.h"
#include "SBSBBTotalShower.h"
#include "SBSBBShower.h"


#include "SBSGEMSpectrometerTracker.h"
#include "SBSBigBite.h"
//#include "SBSGEMStand.h"
//#include "SBSBigBite.h"

void replay_SBSGEM(UInt_t runnum=1342, Long_t nevents=5000, Long_t firstevent=0, const char *fname_prefix="e1209016", UInt_t firstsegment=0, UInt_t maxsegments=1, Int_t pedestalmode=0, Int_t cmplots=1, Int_t parallelmode=0, Int_t singlemode=0, Int_t file_seg=0){

  //  gSystem->Load("libsbs.so");

  int stream = 0;

  if(pedestalmode == 1){
    cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
    cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;
    cout << "-----Running in pedestal mode-----" << endl;
    cout << "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-" << endl;
    cout << "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*" << endl;
  }
  
  SBSBigBite   *sbs = new SBSBigBite("sbs", "Generic apparatus");
  //SBSGEMStand *gems = new SBSGEMStand("gems", "Collection of GEMs in stand");
  SBSGEMSpectrometerTracker *sbsgem = new SBSGEMSpectrometerTracker("gem", "BigBite Hall A GEM data");
    
  sbs->AddDetector(sbsgem);

  //Add trigger TDC info:
  /*  ///////////////////////////////////
  SBSGenericDetector* tdctrig= new SBSGenericDetector("tdctrig","BigBite shower TDC trig");
  tdctrig->SetModeADC(SBSModeADC::kNone);
  tdctrig->SetModeTDC(SBSModeTDC::kTDC);
  tdctrig->SetStoreEmptyElements(kFALSE);
  bb->AddDetector( tdctrig );

  SBSBBTotalShower* ts= new SBSBBTotalShower("ts", "sh", "ps", "BigBite shower");
  ts->SetDataOutputLevel(0);
  bb->AddDetector( ts );
  ts->SetStoreEmptyElements(kFALSE);
  ts->GetShower()->SetStoreEmptyElements(kFALSE);
  ts->GetPreShower()->SetStoreEmptyElements(kFALSE);

  SBSGenericDetector* bbtrig= new SBSGenericDetector("bbtrig","BigBite shower ADC trig");
  bbtrig->SetModeADC(SBSModeADC::kADC);
  bbtrig->SetModeTDC(SBSModeTDC::kTDC);
  bbtrig->SetStoreEmptyElements(kFALSE);
  bb->AddDetector( bbtrig );
  */////////////////////
  //bool pm =  ( pedestalmode != 0 );
  //this will override the database setting:
  ( static_cast<SBSGEMTrackerBase *> (sbsgem) )->SetPedestalMode( pedestalmode );
  ( static_cast<SBSGEMTrackerBase *> (sbsgem) )->SetMakeCommonModePlots( cmplots );  
//
  //  Steering script for Hall A analyzer demo
  //
  
  // Set up the equipment to be analyzed.
  
  // add the two spectrometers with the "standard" configuration
  // (VDC planes, S1, and S2)
  // Collect information about a easily modified random set of channels
  // (see DB_DIR/*/db_D.dat)
  /*
    THaApparatus* DECDAT = new THaDecData("D","Misc. Decoder Data");
    gHaApps->Add( DECDAT );
  */
  

  // Set up the analyzer - we use the standard one,
  // but this could be an experiment-specific one as well.
  // The Analyzer controls the reading of the data, executes
  // tests/cuts, loops over Apparatus's and PhysicsModules,
  // and executes the output routines.
  THaAnalyzer* analyzer = new THaAnalyzer;
  
  gHaApps->Add(sbs);

  // A simple event class to be output to the resulting tree.
  // Creating your own descendant of THaEvent is one way of
  // defining and controlling the output.
  THaEvent* event = new THaEvent;

  TString prefix = gSystem->Getenv("DATA_DIR");
  
  bool segmentexists = true;
  int segment=firstsegment; 

  int lastsegment=firstsegment;
  
  TClonesArray *filelist = new TClonesArray("THaRun",10);

  TDatime now = TDatime();
  
  vector<TString> pathlist;
  pathlist.push_back( prefix );

  if( prefix != "/adaqeb2/data1" )
    pathlist.push_back( "/adaqeb2/data1" );

  if( prefix != "/data/raw" )
    pathlist.push_back( "/data/raw" );

  if( prefix != "/adaq1/data1/sbs" )
    pathlist.push_back( "/adaq1/data1/sbs" );

  if( prefix != "/cache/halla/sbs/raw" )
    pathlist.push_back( "/cache/halla/sbs/raw" );

  for( int i=0; i<pathlist.size(); i++ ){
    cout << "search paths = " << pathlist[i] << endl;
  }

  TDatime RunDate = TDatime(); 

  int max1 = maxsegments;

  int segcounter=0;
  //New: Always add first segment even if not included in request:
  // if( firstsegment > 0 ){
  //   TString codafilename;
  //   codafilename.Form( "%s_%d.evio.%d.%d", fname_prefix, runnum, stream, 0 );
    
  //   TString ftest(fname_prefix);

  //   if( ftest == "bbgem" || ftest == "e1209019_trigtest" ){
  //     codafilename.Form("%s_%d.evio.%d", fname_prefix, runnum, 0 );
  //   }

  //   new( (THaRun*) (*filelist)[segcounter] ) THaRun( pathlist, codafilename.Data(), "GMN run" );

  //   ( (THaRun*) (*filelist)[segcounter] )->SetDataRequired(THaRunBase::kDate|THaRunBase::kRunNumber);
  //   //( (THaRun*) (*filelist)[segcounter] )->Init();
  //   //Not sure if we need to call Init()
  //   ( (THaRun*) (*filelist)[segcounter] )->Init();
  //   RunDate = ( (THaRun*) (*filelist)[segcounter] )->GetDate();

  //   segcounter++;
  //   max1++;
  // }
  

  //This loop adds all file segments found to the list of THaRuns to process:
  while( segcounter < max1 && segment - firstsegment < maxsegments ){
    
    TString codafilename;
    //codafilename.Form( "%s/bbgem_%d.evio.%d", prefix.Data(), runnum, segment );
    // codafilename.Form("%s_%d.evio.%d.%d", fname_prefix, runnum, stream, segment );
    codafilename.Form("%s_SBSGEMs_%d.evio.0.%d", fname_prefix, runnum, segment );
    
    TString ftest(fname_prefix);
    if( ftest == "bbgem" || ftest == "e1209019_trigtest" ){
      codafilename.Form("%s_%d.evio.0.%d", fname_prefix, runnum, segment );
    }

    if( ftest == "eel125"){
      codafilename.Form("%s_%d.evio.%d", fname_prefix, runnum, segment );
    }
    if( ftest == "e1209016"){
      codafilename.Form("%s_SBSGEMs_%d.evio.0.%d", fname_prefix, runnum, segment);
     }
    if( ftest == "e1209016"){
      codafilename.Form("%s_%d.evio.%d.0", fname_prefix, runnum, segment);
      }
    
    
    segmentexists = false;

    cout << "codafilename = " << codafilename << endl;

    for( int ipath=0; ipath<pathlist.size(); ipath++ ){
      TString searchname;
      searchname.Form( "%s/%s", pathlist[ipath].Data(), codafilename.Data() );
      cout << "---------********-----------********-------------" << endl;
      cout << " Search name: " << searchname << endl;
      cout << "---------********-----------********-------------" << endl;
      
      if( !gSystem->AccessPathName( searchname.Data() ) ){
	segmentexists = true;
	break;
      }
    }
   
    if( segmentexists ){
      new( (*filelist)[segcounter] ) THaRun( pathlist, codafilename.Data(), "GMN run" );
      cout << "Added segment " << segment << ", CODA file name = " << codafilename << endl;

      //( (THaRun*) (*filelist)[segcounter] )->SetDate( now );

      //if( stream == 0 && segment == 0 ){
      //( (THaRun*) (*filelist)[segcounter] )->SetDataRequired(THaRunBase::kDate|THaRunBase::kRunNumber);
      //( (THaRun*) (*filelist)[segcounter] )->Init();

      //RunDate = ( (THaRun*) (*filelist)[segcounter] )->GetDate();
      //} else {
      //	( (THaRun*) (*filelist)[segcounter] )->SetDataRequired(0);

      //cout << "Warning: setting date to " << RunDate.AsString() << " for stream " << stream << " segment " << segment 
      //     << endl; 

      //	( (THaRun*) (*filelist)[segcounter] )->SetDate(RunDate);
      //	( (THaRun*) (*filelist)[segcounter] )->SetNumber(UInt_t(runnum));
      //}
      //( (THaRun*) (*filelist)[segcounter] )->SetNumber( runnum );
      //( (THaRun*) (*filelist)[segcounter] )->Init();
    } // else {
    //   THaRun *rtemp = ( (THaRun*) (*filelist)[segcounter-1] ); //make otherwise identical copy of previous run in all respects except coda file name:
    //   new( (*filelist)[segcounter] ) THaRun( *rtemp );
    //   ( (THaRun*) (*filelist)[segcounter] )->SetFilename( codafilename.Data() );
    //   ( (THaRun*) (*filelist)[segcounter] )->SetNumber( runnum );
    //   cout << "Added segment " << segcounter << ", CODA file name = " << codafilename << endl;
    // }
    if( segmentexists ){
      segcounter++;
      lastsegment = segment;
    }
    segment++;
  }

  cout << "n segments to analyze = " << segcounter << endl;

  prefix = gSystem->Getenv("OUT_DIR");

  TString outfilename;
  outfilename.Form( "%s/%s_SBSGEMs_%d.evio.0.%d_%d.root", prefix.Data(), fname_prefix, runnum,
		     firstsegment, lastsegment );
  /*TString fntest(fname_prefix);
  if( fntest == "eel125"){
      outfilename.Form( "%s/%s_replayed_%d_newDB.root", prefix.Data(), fname_prefix, runnum);
      }*/
  
 TString fntest(fname_prefix);
  if( fntest == "eel125"){
    outfilename.Form( "%s/%s_replayed_%d_seg%d_%d.root", prefix.Data(), fname_prefix, runnum, firstsegment, lastsegment );
    } 
  if( parallelmode == 1 && fntest == "eel125"){

    if( singlemode == 1 ){
      outfilename.Form("%s/%s_replayed_%d_seg_%d.root", prefix.Data(), fname_prefix, runnum, file_seg);
    }
    else{
	    outfilename.Form("%s/%s_replayed_%d_seg_%d.root", prefix.Data(), fname_prefix, runnum, firstsegment);
    }	
  }

  if( fntest == "e1209016"){
    outfilename.Form( "%s/%s_SBSGEMs_%d.evio.%d_%d.root", prefix.Data(), fname_prefix, runnum, firstsegment, lastsegment );
    } 
  if( parallelmode == 1 && fntest == "e1209016"){

    if( singlemode == 1 ){
      outfilename.Form("%s/%s_SBSGEMs_%d.evio.0.%d.root", prefix.Data(), fname_prefix, runnum, file_seg);
    }
    else{
      outfilename.Form("%s/%s_SBSGEMs_%d.evio.0.%d.root", prefix.Data(), fname_prefix, runnum, firstsegment);
    } 
  }

  // Define the run(s) that we want to analyze.
  // We just set up one, but this could be many.
  //  THaRun* run = new THaRun( "prod12_4100V_TrigRate25_4.dat" );
  //THaRun* run = new THaRun( "5GEM_sample.dat" );
  //THaRun* run = new THaRun( "/Users/puckett/WORK/GEM_ALIGNMENT/RAWDATA/gem_cleanroom_2811.evio.31" );
  //THaRun* run = new THaRun( codafilename.Data() );
  //THaRun* run = new THaRun( "/Users/puckett/WORK/GEM_ALIGNMENT/RAWDATA/gem_cleanroom_2805.evio.0" );

  

  analyzer->SetVerbosity(2);
  analyzer->SetMarkInterval(100);

  analyzer->EnableBenchmarks();
  
  // Define the analysis parameters
  analyzer->SetEvent( event );
  analyzer->SetOutFile( outfilename.Data() );
  // File to record cuts accounting information
  analyzer->SetSummaryFile("replay_BBGEM.log"); // optional

  prefix = gSystem->Getenv("SBS_REPLAY");
  prefix += "/replay/";

  TString odef_filename = "replay_SBSGEM.odef";
  
  odef_filename.Prepend( prefix );

  analyzer->SetOdefFile( odef_filename );
  
  //analyzer->SetCompressionLevel(0); // turn off compression

  filelist->Compress();

  for( int iseg=0; iseg<filelist->GetEntries(); iseg++ ){
    THaRun *run = ( (THaRun*) (*filelist)[iseg] );

    run->Init();

    if( nevents > 0 ) run->SetLastEvent(nevents); //not sure if this will work as we want it to for multiple file segments chained together

    run->SetFirstEvent( firstevent );
    
    run->SetDataRequired(THaRunBase::kDate|THaRunBase::kRunNumber);

    // run->SetDataRequired(0);
    // run->SetDate(now);
    
    if( run->GetSegment() >= firstsegment && run->GetSegment() - firstsegment < maxsegments ){
      analyzer->Process(run);     // start the actual analysis
    }
  }
}
