void hcalTrigDiff(){

  cout << "Processing hcalTrigDiff_trig.." << endl;

  TStopwatch *st = new TStopwatch();
  st->Start(kTRUE);

  Int_t Ndata_bb_tdctrig_tdcelemID = 0;
  Double_t bb_sh_nclus = 0., sbs_hcal_nclus = 0.;
  Double_t bb_tdctrig_tdc[6] = {0.}, bb_tdctrig_tdcelemID[6] = {0.};

  TH1D *h1_hcalTrigDiff = new TH1D("h1_hcalTrigDiff","HCal/BBCalHi Difference; ns", 1000, -50, 950);

  // Declare branches
  TTree *T = (TTree*) gDirectory->Get("T");
  T->SetBranchStatus("*",0);
  T->SetBranchStatus("bb.tdctrig.tdc",1);
  T->SetBranchStatus("bb.sh.nclus",1);
  T->SetBranchStatus("sbs.hcal.nclus",1);
  T->SetBranchStatus("Ndata.bb.tdctrig.tdcelemID",1);
  T->SetBranchStatus("bb.tdctrig.tdcelemID",1);

  T->SetBranchAddress("bb.tdctrig.tdc", &bb_tdctrig_tdc);
  T->SetBranchAddress("bb.tdctrig.tdcelemID", &bb_tdctrig_tdcelemID);
  T->SetBranchAddress("Ndata.bb.tdctrig.tdcelemID", &Ndata_bb_tdctrig_tdcelemID);
  T->SetBranchAddress("bb.sh.nclus", &bb_sh_nclus);
  T->SetBranchAddress("sbs.hcal.nclus", &sbs_hcal_nclus);

  // Loop through events
  Long64_t nevents = T->GetEntries();
  for(Long64_t nevent=0; nevent<nevents; nevent++){

    T->GetEntry(nevent);

    if( bb_sh_nclus==0 || sbs_hcal_nclus==0 ) continue;

    Double_t bb_time=0.; //L1A
    Double_t hcal_time=0.; //Overlapping regions
    for(Int_t ihit=0; ihit<Ndata_bb_tdctrig_tdcelemID; ihit++){
      if(bb_tdctrig_tdcelemID[ihit]==0) hcal_time=bb_tdctrig_tdc[ihit];
      if(bb_tdctrig_tdcelemID[ihit]==5) bb_time=bb_tdctrig_tdc[ihit];
    }
    h1_hcalTrigDiff->Fill(hcal_time-bb_time);
  }
  
  h1_hcalTrigDiff->SetStats(0);
  h1_hcalTrigDiff->Draw();

  cout << "Processed macro with " << nevents << " entries." << endl;
  cout << "Note: A large spread here indicates that BBCal is setting the trigger time!" << endl;

  st->Stop();
  cout << "CPU time = " << st->CpuTime() << " s " << " Real time = " << st->RealTime() << " s " << endl;
}
