#include <fstream>
#include <sstream>
#include <cmath>

void CoinTimeAnalysis() {

  std::string rootfiles = "/volatile/halla/sbs/koeneman/replays/He3/GEN4b/rootfiles/e1209016_fullreplay*.root";

  TChain *C = new TChain("T");
  C->Add(rootfiles.c_str());
  
  C->SetBranchStatus("*",0);
  C->SetBranchStatus("bb.sh.atimeblk",1);
  C->SetBranchStatus("bb.sh.idblk",1);
  C->SetBranchStatus("sbs.hcal.atimeblk",1);
  C->SetBranchStatus("sbs.hcal.idblk",1);
  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus("g.*",1);
  C->SetBranchStatus("bb.ps.e",1);
  C->SetBranchStatus("bb.tr.*",1);
  C->SetBranchStatus("bb.etot_over_p",1);
  C->SetBranchStatus("e.*",1);

  double bb_sh_atimeblk;
  double sbs_hcal_atimeblk;
  double bb_sh_idblk;
  double sbs_hcal_idblk;
  C->SetBranchAddress("bb.sh.atimeblk", &bb_sh_atimeblk);
  C->SetBranchAddress("bb.sh.idblk", &bb_sh_idblk);
  C->SetBranchAddress("sbs.hcal.atimeblk", &sbs_hcal_atimeblk);
  C->SetBranchAddress("sbs.hcal.idblk", &sbs_hcal_idblk);

  double trn, pse, hcale, vz[100], eop[100], trigbits, w2;
  C->SetBranchAddress("bb.tr.n", &trn);
  C->SetBranchAddress("bb.ps.e", &pse);
  C->SetBranchAddress("sbs.hcal.e", &hcale);
  C->SetBranchAddress("bb.tr.vz", vz);
  C->SetBranchAddress("g.trigbits", &trigbits);
  C->SetBranchAddress("bb.etot_over_p", eop);
  C->SetBranchAddress("e.kine.W2", &w2);

  std::string HCAL_ADCtime_offsets = "HCALt0_GEN4_He3_TOFcal_temp.dat";
  std::string BBSH_ADCtime_offsets = "BBSHt0_GEN4_He3_TOFcal_temp.dat";
  
  int nchan_HCAL = 288;
  int nchan_BBSH = 189;

  std::string line;

  bool found = false;
  double HCALblk[0];
  for (int i = 0; i < nchan_HCAL; i++) {
    HCALblk[i] = i;
  }
  vector<double> HCALoffsets;
  std::ifstream HCALoffsetFile(HCAL_ADCtime_offsets);
  while (std::getline(HCALoffsetFile, line)) {
    if (line.find("sbs.hcal.adc.timeoffset") != std::string::npos) {
      found = true;
      continue;
    }

    if (found) {
      std::stringstream ss(line);
      double val;
      while (ss >> val) {
	HCALoffsets.push_back(val);
	//std::cout << val << std::endl;
      }
    }
  }

  std::cout << HCALoffsets.size() << std::endl;

  found = false;
  double BBSHblk[500];
  for (int i = 0; i < nchan_BBSH; i++) {
    BBSHblk[i] = i;
  }
  vector<double> BBSHoffsets;
  std::ifstream BBSHoffsetFile(BBSH_ADCtime_offsets);
  while (std::getline(BBSHoffsetFile, line)) {
    if (line.find("bb.sh.adc.timeoffset") != std::string::npos) {
      found = true;
      continue;
    }
    
    if (found) {
      std::stringstream ss(line);
      double val;
      while (ss >> val) {
	BBSHoffsets.push_back(val);
	//std::cout << val << std::endl;
      }
    }
  }

  std::cout << BBSHoffsets.size() << std::endl;

  TH1D *cointime_uncorr = new TH1D("cointime_uncorr", "SH HCAL FADC cointime uncorr; t_{BBSH}^{(FADC)}-t_{HCAL}^{(FADC)} (ns); counts", 300, -100.0, 100.0);
  TH1D *cointime = new TH1D("cointime", "SH HCAL FADC cointime;t_{BBSH corr}^{(FADC)}-t_{HCAL corr}^{(FADC)} (ns); counts", 300, -100.0, 100.0);

  /*
  double BBSHoffsetsgraph[500], HCALoffsetsgraph[500];
  for (int i = 0; i < nchan_BBSH; i++) {
    BBSHoffsetsgraph[i] = BBSHoffsets[i];
  }

  for (int i = 0; i < nchan_HCAL; i++) {
    HCALoffsetsgraph[i] = HCALoffsets[i];
  }
  
  TGraph *T0BBSH = new TGraph(189, BBSHoffsetsgraph);
  TGraph *T0HCAL = new TGraph(288, HCALoffsetsgraph);

  TCanvas *graph1 = new TCanvas("graph1", "BBSH offsets", 800, 600);
  TCanvas *graph2 = new TCanvas("graph2", "HCAL offsets", 800, 600);

  graph1->cd();
  T0BBSH->Draw("AP");
  graph2->cd();
  T0HCAL->Draw("AP");
  */
  Long64_t Nevents = C->GetEntries();
  std::cout << Nevents << std::endl;
  Long64_t nevents = 0;
  for (Long64_t i = 0; i < Nevents; i++) {
    C->GetEntry(i);

    int bb_sh_idblki = (int)bb_sh_idblk;
    
    int sbs_hcal_idblki = (int)sbs_hcal_idblk;

    double BBSHoffseti = BBSHoffsets[bb_sh_idblki];
    double HCALoffseti = HCALoffsets[sbs_hcal_idblki];

    
    bool cut = (
		trn>0.0 &&
		pse>0.2 &&
		fabs(vz[0])<0.27 &&
		fabs(eop[0]-1.0)<0.3 &&
		hcale>0.02 &&
		trigbits==4.0
		);

    if (!cut) continue;

    nevents += 1;

    double bb_sh_atimeblk_corri = bb_sh_atimeblk - BBSHoffseti;
    double sbs_hcal_atimeblk_corri = sbs_hcal_atimeblk - HCALoffseti;

    double cointimei = bb_sh_atimeblk_corri + sbs_hcal_atimeblk_corri - 300.0;
    double cointime_uncorri = bb_sh_atimeblk + sbs_hcal_atimeblk - 200.0;

    cointime->Fill(cointimei);
    cointime_uncorr->Fill(cointime_uncorri);
    
  }

  std::cout << nevents << std::endl;

  TCanvas *c1 = new TCanvas("c1", "Coin Time Histogram", 800, 600);
  c1->cd();
  cointime->Draw();

  TCanvas *c2 = new TCanvas("c2", "Coin Time Histogram", 800, 600);
  c2->cd();
  cointime_uncorr->Draw();

  delete C;
}
