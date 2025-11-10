

void coin_time(const std::string rootfile) {
  TFile *f = new TFile(rootfile.c_str(), "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: Bad file!" << std::endl;
    return;
      }

  TH2D *hdt_HCAL_hodo_vs_IDHCAL = (TH2D*) f->Get("hdt_HCAL_hodo_vs_IDHCAL");
  if (!hdt_HCAL_hodo_vs_IDHCAL) {
    std::cerr << "ERROR: Bad histogram name: hdt_HCAL_hodo_vs_IDHCAL" << std::endl;
    return;
  }
  
  TGraphErrors *gT0HCAL = (TGraphErrors*) f->Get("gT0_HCAL");
  if (!gT0HCAL) {
    std::cerr << "ERROR: Bad TGraphErrors name: gT0_HCAL" << std::endl;
    return;
  }

  hdt_HCAL_hodo_vs_IDHCAL->SetDirectory(0);
  f->Close();

  TH2D *hdt_HCAL_hodo_vs_IDHCAL_corr = (TH2D*) hdt_HCAL_hodo_vs_IDHCAL->Clone("hdt_HCAL_hodo_vs_IDHCAL_corr");
  hdt_HCAL_hodo_vs_IDHCAL_corr->Reset();

  int nx = hdt_HCAL_hodo_vs_IDHCAL_corr->GetNbinsX();
  int ny = hdt_HCAL_hodo_vs_IDHCAL_corr->GetNbinsY();
  int n = gT0HCAL->GetN();

  double meanOffset, xjunk;
  for (int i = 0; i < n; i++) {
    gT0HCAL->GetPoint(i, xjunk, meanOffset);
    meanOffset += meanOffset / n;
  }
  

  for (int i = 1; i <= nx; i++) {
    double x = hdt_HCAL_hodo_vs_IDHCAL->GetXaxis()->GetBinCenter(i);
    double y_shift = gT0HCAL->Eval(x);
    // std::cout << y_shift << std::endl;

    for (int j = 1; j <= ny; j++) {
      double y = hdt_HCAL_hodo_vs_IDHCAL->GetYaxis()->GetBinCenter(j);
      double content = hdt_HCAL_hodo_vs_IDHCAL->GetBinContent(i,j);

      double y_corr = y - y_shift + 120.0;

      hdt_HCAL_hodo_vs_IDHCAL_corr->Fill(x,y_corr,content);
      
    }
    
  }

  TCanvas *c1 = new TCanvas("c1", "2D Histogram Correction", 800, 600);
  c1->cd();
  hdt_HCAL_hodo_vs_IDHCAL_corr->Draw("COLZ");

  TCanvas *c2 = new TCanvas("c2", "Proj", 800, 600);
  c2->cd();
  TH1D *hproj = hdt_HCAL_hodo_vs_IDHCAL_corr->ProjectionY("hproj");
  hproj->Draw();
  
  
}
