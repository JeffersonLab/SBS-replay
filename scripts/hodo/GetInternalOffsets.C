#include "gmn_tree.C"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCut.h"
#include <iostream>
#include "TString.h"
#include "TCanvas.h"
#include <vector>
#include "TPad.h"
#include "TSystem.h"
#include "TLine.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TMinuit.h"
#include "TSpectrum.h"
#include "TClonesArray.h"


void GetInternalOffsets( const char *configfilename, const char *outputfilename="HodoInternalOffsets_temp.root", int refID = 40, double reft0 = 0.0, double etof0=15.1667 ){

  TFile *fout = new TFile(outputfilename,"RECREATE");

  ifstream configfile(configfilename);
  if( !configfile ) return;

  TString currentline;

  TChain *C = new TChain("T");

  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add( currentline.Data() );
    }
  }

  //TFile *fout = new TFile(outputfilename,"RECREATE");


  //  C->Add(rootfilename);

  gmn_tree *T = new gmn_tree(C);

  //  double bunchspace = 1.0e9/RF_freq; // ns

  //Total time range is going to be 3*bunchspace, we want about 8 bins/ns so

  //int nbins = ;
  //int nbins=100;
  //  TH2D *hdT_hodo_ID = new TH2D("hdT_hodo_ID", "; hodoscope bar ID ; t_{mean}^{raw}-t_{RF} (ns)", 90,-0.5,89.5, nbins, -1.5*bunchspace, 1.5*bunchspace);

  long nevent = 0;

  C->SetBranchStatus("*",0);
  C->SetBranchStatus("bb.tr.*",1);
  C->SetBranchStatus("bb.hodotdc.*",1);
  C->SetBranchStatus("bb.tdctrig.*",1);
  C->SetBranchStatus("bb.ps.*",1);
  C->SetBranchStatus("bb.sh.*",1);

  //double wL[90];
  //double wR[90];
  //double t0[90];
  //double vscint[90];

  double walk_default = (0.2113+0.2794)/2.0;
  double vscint_default = 0.162;

  vector<double> wL(90,walk_default);
  vector<double> wR(90,walk_default);
  vector<double> t0(90,0.0);
  vector<double> vscint(90,vscint_default);

  //Hodoscope internal alignment does not use nor require RF corrections:

  //Read parameters from initial fit:
  ifstream walkfile("hodowalkcorr_default.dat");
  ifstream vscintfile("hodovscint_default.dat");
  ifstream t0file("hodomeantime_default.dat");

  if( walkfile && vscintfile && t0file ){
    for( int i=0; i<90; i++ ){
      walkfile >> wL[i];
      vscintfile >> vscint[i];
      t0file >> t0[i];
      if( wL[i] < 0 ) wL[i] *= -1.; //the convention in this code requires the walk correction slope to be positive
    }

    for( int i=0; i<90; i++ ){
      walkfile >> wR[i];
      if( wR[i] < 0 ) wR[i] *= -1.;
    }
  }


  // Determine mean time offsets in an (RF-corrected) way;
  // Meantimes are already corrected for time walk
  // The idea is to look at (raw) mean-time differences between neighboring bars in the same cluster.
  // Minimize sum of squared corrected time differences between
  // chi^2 = sum_{i=1}^Nevent [sum_{j=1}^{Nhits_i-1} [sum_{k=j+1}^N_hits_i (t_j^i-t_k^i)^2]]
  // Where t_j^i,corr = t_mean_j^i - t0ji
  // Then dchi^2 / dt0_m = 0 = 2 * sum_i,j,k ( tmean_j^i - tmean_k^i - t0ji + t0ki ) * ( delta mk - delta mj )
  // = 2 * [ sum_i,j (tji-tmi + -t0ji + t0mi) ] - sum_i,k (t_mi - t_ki - t0mi + t0ki ) ] = 0
  // coefficient of
  TMatrixD M(90,90);
  TVectorD b(90);

  for( int i=0; i<90; i++ ){
    b(i) = 0.0;
    for( int j=0; j<90; j++ ){
      M(i,j) = 0.0;
      M(j,i) = 0.0;
    }
  }

  double zhodo = 1.854454; //meters
  double Lbar_hodo = 0.6; //meters

  vector<double> nevent_by_bar(90, 0.0);

  while( C->GetEntry(nevent++)){

    if( nevent % 10000 == 0 ) cout << "event " << nevent << endl;

    if( T->bb_tr_n >= 1. && T->bb_hodotdc_clus_trackindex[0]==0. && T->bb_ps_e>0.2 && T->bb_tdctrig_rftime != 0. && T->bb_tdctrig_rftime < 1.e10 && T->bb_tdctrig_trigtime < 1.e10 ){
      double RFtime = T->bb_tdctrig_rftime + T->bb_tr_vz[0]/0.299792458;

      //RFtime = fmod( RFtime, 160.321 );
      double trigtime = T->bb_tdctrig_trigtime;

      //grab needed track parameters:
      double vz = T->bb_tr_vz[0];
      double pathl = T->bb_tr_pathl[0];
      double etof = pathl/0.299792458-etof0; //electron TOF from vertex to hodo.
      double yhodo = T->bb_tr_y[0]+zhodo*T->bb_tr_ph[0];
      double xhodo = T->bb_tr_x[0]+zhodo*T->bb_tr_th[0];

      double dLEFT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 - yhodo));
      double dRIGHT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 + yhodo));

      double tmean0 = 0.0;

      for( int ibar=0; ibar<T->bb_hodotdc_clus_size[0]; ibar++ ){
	//double tHODO = T->bb_hodotdc_clus_bar_tdc_meantime[ibar];
	double tleft = T->bb_hodotdc_clus_bar_tdc_tleft[ibar];
	double tright = T->bb_hodotdc_clus_bar_tdc_tright[ibar];
	double totleft = T->bb_hodotdc_clus_bar_tdc_totleft[ibar];
	double totright = T->bb_hodotdc_clus_bar_tdc_totright[ibar];

	//double tmean = T->bb_hodotdc_clus_bar_tdc_meantime[ibar];
	//double ypos = T->bb_hodotdc_clus_bar_tdc_timehitpos[ibar];
	double xpos = T->bb_hodotdc_clus_bar_tdc_vpos[ibar];

	double IDHODO = T->bb_hodotdc_clus_bar_tdc_id[ibar];

	double tleft_corr = tleft - etof + wL[IDHODO]*totleft - t0[IDHODO] - dLEFT/vscint[IDHODO];
	double tright_corr = tright - etof + wR[IDHODO]*totright - t0[IDHODO] - dRIGHT/vscint[IDHODO];

	double tHODO = 0.5*(tleft_corr + tright_corr);

	//t0 and tof cancel in L/R difference:
	double tdiff = tleft + wL[IDHODO]*totleft - (tright + wR[IDHODO]*totright);
	//tL - tR = -2y/vscint
	//Recalculate ypos based on fit results:
	double ypos = -0.5*vscint[IDHODO]*tdiff;

	//	double deltat = tHODO - RFtime;
	//deltat -= bunchspace * ( floor(deltat/bunchspace) + 0.5 );

	double t_i = tHODO;
	if( ibar == 0 ){
	  tmean0 = t_i;
	}

	bool goodhit_i = totleft >= 6.0 && totleft <= 32.0 &&
	  totright >= 6.0 && totright <= 32.0 &&
	  fabs(totleft - totright) <= 8.0 &&
	  fabs(ypos - yhodo) <= 0.12 &&
	  fabs(xpos - xhodo) <= 0.055 &&
	  fabs(tHODO-tmean0)<=5.5;

	if( goodhit_i ){
	  for( int jbar=ibar+1; jbar<T->bb_hodotdc_clus_size[0]; jbar++ ){
	    double tleft_j = T->bb_hodotdc_clus_bar_tdc_tleft[jbar];
	    double tright_j = T->bb_hodotdc_clus_bar_tdc_tright[jbar];
	    double totleft_j = T->bb_hodotdc_clus_bar_tdc_totleft[jbar];
	    double totright_j = T->bb_hodotdc_clus_bar_tdc_totright[jbar];

	    //double tmean = T->bb_hodotdc_clus_bar_tdc_meantime[jbar];
	    //double ypos_j = T->bb_hodotdc_clus_bar_tdc_timehitpos[jbar];
	    double xpos_j = T->bb_hodotdc_clus_bar_tdc_vpos[jbar];

	    double IDHODO_j = T->bb_hodotdc_clus_bar_tdc_id[jbar];

	    double tleft_corr_j = tleft_j - etof + wL[IDHODO_j]*totleft_j - t0[IDHODO_j] - dLEFT/vscint[IDHODO_j];
	    double tright_corr_j = tright_j - etof + wR[IDHODO_j]*totright_j - t0[IDHODO_j] - dRIGHT/vscint[IDHODO_j];

	    double tHODO_j = 0.5*(tleft_corr_j + tright_corr_j);
	    double tdiff_j = tleft_j + wL[IDHODO_j]*totleft_j - (tright_j + wR[IDHODO_j]*totright_j);
	    double ypos_j = -0.5*vscint[IDHODO_j]*tdiff_j;

	    double t_j = tHODO_j;
	    bool goodhit_j = totleft_j >= 6.0 && totleft_j <= 32.0 &&
	      totright_j >= 6.0 && totright_j <= 32.0 &&
	      fabs(totleft_j - totright_j) <= 8.0 &&
	      fabs(ypos_j - yhodo) <= 0.12 &&
	      fabs(xpos_j - xhodo) <= 0.055 &&
	      fabs(tHODO_j-tmean0)<=5.5;

	    //t_corr = tmean - t0;
	    // t_corr_j - t_corr_i = (tmeanj - t0j - tmeani + t0i) -->
	    // we want RHS = tj-ti = t0j - t0i = LHS
	    if( goodhit_j ){
	      nevent_by_bar[IDHODO] += 1.0;
	      nevent_by_bar[IDHODO_j] += 1.0;

	      b(IDHODO) += -1.0*(t_j - t_i);
	      b(IDHODO_j) += 1.0*(t_j - t_i);
	      M(IDHODO,IDHODO) += 1.0;
	      M(IDHODO_j,IDHODO_j) += 1.0;
	      M(IDHODO,IDHODO_j) += -1.0;
	      M(IDHODO_j,IDHODO) += -1.0;
	    }
	  }
	}

	//double deltat = fmod(tHODO - (RFtime - trigtime - bunchspace * floor( ( RFtime-trigtime ) / bunchspace ) ), bunchspace);

	// hdT_hodo_ID->Fill( IDHODO, deltat );

	// hdT_hodo_ID->Fill( IDHODO, deltat - bunchspace );
	// hdT_hodo_ID->Fill( IDHODO, deltat + bunchspace );
	// for( int jbar=ibar+1; jbar<T->bb_hodotdc_clus_size[0]; jbar++ ){
	//   double tj = T->bb_hodotdc_clus_bar_tdc_meantime[jbar];
	//   double IDj = T->bb_hodotdc_clus_bar_tdc_id[jbar];

	//   double deltat_ij = tHODO - tj;

	//   //coefficient of t0i in equation for t0i
	//   // two cases: in both we have (t_ji - t_ki = t_0ji - t_0ki)
	//   // if m = j, coefficient of t_0ji in m equation is +1
	//   // OTOH if m = k, coefficient of t_0ki in m equation is -1
	//   M(IDHODO,IDj) += 1.0;
	//   M(IDj,IDHODO) += 1.0;
	//   M(IDHODO,IDHODO) += -1.0;
	//   M(IDj,IDj) += -1.0;
	//   b(IDHODO) += deltat_ij;
	//   b(IDj) += deltat_ij;
	// }
      }
    }
  }

  for( int i=0; i<90; i++ ){
    if( nevent_by_bar[i] <= 100. ){
      cout << "Warning: bar " << i << " fewer than 100 events... zeroing out offsets"
	   << endl;
      b(i) = 0.0;
      M(i,i) = 1.0;
      for( int j=0; j<90; j++ ){
	if( j != i ){
	  M(i,j) = M(j,i) = 0.0;
	}
      }
    }
  }

  TDecompSVD A(M);
  A.Solve(b);

  b.Print();

  double corr = reft0 - b(refID);

  for( int i=0; i<90; i++ ){
    // we WANT b(refID) = reft0
    // b(refID) += reft0 - b(refID)
    b(i) += corr;
    t0[i] += b(i);
  }

  TString dbfilename = outputfilename;
  dbfilename.ReplaceAll(".root",".dat");

  ofstream dbfile( dbfilename.Data() );

  dbfile << "bb.hodotdc.meantime_offset = " << endl;
  for( int i=0; i<90; i++ ){
    if( i>0 && i%9 == 0 ){
      dbfile << endl;
    }
    TString entry;
    dbfile << entry.Format("%15.6g", t0[i]) << " ";

  }

  //Do a second loop over the data and fill some histos:
  nevent = 0;

  TH1D *hdt_ji = new TH1D("hdt_ji", "All hits corrected time difference; t_{j}-t_{i} (ns);",600,-30,30);

  TH2D *hdt_ji_IDj = new TH2D("hdt_ji_IDj", "All hits; ID j; t_{j}-t_{i} (ns);",90,-0.5,89.5,600,-30,30);
  TH2D *hdt_ji_IDi = new TH2D("hdt_ji_IDi", "All hits; ID i; t_{j}-t_{i} (ns);",90,-0.5,89.5,600,-30,30);

  while( C->GetEntry(nevent++)){

    if( nevent % 10000 == 0 ) cout << "event " << nevent << endl;

    if( T->bb_tr_n >= 1. && T->bb_hodotdc_clus_trackindex[0]==0. && T->bb_ps_e>0.2 && T->bb_tdctrig_rftime != 0. && T->bb_tdctrig_rftime < 1.e10 && T->bb_tdctrig_trigtime < 1.e10 ){
      double RFtime = T->bb_tdctrig_rftime + T->bb_tr_vz[0]/0.299792458;

      //RFtime = fmod( RFtime, 160.321 );
      double trigtime = T->bb_tdctrig_trigtime;

      //grab needed track parameters:
      double vz = T->bb_tr_vz[0];
      double pathl = T->bb_tr_pathl[0];
      double etof = pathl/0.299792458-etof0; //electron TOF from vertex to hodo.
      double yhodo = T->bb_tr_y[0]+zhodo*T->bb_tr_ph[0];
      double xhodo = T->bb_tr_x[0]+zhodo*T->bb_tr_th[0];

      double dLEFT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 - yhodo));
      double dRIGHT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 + yhodo));

      double tmean0 = 0.0;

      for( int ibar=0; ibar<T->bb_hodotdc_clus_size[0]; ibar++ ){
	//double tHODO = T->bb_hodotdc_clus_bar_tdc_meantime[ibar];
	double tleft = T->bb_hodotdc_clus_bar_tdc_tleft[ibar];
	double tright = T->bb_hodotdc_clus_bar_tdc_tright[ibar];
	double totleft = T->bb_hodotdc_clus_bar_tdc_totleft[ibar];
	double totright = T->bb_hodotdc_clus_bar_tdc_totright[ibar];

	//double tmean = T->bb_hodotdc_clus_bar_tdc_meantime[ibar];
	//double ypos = T->bb_hodotdc_clus_bar_tdc_timehitpos[ibar];
	double xpos = T->bb_hodotdc_clus_bar_tdc_vpos[ibar];

	double IDHODO = T->bb_hodotdc_clus_bar_tdc_id[ibar];

	double tleft_corr = tleft - etof + wL[IDHODO]*totleft - t0[IDHODO] - dLEFT/vscint[IDHODO];
	double tright_corr = tright - etof + wR[IDHODO]*totright - t0[IDHODO] - dRIGHT/vscint[IDHODO];

	double tHODO = 0.5*(tleft_corr + tright_corr);

	//t0 and tof cancel in L/R difference:
	double tdiff = tleft + wL[IDHODO]*totleft - (tright + wR[IDHODO]*totright);
	//tL - tR = -2y/vscint
	//Recalculate ypos based on fit results:
	double ypos = -0.5*vscint[IDHODO]*tdiff;

	//double deltat = tHODO - RFtime;
	//deltat -= bunchspace * ( floor(deltat/bunchspace) + 0.5 );

	double t_i = tHODO;
	if( ibar == 0 ){
	  tmean0 = t_i;
	}



	//	if( goodhit_i ){
	for( int jbar=ibar+1; jbar<T->bb_hodotdc_clus_size[0]; jbar++ ){
	  double tleft_j = T->bb_hodotdc_clus_bar_tdc_tleft[jbar];
	  double tright_j = T->bb_hodotdc_clus_bar_tdc_tright[jbar];
	  double totleft_j = T->bb_hodotdc_clus_bar_tdc_totleft[jbar];
	  double totright_j = T->bb_hodotdc_clus_bar_tdc_totright[jbar];

	  //double tmean = T->bb_hodotdc_clus_bar_tdc_meantime[jbar];
	  //double ypos_j = T->bb_hodotdc_clus_bar_tdc_timehitpos[jbar];
	  double xpos_j = T->bb_hodotdc_clus_bar_tdc_vpos[jbar];

	  double IDHODO_j = T->bb_hodotdc_clus_bar_tdc_id[jbar];

	  double tleft_corr_j = tleft_j - etof + wL[IDHODO_j]*totleft_j - t0[IDHODO_j] - dLEFT/vscint[IDHODO_j];
	  double tright_corr_j = tright_j - etof + wR[IDHODO_j]*totright_j - t0[IDHODO_j] - dRIGHT/vscint[IDHODO_j];

	  double tHODO_j = 0.5*(tleft_corr_j + tright_corr_j);
	  double tdiff_j = tleft_j + wL[IDHODO_j]*totleft_j - (tright_j + wR[IDHODO_j]*totright_j);
	  double ypos_j = -0.5*vscint[IDHODO_j]*tdiff_j;

	  double t_j = tHODO_j;

	  hdt_ji->Fill( t_j - t_i );
	  hdt_ji_IDj->Fill( IDHODO_j, t_j - t_i );
	  hdt_ji_IDi->Fill( IDHODO, t_j - t_i );
	}
      }
    }
  }

  fout->Write();

}
