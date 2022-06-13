void GEM_hit_map(int layer){

  gStyle->SetPalette(1);

  double width = .60;
  double height = 2.0;
  double edge = 0.05;

  if(layer < 4){
    width = 0.4;
    height = 1.5;
  }

  TTree *t = (TTree*) gFile->Get("T");

  TH2F *hit2d = new TH2F(Form("hit2d_layer%d",layer),Form("Layer %i;y (m); x (m)",layer),120,-width/2 - edge,width/2 + edge,120,-height/2 - edge,height/2 + edge);

  
  t->Draw( Form("bb.gem.hit.xglobal:bb.gem.hit.yglobal>>hit2d_layer%d",layer),Form("bb.gem.hit.layer == %d&& bb.gem.track.ntrack==1&&bb.gem.track.nhits[bb.gem.hit.trackindex]>=3",layer),"colz" );

  gStyle->SetOptStat(0);  
 
}
