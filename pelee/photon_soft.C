// The Sigma0 -> Lambda gamma photon carries only 74.5 MeV in the Sigma0 rest frame,
// so the SIGNAL photon is soft and the pi0 photons that dominate the background are
// hard. Test cuts in that direction, and the two-sided windows a soft photon implies.
void photon_soft(){
  TFile f("BDTInput.root");
  TTree* S=(TTree*)f.Get("TreeS"); TTree* B=(TTree*)f.Get("TreeB");
  auto scan=[&](const char* br,const char* nm,int nb,double lo,double hi,bool lowIsSig){
    TH1D hs("hs","",nb,lo,hi), hb("hb","",nb,lo,hi);
    S->Draw(Form("%s>>hs",br),"","goff"); B->Draw(Form("%s>>hb",br),"","goff");
    double ns=hs.Integral(0,nb+1), nbk=hb.Integral(0,nb+1);
    double base=ns/std::sqrt(ns+nbk), best=0,bcut=0,be=0,br_=0;
    for(int i=1;i<=nb;++i){ double x=hs.GetBinLowEdge(i);
      double s=lowIsSig?hs.Integral(0,i):hs.Integral(i,nb+1);
      double b=lowIsSig?hb.Integral(0,i):hb.Integral(i,nb+1);
      if(s<0.5*ns||s+b<=0) continue;
      double fom=s/std::sqrt(s+b);
      if(fom>best){best=fom;bcut=x;be=s/ns;br_=1.-b/nbk;} }
    printf("  %-18s %s %6.2f : eff %5.1f%%  rej %5.1f%%  FoM x%.3f\n",
      nm, lowIsSig?"<":">", bcut, 100*be, 100*br_, best/base);
    return best/base;
  };
  printf("\n  cuts in the SOFT-photon direction (signal = low)\n");
  scan("Shower1Energy","energy [GeV]",300,0,1.0,true);
  scan("Shower1NHitsY","NHits (Y)",   400,0,400,true);
  printf("\n  two-sided windows (a soft photon still needs to be reconstructed)\n");
  auto win=[&](const char* br,const char* nm,double lo,double hi,int nb,double a,double b_){
    TH1D hs("hs2","",nb,a,b_), hb("hb2","",nb,a,b_);
    S->Draw(Form("%s>>hs2",br),"","goff"); B->Draw(Form("%s>>hb2",br),"","goff");
    double ns=hs.Integral(0,nb+1), nbk=hb.Integral(0,nb+1);
    double base=ns/std::sqrt(ns+nbk);
    int i0=hs.FindBin(lo), i1=hs.FindBin(hi);
    double s=hs.Integral(i0,i1), bb=hb.Integral(i0,i1);
    printf("  %-18s [%.3f,%.3f] : eff %5.1f%%  rej %5.1f%%  FoM x%.3f\n",
      nm, lo, hi, 100*s/ns, 100*(1.-bb/nbk), (s/std::sqrt(s+bb))/base);
  };
  for(double hi : {0.10,0.15,0.20,0.30}) win("Shower1Energy","energy window",0.005,hi,300,0,1.0);
  printf("\n  combined: soft photon AND pi0 veto\n");
  // pi0 veto: reject a reconstructed di-photon mass near 135 MeV
  auto comb=[&](const char* cut,const char* nm){
    double s=S->GetEntries(cut), b=B->GetEntries(cut);
    double ns=S->GetEntries(), nbk=B->GetEntries();
    double base=ns/std::sqrt(ns+nbk);
    printf("  %-46s eff %5.1f%%  rej %5.1f%%  FoM x%.3f\n",
      nm, 100*s/ns, 100*(1.-b/nbk), (s+b>0)?(s/std::sqrt(s+b))/base:0);
  };
  comb("Shower1Energy>0.005 && Shower1Energy<0.20","E in [5,200] MeV");
  comb("Pi0Mass<0 || Pi0Mass>0.16 || Pi0Mass<0.11","pi0 mass veto [110,160] MeV");
  comb("Shower1Energy>0.005 && Shower1Energy<0.20 && (Pi0Mass<0 || Pi0Mass>0.16 || Pi0Mass<0.11)","both");
}
