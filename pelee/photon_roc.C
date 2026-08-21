// Signal vs background separation for the photon-candidate variables, using the
// actual selected samples (TreeS = 650 signal, TreeB = 69557 background). This is
// the quantity that matters for a cut, independent of the unreliable non-muon truth
// matching in the NoCR signal sample.
void photon_roc(){
  TFile f("BDTInput.root");
  TTree* S=(TTree*)f.Get("TreeS"); TTree* B=(TTree*)f.Get("TreeB");
  struct V { const char* br; const char* nm; int nb; double lo,hi; bool lowIsSig; };
  std::vector<V> vs = {
    {"Shower1NHitsY",   "NHits (Y)",         200, 0, 400,  false},
    {"Shower1VtxDist",  "vtx dist [cm]",     200, 0, 80,   false},
    {"PhotonPointAngle","point angle [rad]", 200, 0, 3.15, true },
    {"Shower1DeDxY",    "dE/dx [MeV/cm]",    200, 0, 15,   true },
    {"Shower1Energy",   "energy [GeV]",      200, 0, 1.0,  false},
    {"PhotonImpact",    "impact [cm]",       200, 0, 40,   true },
  };
  printf("\n  %-20s %8s %8s | %s\n","variable","sig mean","bkg mean","best cut: eff_S / rej_B / FoM gain");
  for(auto&v:vs){
    TH1D hs("hs","",v.nb,v.lo,v.hi), hb("hb","",v.nb,v.lo,v.hi);
    S->Draw(Form("%s>>hs",v.br),"","goff"); B->Draw(Form("%s>>hb",v.br),"","goff");
    double ns=hs.Integral(0,v.nb+1), nb=hb.Integral(0,v.nb+1);
    if(ns<=0||nb<=0){ printf("  %-20s (empty)\n",v.nm); continue; }
    // scan every threshold, keep the one maximising S/sqrt(S+B) relative to no cut
    double base = ns/std::sqrt(ns+nb), best=0, bcut=0, beffs=0, brejb=0;
    for(int i=1;i<=v.nb;++i){
      double x=hs.GetBinLowEdge(i);
      double s = v.lowIsSig ? hs.Integral(0,i) : hs.Integral(i,v.nb+1);
      double b = v.lowIsSig ? hb.Integral(0,i) : hb.Integral(i,v.nb+1);
      if(s<0.5*ns) continue;             // never give up more than half the signal
      if(s+b<=0) continue;
      double fom=s/std::sqrt(s+b);
      if(fom>best){ best=fom; bcut=x; beffs=s/ns; brejb=1.-b/nb; }
    }
    printf("  %-20s %8.2f %8.2f | %s %.2f : eff %5.1f%%  rej %5.1f%%  FoM x%.3f\n",
      v.nm, hs.GetMean(), hb.GetMean(), v.lowIsSig?"<":">", bcut, 100*beffs, 100*brejb, best/base);
  }
}
