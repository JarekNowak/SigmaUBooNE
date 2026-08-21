// What is the Sigma0 photon candidate actually made of, and what separates a real
// photon from the impostors? Background events only: the NoCR signal sample's
// non-muon truth matching is unreliable (see the note).
void photon_study(){
  TFile f("BDTInput.root");
  TTree* B=(TTree*)f.Get("TreeB");
  float pdg, dedx, vtxd, nhit, ene, ang, imp, w;
  B->SetBranchAddress("Shr1BktPdg",&pdg);
  B->SetBranchAddress("Shower1DeDxY",&dedx);
  B->SetBranchAddress("Shower1VtxDist",&vtxd);
  B->SetBranchAddress("Shower1NHitsY",&nhit);
  B->SetBranchAddress("Shower1Energy",&ene);
  B->SetBranchAddress("PhotonPointAngle",&ang);
  B->SetBranchAddress("PhotonImpact",&imp);
  B->SetBranchAddress("Weight",&w);
  std::map<int,long> comp;
  Long64_t N=B->GetEntries();
  for(Long64_t i=0;i<N;++i){ B->GetEntry(i); comp[(int)pdg]++; }
  printf("\n  photon-candidate truth composition, background (%lld events)\n", N);
  std::vector<std::pair<long,int>> v; for(auto&p:comp) v.push_back({p.second,p.first});
  std::sort(v.rbegin(),v.rend());
  auto name=[](int p)->const char*{ switch(abs(p)){
    case 22:return "photon"; case 2212:return "proton"; case 211:return "pi+-";
    case 13:return "muon"; case 11:return "electron"; case 2112:return "neutron";
    case 321:return "K+-"; case 0:return "unmatched"; default:return "other";} };
  for(size_t k=0;k<v.size()&&k<8;++k)
    printf("    %-10s (%6d) : %7ld  %5.1f%%\n", name(v[k].second), v[k].second, v[k].first, 100.*v[k].first/N);

  // Separation power of each handle: true photon vs proton impostor
  struct H { const char* nm; TH1D* g; TH1D* p; };
  std::vector<H> hs = {
    {"dE/dx [MeV/cm]",   new TH1D("g0","",60,0,12),  new TH1D("p0","",60,0,12)},
    {"vtx dist [cm]",    new TH1D("g1","",60,0,60),  new TH1D("p1","",60,0,60)},
    {"NHits (Y)",        new TH1D("g2","",60,0,300), new TH1D("p2","",60,0,300)},
    {"energy [GeV]",     new TH1D("g3","",60,0,0.6), new TH1D("p3","",60,0,0.6)},
    {"point angle [rad]",new TH1D("g4","",60,0,1.6), new TH1D("p4","",60,0,1.6)},
    {"impact [cm]",      new TH1D("g5","",60,0,30),  new TH1D("p5","",60,0,30)},
  };
  for(Long64_t i=0;i<N;++i){ B->GetEntry(i);
    bool isg = (abs((int)pdg)==22), isp = (abs((int)pdg)==2212);
    if(!isg && !isp) continue;
    double vals[6]={dedx,vtxd,nhit,ene,ang,imp};
    for(int k=0;k<6;++k){ (isg?hs[k].g:hs[k].p)->Fill(vals[k]); } }
  printf("\n  true photon vs proton impostor: means and separation\n");
  printf("  %-20s %10s %10s %10s\n","variable","photon","proton","KS dist");
  for(auto&h:hs){
    if(h.g->Integral()<=0||h.p->Integral()<=0){ printf("  %-20s  (empty)\n",h.nm); continue; }
    h.g->Scale(1./h.g->Integral()); h.p->Scale(1./h.p->Integral());
    double ks=h.g->KolmogorovTest(h.p,"M");
    printf("  %-20s %10.3f %10.3f %10.3f\n", h.nm, h.g->GetMean(), h.p->GetMean(), ks);
  }
}
