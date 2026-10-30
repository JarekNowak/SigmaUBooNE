// DupGen.C -- count distinct generator events in a hypana list, keyed on (lepton momentum,
// hyperon momentum, Q2) to 1e-6. The previous GENIE hyperon sample (lists/sig_genie.list)
// gives 754 distinct numubar Sigma0 among 15,341 (~20 copies each); the QE sample none repeated.
// Run:  root -l -b -q 'DupGen.C("lists/sig_genie.list")'
void DupGen(const char* list, int maxFiles=0){
  std::ifstream in(list); std::string l; std::vector<std::string> F; while(std::getline(in,l)) if(l.size()&&(maxFiles<=0||(int)F.size()<maxFiles)) F.push_back(l);
  std::map<std::tuple<long,long,long>,int> seen; long n=0, nsig=0;
  for(auto& fn:F){ TFile f(fn.c_str()); TTree*t=(TTree*)f.Get("ana/OutputTree"); if(!t) continue;
    double pm,ph,q2; int hyp,nup; t->SetBranchStatus("*",0);
    for(const char*b:{"mc_lepton_mom","mc_hyperon_mom","mc_nu_q2","mc_hyperon_pdg","mc_nu_pdg"}) t->SetBranchStatus(b,1);
    t->SetBranchAddress("mc_lepton_mom",&pm); t->SetBranchAddress("mc_hyperon_mom",&ph); t->SetBranchAddress("mc_nu_q2",&q2);
    t->SetBranchAddress("mc_hyperon_pdg",&hyp); t->SetBranchAddress("mc_nu_pdg",&nup);
    for(Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i); n++; if(hyp!=3212||nup!=-14) continue; nsig++;
      seen[{std::lround(pm*1e6),std::lround(ph*1e6),std::lround(q2*1e6)}]++; } }
  std::map<int,int> mult; for(auto&kv:seen) mult[kv.second]++;
  printf("%s: %ld events, %ld numubar Sigma0, %zu distinct (%.2f copies each)\n",list,n,nsig,seen.size(),double(nsig)/seen.size());
  for(auto&kv:mult) printf("   seen %dx: %d\n",kv.first,kv.second);
}
