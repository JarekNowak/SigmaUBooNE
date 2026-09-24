// What hyperon species does the GENIE run3b sample contain, and do we have
// enough direct Lambda CCQE events to reproduce their selection?
#include "HypInput.h"
#include <fstream>
#include <map>
void HypComp(int maxFiles=400){
  TChain ch("ana/OutputTree");
  std::ifstream in("lists/sig_genie.list"); std::string l; int k=0;
  while(std::getline(in,l)) if(l.size()&&(maxFiles<=0||k<maxFiles)){ch.Add(l.c_str());k++;}
  HypInput H; if(!H.attach(&ch)){printf("no attach\n");return;}
  std::vector<int>* dpdg=0; std::vector<double>* dmom=0;
  if(ch.GetBranch("mc_decay_pdg")) ch.SetBranchAddress("mc_decay_pdg",&dpdg);
  if(ch.GetBranch("mc_decay_mom")) ch.SetBranchAddress("mc_decay_mom",&dmom);
  std::map<int,long> hyp; long N=0, lamCC=0, lamPPi=0, lamThresh=0, sigCC=0;
  for(Long64_t i=0;i<ch.GetEntries();++i){ ch.GetEntry(i); H.convert(); N++;
    hyp[H.h_hyperon_pdg]++;
    const bool ccbar = (H.nu_pdg==-14 && H.ccnc==0);
    if(H.h_hyperon_pdg==3122 && ccbar){ lamCC++;
      bool hasP=false,hasPi=false; double pP=0,pPi=0;
      if(dpdg&&dmom) for(size_t j=0;j<dpdg->size()&&j<dmom->size();++j){
        if(dpdg->at(j)==2212){ hasP=true; pP=dmom->at(j); }
        if(dpdg->at(j)==-211){ hasPi=true; pPi=dmom->at(j); } }
      if(hasP&&hasPi){ lamPPi++;
        if(pP>0.300 && pPi>0.100) lamThresh++; } }
    if(H.h_hyperon_pdg==3212 && ccbar) sigCC++;
  }
  printf("\n  %lld events in %d files\n",(Long64_t)N,k);
  printf("  hyperon species:\n");
  for(auto&kv:hyp) printf("    pdg %6d : %7ld  (%.1f%%)\n",kv.first,kv.second,100.*kv.second/N);
  printf("\n  their signal definition, step by step:\n");
  printf("    Lambda + anti-numu CC            : %ld\n",lamCC);
  printf("    ... decaying to p pi-            : %ld  (%.1f%% of the above; BF is 63.9%%)\n",
         lamPPi, lamCC?100.*lamPPi/lamCC:0.);
  printf("    ... p>300 MeV/c and pi>100 MeV/c : %ld  (%.1f%%)\n",
         lamThresh, lamPPi?100.*lamThresh/lamPPi:0.);
  printf("\n  for scale, Sigma0 + anti-numu CC   : %ld\n",sigCC);
}
