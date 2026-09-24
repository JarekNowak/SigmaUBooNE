// Composition of the NuWro sample: is it Sigma0-filtered? If so its Lambda
// content is suppressed by construction and cannot be compared to a Lambda
// CCQE prediction.
#include "HypInput.h"
#include <fstream>
#include <map>
void NuwroComp(int maxFiles=300){
  TChain ch("ana/OutputTree");
  std::ifstream in("lists/sig_nuwro.list"); std::string l; int k=0;
  while(std::getline(in,l)) if(l.size()&&(maxFiles<=0||k<maxFiles)){ch.Add(l.c_str());k++;}
  HypInput H; if(!H.attach(&ch)){printf("no attach\n");return;}
  std::map<int,long> hyp; long N=0,lam=0,sig=0;
  for(Long64_t i=0;i<ch.GetEntries();++i){ ch.GetEntry(i); H.convert(); N++;
    hyp[H.h_hyperon_pdg]++;
    if(H.nu_pdg==-14&&H.ccnc==0){ if(H.h_hyperon_pdg==3122) lam++; if(H.h_hyperon_pdg==3212) sig++; } }
  printf("\n  NuWro sample: %lld events in %d files\n",(Long64_t)N,k);
  for(auto&kv:hyp) printf("    pdg %6d : %7ld  (%.1f%%)\n",kv.first,kv.second,100.*kv.second/N);
  printf("  Lambda CC %ld , Sigma0 CC %ld  -> ratio %.2f\n",lam,sig,sig?double(lam)/sig:0);
}
