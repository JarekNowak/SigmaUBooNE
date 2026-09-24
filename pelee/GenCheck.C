// Which generator produced the hypana signal sample? mc_mode carries the
// generator's own interaction-mode string, whose vocabulary differs between
// GENIE and NuWro.
#include "HypInput.h"
#include <fstream>
#include <map>
void GenCheck(int maxFiles=60){
  TChain ch("ana/OutputTree");
  std::ifstream in("lists/sig_genie.list"); std::string l; int k=0;
  while(std::getline(in,l)) if(l.size()&&(maxFiles<=0||k<maxFiles)){ch.Add(l.c_str());k++;}
  std::string* mode=0; if(ch.GetBranch("mc_mode")) ch.SetBranchAddress("mc_mode",&mode);
  std::string* ccnc=0; if(ch.GetBranch("mc_ccnc")) ch.SetBranchAddress("mc_ccnc",&ccnc);
  std::map<std::string,long> m, c;
  for(Long64_t i=0;i<ch.GetEntries();++i){ ch.GetEntry(i);
    if(mode) m[*mode]++; if(ccnc) c[*ccnc]++; }
  printf("\n  mc_mode values over %lld events in %d files:\n",(Long64_t)ch.GetEntries(),k);
  for(auto&kv:m) printf("    %-14s %7ld\n",kv.first.c_str(),kv.second);
  printf("  mc_ccnc values:\n");
  for(auto&kv:c) printf("    %-14s %7ld\n",kv.first.c_str(),kv.second);
}
