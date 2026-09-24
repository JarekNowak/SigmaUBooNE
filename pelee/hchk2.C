#include "HypInput.h"
#include "SigmaFV.h"
#include <fstream>
void hchk2(){
  TChain ch("ana/OutputTree");
  std::ifstream in("/tmp/claude-661400007/-home-t2k-nowak-MicroBooNE-Sigma/4861d86a-5f29-4f62-92f1-f0c7167b5758/scratchpad/b40.list");
  std::string l; while(std::getline(in,l)) if(l.size()) ch.Add(l.c_str());
  HypInput H; if(!H.attach(&ch)){printf("no attach\n");return;}
  double topo=-99; ch.SetBranchAddress("topological_score",&topo);
  Long64_t N=ch.GetEntries();
  long fv=0, fv_noscore=0, fv_lowscore=0, fv_pass=0, nosc=0;
  for(Long64_t i=0;i<N;++i){ ch.GetEntry(i); H.convert();
    bool infv=SigmaFV::isFiducial(H.vtx_x,H.vtx_y,H.vtx_z);
    if(topo<=0) nosc++;
    if(!infv) continue; fv++;
    if(topo<=0)   fv_noscore++;
    if(topo<0.20) fv_lowscore++;
    else          fv_pass++; }
  printf("  hypana overlay, %lld events\n",N);
  printf("    topological_score <= 0 (no nu slice), all events : %ld  (%.1f%%)\n",nosc,100.*nosc/N);
  printf("    pass FiducialV                                    : %ld  (%.1f%%)\n",fv,100.*fv/N);
  printf("      of those, topo <= 0                             : %ld  (%.1f%% of FV)\n",fv_noscore,100.*fv_noscore/fv);
  printf("      of those, topo < 0.20 (fail CosmicVeto)         : %ld  (%.1f%% of FV)\n",fv_lowscore,100.*fv_lowscore/fv);
  printf("      surviving CosmicVeto                            : %ld  (%.1f%% of FV)\n",fv_pass,100.*fv_pass/fv);
}
