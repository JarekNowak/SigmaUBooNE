#include "HypInput.h"
#include "SigmaFV.h"
#include <fstream>
void hchk(){
  TChain ch("ana/OutputTree");
  std::ifstream in("/tmp/claude-661400007/-home-t2k-nowak-MicroBooNE-Sigma/4861d86a-5f29-4f62-92f1-f0c7167b5758/scratchpad/b40.list");
  std::string l; while(std::getline(in,l)) if(l.size()) ch.Add(l.c_str());
  HypInput H; if(!H.attach(&ch)){printf("no attach\n");return;}
  Long64_t N=ch.GetEntries(); long nopfp=0, fv=0, zero=0;
  for(Long64_t i=0;i<N;++i){ ch.GetEntry(i); H.convert();
    if(H.trk_score.size()==0) nopfp++;
    if(H.vtx_x==0&&H.vtx_y==0&&H.vtx_z==0) zero++;
    if(SigmaFV::isFiducial(H.vtx_x,H.vtx_y,H.vtx_z)) fv++; }
  printf("  hypana overlay, %lld events in 40 files\n",N);
  printf("    zero PFPs                  : %ld  (%.1f%%)\n",nopfp,100.*nopfp/N);
  printf("    reco vtx exactly 0,0,0     : %ld  (%.1f%%)\n",zero,100.*zero/N);
  printf("    pass FiducialV             : %ld  (%.1f%%)\n",fv,100.*fv/N);
}
