#include "HypInput.h"
#include "SigmaFV.h"
#include <fstream>
void hchk3(){
  TChain ch("ana/OutputTree");
  std::ifstream in("/tmp/claude-661400007/-home-t2k-nowak-MicroBooNE-Sigma/4861d86a-5f29-4f62-92f1-f0c7167b5758/scratchpad/b40.list");
  std::string l; while(std::getline(in,l)) if(l.size()) ch.Add(l.c_str());
  HypInput H; if(!H.attach(&ch)){printf("no attach\n");return;}
  int nsl=-1, pnsl=-99, fmsl=-99; double topo=-99;
  ch.SetBranchAddress("n_slices",&nsl);
  ch.SetBranchAddress("pandora_nu_slice_ID",&pnsl);
  ch.SetBranchAddress("flash_match_nu_slice_ID",&fmsl);
  ch.SetBranchAddress("topological_score",&topo);
  Long64_t N=ch.GetEntries();
  long z_nsl=0,z_pn=0,z_fm=0,both=0, fv_pn=0, fv=0;
  for(Long64_t i=0;i<N;++i){ ch.GetEntry(i); H.convert();
    if(nsl<=0) z_nsl++;
    if(pnsl<0) z_pn++;
    if(fmsl<0) z_fm++;
    bool infv=SigmaFV::isFiducial(H.vtx_x,H.vtx_y,H.vtx_z);
    if(infv){ fv++; if(pnsl>=0) fv_pn++; } }
  printf("  hypana overlay, %lld events\n",N);
  printf("    n_slices <= 0                 : %ld  (%.1f%%)\n",z_nsl,100.*z_nsl/N);
  printf("    pandora_nu_slice_ID < 0       : %ld  (%.1f%%)\n",z_pn,100.*z_pn/N);
  printf("    flash_match_nu_slice_ID < 0   : %ld  (%.1f%%)\n",z_fm,100.*z_fm/N);
  printf("    pass FV                       : %ld  (%.1f%%)\n",fv,100.*fv/N);
  printf("    pass FV AND pandora slice >=0 : %ld  (%.1f%% of all, %.1f%% of FV)\n",
         fv_pn,100.*fv_pn/N,100.*fv_pn/fv);
}
