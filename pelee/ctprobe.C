#include "HypInput.h"
#include <fstream>
void ctprobe(){
  TChain ch("ana/OutputTree");
  std::ifstream in("lists/sig_genie.list"); std::string l; int k=0;
  while(std::getline(in,l)&&k<3) if(l.size()){ch.Add(l.c_str());k++;}
  HypInput H; if(!H.attach(&ch)){printf("no attach\n");return;}
  std::vector<std::vector<float>> *w[3]={0,0,0};
  std::vector<std::vector<float>> *adc[3]={0,0,0};
  std::vector<std::vector<int>>   *pdg[3]={0,0,0};
  std::vector<int>* vw=0;
  ch.SetBranchAddress("ct_test_window_plane0",&w[0]);
  ch.SetBranchAddress("ct_test_window_plane1",&w[1]);
  ch.SetBranchAddress("ct_test_window_plane2",&w[2]);
  ch.SetBranchAddress("ct_test_adc_plane0",&adc[0]);
  ch.SetBranchAddress("ct_test_adc_plane1",&adc[1]);
  ch.SetBranchAddress("ct_test_adc_plane2",&adc[2]);
  ch.SetBranchAddress("ct_test_pdg_plane0",&pdg[0]);
  ch.SetBranchAddress("ct_test_pdg_plane1",&pdg[1]);
  ch.SetBranchAddress("ct_test_pdg_plane2",&pdg[2]);
  ch.SetBranchAddress("ct_test_primary_vtx_wires",&vw);
  int shown=0;
  for(Long64_t i=0;i<ch.GetEntries()&&shown<4;++i){ ch.GetEntry(i); H.convert();
    if(!w[2]||w[2]->empty()) continue;
    printf("\n evt %lld  vtx=(%.1f,%.1f,%.1f)  vtxWires=",i,H.vtx_x,H.vtx_y,H.vtx_z);
    if(vw) for(auto x:*vw) printf(" %d",x);
    printf("\n");
    for(int p=0;p<3;++p){
      if(!w[p]) continue;
      size_t nw=w[p]->size(); long nhit=0; float tmin=1e9,tmax=-1e9;
      for(auto&col:*w[p]){ nhit+=col.size(); for(float t:col){ if(t<tmin)tmin=t; if(t>tmax)tmax=t; } }
      printf("   plane%d: %zu wires in window, %ld hits, t in [%.0f,%.0f]",p,nw,nhit,tmin,tmax);
      if(adc[p]&&adc[p]->size()==nw) printf("  adc ok");
      if(pdg[p]&&pdg[p]->size()==nw) printf("  pdg ok");
      printf("\n"); }
    shown++; }
}
