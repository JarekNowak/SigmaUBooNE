#include "HypInput.h"
#include <fstream>
void ctcalib(){
  gInterpreter->GenerateDictionary("vector<vector<int> >","vector");
  TChain ch("ana/OutputTree");
  std::ifstream in("lists/sig_genie.list"); std::string l; int k=0;
  while(std::getline(in,l)&&k<25) if(l.size()){ch.Add(l.c_str());k++;}
  HypInput H; if(!H.attach(&ch)){printf("no attach\n");return;}
  std::vector<std::vector<float>> *w[3]={0,0,0};
  std::vector<int>* vw=0;
  ch.SetBranchAddress("ct_test_window_plane0",&w[0]);
  ch.SetBranchAddress("ct_test_window_plane1",&w[1]);
  ch.SetBranchAddress("ct_test_window_plane2",&w[2]);
  ch.SetBranchAddress("ct_test_primary_vtx_wires",&vw);
  // wire = (A*y + B*z + C)/0.3
  const double A[3]={-0.86603,0.86603,0.0}, B[3]={0.5,0.5,1.0}, C[3]={101.475,99.748,0.0};
  std::vector<double> xs, ts;
  long nchk=0, nok=0;
  for(Long64_t i=0;i<ch.GetEntries();++i){ ch.GetEntry(i); H.convert();
    if(!vw||vw->size()<3||!w[2]||w[2]->empty()) continue;
    // (a) verify the wire transform against the stored vertex wire
    for(int p=0;p<3;++p){
      double pred=(A[p]*H.vtx_y + B[p]*H.vtx_z + C[p])/0.3;
      nchk++; if(fabs(pred - vw->at(p))<2.0) nok++;
    }
    // (b) vertex tick: median time of hits within 3 wires of the vertex wire, plane 2
    int nw=w[2]->size(); int centre=nw/2;
    std::vector<float> tt;
    for(int d=-3;d<=3;++d){ int idx=centre+d; if(idx<0||idx>=nw) continue;
      for(float t:(*w[2])[idx]) tt.push_back(t); }
    if(tt.size()<3) continue;
    std::sort(tt.begin(),tt.end());
    xs.push_back(H.vtx_x); ts.push_back(tt[tt.size()/2]);
  }
  printf("\n  wire transform: %ld/%ld vertex wires reproduced within 2 wires (%.1f%%)\n",
         nok,nchk,100.*nok/nchk);
  // linear fit tick = m*x + c
  int n=xs.size(); double sx=0,st=0,sxx=0,sxt=0;
  for(int i=0;i<n;++i){ sx+=xs[i]; st+=ts[i]; sxx+=xs[i]*xs[i]; sxt+=xs[i]*ts[i]; }
  double m=(n*sxt-sx*st)/(n*sxx-sx*sx), c=(st-m*sx)/n;
  double rms=0; for(int i=0;i<n;++i){ double d=ts[i]-(m*xs[i]+c); rms+=d*d; }
  rms=sqrt(rms/n);
  printf("  tick(x) fit on %d events: tick = %.4f*x + %.1f   rms %.1f ticks (%.2f cm)\n",
         n,m,c,rms,rms*0.3/ (m*0.3) *0.3);
  printf("  implied drift: %.4f cm/tick\n", 1.0/m);
}
