#include "HypInput.h"
#include <fstream>
void lamflight(){
  TChain ch("ana/OutputTree");
  std::ifstream in("lists/sig_genie.list"); std::string l; int k=0;
  while(std::getline(in,l)&&k<300) if(l.size()){ch.Add(l.c_str());k++;}
  HypInput H; if(!H.attach(&ch)){printf("no attach\n");return;}
  double lsx=0,lsy=0,lsz=0,lex=0,ley=0,lez=0,lmom=0;
  auto B=[&](const char*n,void*p){ if(ch.GetBranch(n)) ch.SetBranchAddress(n,p); };
  B("mc_sigmazero_lambda_start_x",&lsx); B("mc_sigmazero_lambda_start_y",&lsy);
  B("mc_sigmazero_lambda_start_z",&lsz); B("mc_sigmazero_lambda_end_x",&lex);
  B("mc_sigmazero_lambda_end_y",&ley);   B("mc_sigmazero_lambda_end_z",&lez);
  B("mc_sigmazero_lambda_mom",&lmom);
  std::vector<double> d;
  for(Long64_t i=0;i<ch.GetEntries();++i){ ch.GetEntry(i); H.convert();
    if(!(H.h_hyperon_pdg==3212&&H.nu_pdg==-14&&H.ccnc==0)) continue;
    double dx=lex-lsx, dy=ley-lsy, dz=lez-lsz;
    double r=sqrt(dx*dx+dy*dy+dz*dz);
    if(r>0&&r<1000) d.push_back(r); }
  if(d.empty()){ printf("  no flight distances\n"); return; }
  std::sort(d.begin(),d.end());
  double mean=0; for(double x:d) mean+=x; mean/=d.size();
  auto q=[&](double f){ return d[(size_t)(f*(d.size()-1))]; };
  printf("\n  true Lambda flight distance (Sigma0 sample), %zu events\n",d.size());
  printf("    mean %.2f cm   median %.2f cm\n",mean,q(0.5));
  printf("    quartiles: 25%% %.2f   75%% %.2f   90%% %.2f cm\n",q(0.25),q(0.75),q(0.90));
  const double cuts[]={1,2,3,5,10};
  for(double c:cuts){ size_t n=0; for(double x:d) if(x>c) n++;
    printf("    flight > %4.1f cm : %5.1f%%\n",c,100.*n/d.size()); }
}
