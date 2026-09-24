// Pure truth: where does the Sigma0 photon convert relative to the Lambda flight?
// If it converts inside the gap between vertex and Lambda decay, it can bridge them.
#include "HypInput.h"
#include <fstream>
void PhotonGeom(int maxFiles=300){
  TChain ch("ana/OutputTree");
  std::ifstream in("lists/sig_genie.list"); std::string l; int k=0;
  while(std::getline(in,l)) if(l.size()&&(maxFiles<=0||k<maxFiles)){ch.Add(l.c_str());k++;}
  HypInput H; if(!H.attach(&ch)){printf("no attach\n");return;}
  double gsx=0,gsy=0,gsz=0,gex=0,gey=0,gez=0,lsx=0,lsy=0,lsz=0,lex=0,ley=0,lez=0,gmom=0;
  auto B=[&](const char*n,void*q){ if(ch.GetBranch(n)) ch.SetBranchAddress(n,q); };
  B("mc_sigmazero_photon_start_x",&gsx);B("mc_sigmazero_photon_start_y",&gsy);B("mc_sigmazero_photon_start_z",&gsz);
  B("mc_sigmazero_photon_end_x",&gex);  B("mc_sigmazero_photon_end_y",&gey);  B("mc_sigmazero_photon_end_z",&gez);
  B("mc_sigmazero_lambda_start_x",&lsx);B("mc_sigmazero_lambda_start_y",&lsy);B("mc_sigmazero_lambda_start_z",&lsz);
  B("mc_sigmazero_lambda_end_x",&lex);  B("mc_sigmazero_lambda_end_y",&ley);  B("mc_sigmazero_lambda_end_z",&lez);
  B("mc_sigmazero_photon_mom",&gmom);
  std::vector<double> conv, cosang, pmom; long nInGap=0, nLong=0, nLongGap=0, N=0;
  for(Long64_t ie=0;ie<ch.GetEntries();++ie){ ch.GetEntry(ie); H.convert();
    if(!(H.h_hyperon_pdg==3212&&H.nu_pdg==-14&&H.ccnc==0)) continue;
    TVector3 g0(gsx,gsy,gsz), g1(gex,gey,gez), l0(lsx,lsy,lsz), l1(lex,ley,lez);
    const double dg=(g1-g0).Mag(), dl=(l1-l0).Mag();
    if(!(dg>0&&dg<2000&&dl>0&&dl<1000)) continue;
    N++; conv.push_back(dg); pmom.push_back(gmom);
    TVector3 ug=(g1-g0).Unit(), ul=(l1-l0).Unit(); cosang.push_back(ug.Dot(ul));
    if(dg<dl) nInGap++;                       // converts before the Lambda decays
    if(dl>5){ nLong++; if(dg<dl) nLongGap++; }
  }
  auto med=[&](std::vector<double> v){ std::sort(v.begin(),v.end()); return v[v.size()/2]; };
  printf("\n  Sigma0 photon geometry, %ld signal events (truth)\n",N);
  printf("    photon conversion distance : median %.1f cm   (LAr conversion length ~18 cm)\n",med(conv));
  printf("    photon momentum            : median %.0f MeV/c\n",1000*med(pmom));
  printf("    cos(photon, Lambda)        : median %.2f  (both boosted forward)\n",med(cosang));
  printf("    converts BEFORE Lambda decays        : %.1f%%\n",100.*nInGap/N);
  printf("    among Lambda flight > 5 cm (n=%ld)   : %.1f%% convert inside the gap\n",nLong,nLong?100.*nLongGap/nLong:0.);
}
