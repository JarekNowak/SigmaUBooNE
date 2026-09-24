// Lambda mass on TRUTH-CORRECT (p,pi-) pairs in signal, three ways:
//   current : proton p from trk_energy_proton treated as KE ; pion p = muon range mom
//   fixA    : proton p = trk_energy_proton used directly as momentum ; pion as current
//   fixB    : fixA plus pion p from the fitted range relation 82.438*R^0.2593 MeV/c
#include "HypInput.h"
#include <fstream>
void MassFixC(int maxFiles=800){
  TChain ch("ana/OutputTree");
  std::ifstream in("lists/sig_genie.list"); std::string l; int k=0;
  while(std::getline(in,l)) if(l.size()&&(maxFiles<=0||k<maxFiles)){ch.Add(l.c_str());k++;}
  HypInput H; if(!H.attach(&ch)){printf("no attach\n");return;}
  std::vector<int>* tpdg=0; if(ch.GetBranch("pfp_true_pdg")) ch.SetBranchAddress("pfp_true_pdg",&tpdg);
  const double MP=0.938272, MPI=0.139570;
  std::vector<double> m0,mA,mB,mC;
  for(Long64_t ie=0;ie<ch.GetEntries();++ie){ ch.GetEntry(ie); H.convert();
    if(!(H.h_hyperon_pdg==3212&&H.nu_pdg==-14&&H.ccnc==0)) continue;
    if(!tpdg) continue;
    const int n=(int)H.trk_score.size(); int pr=-1,pi=-1;
    for(int i=0;i<n&&i<(int)tpdg->size();++i){ if(H.trk_score[i]<=0.5) continue;
      if(tpdg->at(i)==2212&&pr<0) pr=i; if(tpdg->at(i)==-211&&pi<0) pi=i; }
    if(pr<0||pi<0) continue;
    TVector3 dp(H.trk_dx[pr],H.trk_dy[pr],H.trk_dz[pr]), di(H.trk_dx[pi],H.trk_dy[pi],H.trk_dz[pi]);
    if(dp.Mag()<=0||di.Mag()<=0) continue; dp=dp.Unit(); di=di.Unit();
    const double Ep_branch=H.trk_eproton[pr], ppi_mu=H.trk_rmumom[pi], Rpi=H.trk_len[pi];
    if(!(Ep_branch>0&&ppi_mu>0&&Rpi>0)) continue;
    auto mass=[&](double pp,double ppi){
      const double Ep=std::sqrt(pp*pp+MP*MP), Epi=std::sqrt(ppi*ppi+MPI*MPI);
      TVector3 P=dp*pp+di*ppi; const double m2=(Ep+Epi)*(Ep+Epi)-P.Mag2();
      return m2>0?std::sqrt(m2):-1.0; };
    // current: treat branch as KE
    const double E0=Ep_branch+MP, pp0=std::sqrt(E0*E0-MP*MP);
    const double ppiB=(82.438*std::pow(Rpi,0.2593))/1000.0;
    // fixC: pion momentum = 1.2494 * muon-hypothesis range momentum. The factor is
    // median(p_true/p_mu) on 277 true pi- (RangeFit.csv): zero residual bias, 7.8%
    // half-IQR, vs 8.2% for the pure power-law fit -- keeps the tight muon spread.
    double a=mass(pp0,ppi_mu), b=mass(Ep_branch,ppi_mu), c=mass(Ep_branch,ppiB), d=mass(Ep_branch,1.2494*ppi_mu);
    if(a>0) m0.push_back(a); if(b>0) mA.push_back(b); if(c>0) mB.push_back(c); if(d>0) mC.push_back(d);
  }
  auto rep=[&](const char* nm,std::vector<double>& v){
    std::sort(v.begin(),v.end()); if(v.size()<10){printf("  %-8s too few\n",nm);return;}
    long in20=0,in50=0; for(double x:v){ if(fabs(x-1.1157)<0.02) in20++; if(fabs(x-1.1157)<0.05) in50++; }
    printf("  %-8s n=%4zu  median %.4f  IQR/2 %.4f  within 20 MeV %5.1f%%  within 50 MeV %5.1f%%\n",
      nm,v.size(),v[v.size()/2],(v[3*v.size()/4]-v[v.size()/4])/2,100.*in20/v.size(),100.*in50/v.size()); };
  printf("\n  Lambda mass on truth-correct p/pi- pairs (true Lambda = 1.1157 GeV)\n");
  rep("current",m0); rep("fixA",mA); rep("fixB",mB); rep("fixC",mC);
}
