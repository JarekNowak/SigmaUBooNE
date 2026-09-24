// RangeFit.C -- fit p(range) for true protons and true pions in our own sample,
// and compare against the estimators the selection currently uses. The pion
// momentum is presently taken from trk_range_muon_mom, i.e. the MUON hypothesis,
// which is a known bias with no pion-hypothesis branch available.
#include "HypInput.h"
#include "SigmaFV.h"
#include <fstream>
void RangeFit(int maxFiles=800){
  TChain ch("ana/OutputTree");
  std::ifstream in("lists/sig_genie.list"); std::string l; int k=0;
  while(std::getline(in,l)) if(l.size()&&(maxFiles<=0||k<maxFiles)){ch.Add(l.c_str());k++;}
  HypInput H; if(!H.attach(&ch)){printf("no attach\n");return;}
  std::vector<double> *tpx=0,*tpy=0,*tpz=0,*tex=0,*tey=0,*tez=0;
  std::vector<int>* tpdg=0;
  auto B=[&](const char*n,void*p){ if(ch.GetBranch(n)) ch.SetBranchAddress(n,p); };
  B("pfp_true_px",&tpx); B("pfp_true_py",&tpy); B("pfp_true_pz",&tpz);
  B("pfp_true_pdg",&tpdg);
  B("trk_end_x",&tex); B("trk_end_y",&tey); B("trk_end_z",&tez);
  FILE* fo=fopen("RangeFit.csv","w"); fprintf(fo,"pdg,range,ptrue,pcur\n");
  long np=0,npi=0;
  for(Long64_t ie=0; ie<ch.GetEntries(); ++ie){ ch.GetEntry(ie); H.convert();
    if(!(H.h_hyperon_pdg==3212&&H.nu_pdg==-14&&H.ccnc==0)) continue;
    const int n=(int)H.trk_score.size();
    for(int i=0;i<n;++i){
      if(H.trk_score[i]<=0.5) continue;                       // tracks only
      if(!tpdg||i>=(int)tpdg->size()) continue;
      const int pdg=tpdg->at(i);
      if(pdg!=2212 && pdg!=-211) continue;
      // require a contained end point: range-momentum only holds for stoppers
      if(tex&&i<(int)tex->size() &&
         !SigmaFV::isFiducial(tex->at(i),tey->at(i),tez->at(i))) continue;
      const double R=H.trk_len[i]; if(!(R>0.3&&R<300)) continue;
      if(!tpx||i>=(int)tpx->size()) continue;
      const double pt=std::sqrt(tpx->at(i)*tpx->at(i)+tpy->at(i)*tpy->at(i)+tpz->at(i)*tpz->at(i));
      if(!(pt>0&&pt<5)) continue;
      double pcur=-1;
      if(pdg==2212){ const double KE=H.trk_eproton[i];           // current: proton KE
        if(KE>0){ const double E=KE+0.938272; pcur=std::sqrt(E*E-0.938272*0.938272); } np++; }
      else { pcur=H.trk_rmumom[i]; npi++; }                      // current: MUON range mom
      fprintf(fo,"%d,%.4f,%.5f,%.5f\n",pdg,R,pt,pcur);
    }
  }
  fclose(fo);
  printf("  wrote RangeFit.csv : %ld true protons, %ld true pions\n",np,npi);
}
