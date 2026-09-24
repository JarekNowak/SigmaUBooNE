// pion_purity_scan.C -- the pion candidate is only 34.2% correct and its largest
// contaminant is PHOTONS (25%), i.e. shower-like PFPs passing the 0.5 track/shower
// split. Scan a dedicated track-score requirement for the pion candidate, and test
// the discriminants hypana provides that PeLEE did not (wiggliness, 3-plane dE/dx).
#include "HypInput.h"
#include "SigmaFV.h"
#include "TVector3.h"
#include <fstream>
static const double MUON_TRK_SCORE=0.8, MUON_TRK_DIST=20.0, MUON_TRK_LEN=10.0, MUON_LLR=0.7;
static const int NSHOWERS_MAX=4, NTRACKS_MAX=5;
static const double TRK_SHR_SPLIT=0.5, TOPO_MIN=0.20, LLR_SCALE=0.591, DCHI_SCALE=124.5;

void pion_purity_scan(const char* list="/home/t2k/nowak/.claude/jobs/a2cc9f4f/tmp/nuwro_hypana.list"){
  TChain ch("ana/OutputTree");
  std::ifstream in(list); std::string l; while(std::getline(in,l)) if(l.size()) ch.Add(l.c_str());
  double topo=0; ch.SetBranchAddress("topological_score",&topo);
  HypInput H; if(!H.attach(&ch)){printf(" not hypana\n");return;}

  const std::vector<double> SCAN={0.50,0.60,0.70,0.80,0.85,0.90,0.95};
  const int NS=SCAN.size();
  std::vector<long> sig(NS,0), bkg(NS,0), pok(NS,0);
  long sig0=0,bkg0=0,pok0=0;   // baseline (0.5)
  // wiggliness of true pi- vs true gamma among pion candidates
  TH1D* wPi=new TH1D("wPi",";rms deflection",60,0,1.0);
  TH1D* wGa=new TH1D("wGa",";rms deflection",60,0,1.0);

  // bound after attach() above: binding before it leaves the vector empty on a TChain
  std::vector<double>* rmsdef=nullptr;
  ch.SetBranchAddress("trk_rms_deflection",&rmsdef);

  for(Long64_t ie=0;ie<ch.GetEntries();++ie){
    ch.GetEntry(ie); H.convert();
    const bool isSig=(H.h_hyperon_pdg==3212 && H.nu_pdg==-14 && H.ccnc==0);
    if(!SigmaFV::isFiducial(H.vtx_x,H.vtx_y,H.vtx_z)) continue;
    if(topo<TOPO_MIN) continue;
    const int n=(int)H.trk_score.size();
    std::vector<TVector3> st(n); std::vector<double> dv(n,1e9);
    TVector3 vtx(H.vtx_x,H.vtx_y,H.vtx_z);
    auto isTrk=[&](int i){ return H.trk_score[i]>TRK_SHR_SPLIT; };
    for(int i=0;i<n;++i){ if(isTrk(i)) st[i].SetXYZ(H.trk_sx[i],H.trk_sy[i],H.trk_sz[i]);
      else st[i].SetXYZ(H.shr_sx[i],H.shr_sy[i],H.shr_sz[i]); dv[i]=(st[i]-vtx).Mag(); }
    int ntr=0,nsh=0; for(int i=0;i<n;++i)(isTrk(i)?ntr:nsh)++;
    if(!(nsh>0&&nsh<NSHOWERS_MAX)) continue;
    if(!(ntr<NTRACKS_MAX)) continue;
    int mu=-1; for(int i=0;i<n;++i)
      if(isTrk(i)&&H.trk_score[i]>MUON_TRK_SCORE&&dv[i]<MUON_TRK_DIST&&H.trk_len[i]>MUON_TRK_LEN&&H.trk_llr[i]>MUON_LLR)
        if(mu==-1||H.trk_len[i]>H.trk_len[mu]) mu=i;
    if(mu<0) continue;
    auto D=[&](int i){ double d=H.trk_llr[i]/LLR_SCALE; double dc=H.trk_chipr[i]-H.trk_chipi[i];
      if(std::isfinite(dc)&&H.trk_chipr[i]>-900&&H.trk_chipi[i]>-900) d+=dc/DCHI_SCALE; return d; };
    int pr=-1; for(int i=0;i<n;++i){ if(i==mu||!isTrk(i))continue; if(pr==-1||D(i)<D(pr))pr=i; }
    if(pr<0) continue;
    // scan: pion candidate restricted to PFPs with trk_score above threshold
    for(int s=0;s<NS;++s){
      int pi=-1;
      for(int i=0;i<n;++i){ if(i==mu||i==pr)continue;
        if(H.trk_score[i]<=SCAN[s]) continue;
        if(pi==-1||D(i)>D(pi)) pi=i; }
      if(pi<0) continue;
      if(isSig){ sig[s]++; if(H.bkt_pdg[pi]==-211) pok[s]++; } else bkg[s]++;
      if(s==0&&isSig&&rmsdef&&pi<(int)rmsdef->size()){
        if(H.bkt_pdg[pi]==-211) wPi->Fill(rmsdef->at(pi));
        if(H.bkt_pdg[pi]==22)   wGa->Fill(rmsdef->at(pi));
      }
    }
  }
  printf("\n  === pion-candidate track-score scan ===\n");
  printf("  (baseline is the current 0.50 split)\n\n");
  printf("  %8s %8s %9s %10s %10s %9s\n","minScore","signal","bkg","pi- purity","sig rel","FoM rel");
  double base=sig[0]/std::sqrt((double)sig[0]+bkg[0]);
  for(int s=0;s<NS;++s){
    if(sig[s]+bkg[s]<=0) continue;
    double fom=sig[s]/std::sqrt((double)sig[s]+bkg[s]);
    printf("  %8.2f %8ld %9ld %9.1f%% %10.3f %9.3f\n",SCAN[s],sig[s],bkg[s],
      sig[s]?100.*pok[s]/sig[s]:0., (double)sig[s]/sig[0], fom/base);
  }
  printf("\n  === wiggliness (trk_rms_deflection) of pion candidates, signal ===\n");
  printf("  true pi- : n=%.0f mean %.4f\n", wPi->Integral(), wPi->GetMean());
  printf("  true gam : n=%.0f mean %.4f\n", wGa->Integral(), wGa->GetMean());
}
