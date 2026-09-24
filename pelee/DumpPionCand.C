// DumpPionCand.C -- one pass over the hypana signal sample dumping every non-muon
// track in events that survive to the proton/pion assignment stage, with its truth
// pdg and every discriminant either ntuple format offers. Lets the assignment be
// re-optimised offline without re-reading 4851 files each time.
#include "HypInput.h"
#include "SigmaFV.h"
#include "TVector3.h"
#include <fstream>
static const double MUON_TRK_SCORE=0.8,MUON_TRK_DIST=20.0,MUON_TRK_LEN=10.0,MUON_LLR=0.7;
static const int NSHOWERS_MAX=4,NTRACKS_MAX=5;
static const double TRK_SHR_SPLIT=0.5,TOPO_MIN=0.20;

void DumpPionCand(const char* list="lists/sig_genie.list", const char* out="PionCand.root"){
  TChain ch("ana/OutputTree");
  std::ifstream in(list); std::string l;
  while(std::getline(in,l)) if(l.size()) ch.Add(l.c_str());
  double topo=0; ch.SetBranchAddress("topological_score",&topo);
  HypInput H; if(!H.attach(&ch)){printf("  not hypana\n");return;}
  // NOTE: these must be bound AFTER HypInput::attach. Binding them first leaves the
  // vectors empty on a TChain -- attach() re-binds the branch list and drops them.
  std::vector<double> *rms=nullptr,*d3=nullptr,*dedx2=nullptr;
  ch.SetBranchAddress("trk_rms_deflection",&rms);
  ch.SetBranchAddress("trk_three_plane_dedx",&d3);
  if(ch.GetBranch("trk_mean_dedx_plane2")) ch.SetBranchAddress("trk_mean_dedx_plane2",&dedx2);

  TFile fo(out,"RECREATE");
  TTree* T=new TTree("cand","non-muon track candidates in surviving signal events");
  Int_t   ev=0, pdg=0, ntrk=0, nhy=0;
  Float_t score=0,llr=0,chipr=0,chipi=0,len=0,dvtx=0,wig=0,dedx3=0,dedxY=0;
  T->Branch("ev",&ev);        T->Branch("pdg",&pdg);      T->Branch("ntrk",&ntrk);
  T->Branch("score",&score);  T->Branch("llr",&llr);      T->Branch("chipr",&chipr);
  T->Branch("chipi",&chipi);  T->Branch("len",&len);      T->Branch("dvtx",&dvtx);
  T->Branch("wig",&wig);      T->Branch("dedx3",&dedx3);  T->Branch("dedxY",&dedxY);
  T->Branch("nhy",&nhy);

  long nev=0;
  Long64_t N=ch.GetEntries();
  for(Long64_t ie=0;ie<N;++ie){
    ch.GetEntry(ie); H.convert();
    if(!(H.h_hyperon_pdg==3212 && H.nu_pdg==-14 && H.ccnc==0)) continue;   // signal only
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
    int mu=-1; for(int i=0;i<n;++i)
      if(isTrk(i)&&H.trk_score[i]>MUON_TRK_SCORE&&dv[i]<MUON_TRK_DIST&&
         H.trk_len[i]>MUON_TRK_LEN&&H.trk_llr[i]>MUON_LLR)
        if(mu==-1||H.trk_len[i]>H.trk_len[mu]) mu=i;
    if(mu<0) continue;
    if(!(ntr<NTRACKS_MAX)) continue;          // same preselection as the BDT tree
    // need at least two non-muon tracks for a p/pi pair to exist
    int navail=0; for(int i=0;i<n;++i) if(i!=mu&&isTrk(i)) navail++;
    if(navail<2) continue;
    ntrk=navail;
    for(int i=0;i<n;++i){
      if(i==mu||!isTrk(i)) continue;
      ev=(int)nev; pdg=H.bkt_pdg[i];
      score=H.trk_score[i]; llr=H.trk_llr[i]; chipr=H.trk_chipr[i]; chipi=H.trk_chipi[i];
      len=H.trk_len[i]; dvtx=dv[i]; nhy=H.pfnplanehits_Y.size()>(size_t)i?H.pfnplanehits_Y[i]:0;
      wig   = (rms  &&i<(int)rms->size()  &&std::isfinite(rms->at(i)))   ? rms->at(i)   : -1;
      dedx3 = (d3   &&i<(int)d3->size()   &&std::isfinite(d3->at(i)))    ? d3->at(i)    : -1;
      dedxY = (dedx2&&i<(int)dedx2->size()&&std::isfinite(dedx2->at(i))) ? dedx2->at(i) : -1;
      T->Fill();
    }
    nev++;
  }
  fo.cd(); T->Write();
  printf("\n  dumped %ld signal events, %lld candidate tracks -> %s\n",nev,T->GetEntries(),out);
  fo.Close();
}
