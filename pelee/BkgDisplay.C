// ============================================================================
//  BkgDisplay.C -- candidate-level event displays for the BACKGROUND events that
//  survive the Sigma0 preselection, drawn from the hypana (HyperonProduction)
//  overlay sample.
//
//  WHY this and not EventDisplay.C: the PeLEE truth record stops at the primary
//  vertex, so it cannot say what each candidate really is. hypana backtracks every
//  PFP, so each drawn candidate can be labelled with its TRUE particle. That turns
//  the display into an answer to the question the cut flow raises -- why does a
//  background event pass a Sigma0 selection? -- rather than just a picture.
//
//  Each panel shows the muon/proton/pion track sticks and the photon shower, and
//  is annotated with the true PDG of every candidate plus the reconstructed
//  Lambda and Sigma0 masses. Candidates are assigned with PidAssign.h, i.e. the
//  SAME rule the selection uses, so the picture matches the analysis.
//
//  By default it prefers the events that actually do the damage: background whose
//  reconstructed Lambda mass lands inside the signal window (fakeOnly=true). Those
//  are the events the BDT cannot reject. Set fakeOnly=false to sample all comers.
//
//  Run:  root -l -b -q 'BkgDisplay.C'                       // 2 per mode, fakes
//        root -l -b -q 'BkgDisplay.C("lists/bkg_overlay.list",2,false)'
//  Out:  figs/eventdisplays/hypbkg_<mode>_run<r>_evt<e>.png
// ============================================================================
#include "HypInput.h"
#include "PidAssign.h"
#include "SigmaFV.h"
#include "DisplayUtil.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TChain.h"
#include <fstream>
#include <map>
#include <string>

static const double MUON_TRK_SCORE=0.8, MUON_TRK_DIST=20.0, MUON_TRK_LEN=10.0, MUON_LLR=0.7;
static const int    NSHOWERS_MAX=4, NTRACKS_MAX=5;
static const double TRK_SHR_SPLIT=0.5, TOPO_MIN=0.20;
static const double MP=0.938272, MPI=0.139570;
static const double M_LAMBDA=1.115683, LAM_WIN=0.05;   // "fake Lambda" window

void BkgDisplay(const char* list="lists/bkg_overlay.list",
                int perMode=2, bool fakeOnly=true, int maxFiles=400)
{
  gStyle->SetOptStat(0);
  gSystem->Exec("mkdir -p figs/eventdisplays");

  TChain ch("ana/OutputTree");
  std::ifstream in(list); std::string l; int nf=0;
  while(std::getline(in,l)) if(l.size() && (maxFiles<=0 || nf<maxFiles)){ ch.Add(l.c_str()); nf++; }
  printf("  chained %d files from %s\n", nf, list);

  double topo=0; ch.SetBranchAddress("topological_score",&topo);
  HypInput H;
  if(!H.attach(&ch)){ printf("  not a hypana tree\n"); return; }
  // neutrino-slice requirement, matching the selection (see HypCutFlow.C)
  int fmSliceID=-1; const bool hasSlice = ch.GetBranch("flash_match_nu_slice_ID")!=nullptr;
  if(hasSlice) ch.SetBranchAddress("flash_match_nu_slice_ID",&fmSliceID);

  std::map<std::string,int> got;
  int ndrawn=0;
  const Long64_t N=ch.GetEntries();
  printf("  scanning %lld events (want %d per mode, fakeOnly=%d)\n",N,perMode,(int)fakeOnly); fflush(stdout);

  for(Long64_t ie=0; ie<N; ++ie){
    ch.GetEntry(ie); H.convert();

    // ── background only: drop any true Sigma0 CC anti-numu signal ────────────
    if(H.h_hyperon_pdg==3212 && H.nu_pdg==-14 && H.ccnc==0) continue;

    // ── the selection, up to NTracksCut ──────────────────────────────────────
    if(hasSlice && fmSliceID<0) continue;
    if(!SigmaFV::isFiducial(H.vtx_x,H.vtx_y,H.vtx_z)) continue;
    if(topo<TOPO_MIN) continue;
    const int n=(int)H.trk_score.size();
    if(n<1) continue;
    TVector3 vtx(H.vtx_x,H.vtx_y,H.vtx_z);
    auto isTrk=[&](int i){ return H.trk_score[i]>TRK_SHR_SPLIT; };
    std::vector<TVector3> start(n); std::vector<double> dvtx(n,1e9);
    for(int i=0;i<n;++i){
      if(isTrk(i)) start[i].SetXYZ(H.trk_sx[i],H.trk_sy[i],H.trk_sz[i]);
      else         start[i].SetXYZ(H.shr_sx[i],H.shr_sy[i],H.shr_sz[i]);
      dvtx[i]=(start[i]-vtx).Mag();
    }
    int ntr=0,nsh=0; for(int i=0;i<n;++i)(isTrk(i)?ntr:nsh)++;
    if(!(nsh>0&&nsh<NSHOWERS_MAX)) continue;
    if(!(ntr<NTRACKS_MAX)) continue;

    int mu=-1;
    for(int i=0;i<n;++i)
      if(isTrk(i)&&H.trk_score[i]>MUON_TRK_SCORE&&dvtx[i]<MUON_TRK_DIST&&
         H.trk_len[i]>MUON_TRK_LEN&&H.trk_llr[i]>MUON_LLR)
        if(mu==-1||H.trk_len[i]>H.trk_len[mu]) mu=i;
    if(mu<0) continue;

    std::vector<int> cand;
    for(int i=0;i<n;++i) if(i!=mu&&isTrk(i)) cand.push_back(i);
    if(cand.size()<2) continue;
    std::vector<double> dvv(n), nhv(n);
    for(int i=0;i<n;++i){ dvv[i]=dvtx[i];
      nhv[i]=(i<(int)H.pfnplanehits_Y.size())?H.pfnplanehits_Y[i]:0; }
    int pr=-1, pi=-1;
    if(!PidAssign::assign(cand,H.trk_score,H.trk_llr,H.trk_chipr,H.trk_chipi,
                          H.trk_len,dvv,nhv,pr,pi)) continue;

    int s1=-1; double b1=1e9;
    for(int i=0;i<n;++i){ if(i==mu||i==pr||i==pi||isTrk(i)) continue;
      if(dvtx[i]<b1){ b1=dvtx[i]; s1=i; } }
    if(s1<0) continue;

    // ── reconstructed Lambda / Sigma0 masses (same recipe as the BDT tree) ───
    auto dirOf=[&](int i){ return TVector3(H.trk_dx[i],H.trk_dy[i],H.trk_dz[i]); };
    TVector3 dp=dirOf(pr), di=dirOf(pi);
    double mL=-1, mS=-1;
    const double KEp=H.trk_eproton[pr], ppi=H.trk_rmumom[pi];
    if(KEp>0&&ppi>0){
      const double Ep=KEp+MP, pp=(Ep>MP)?std::sqrt(Ep*Ep-MP*MP):0.0;
      if(pp>0){
        TVector3 Plam=dp*pp+di*ppi;
        const double E=Ep+std::sqrt(ppi*ppi+MPI*MPI);
        const double m2=E*E-Plam.Mag2(); if(m2>0) mL=std::sqrt(m2);
        const double Eg=H.shr_energy_y[s1]/1000.0;
        TVector3 ds(H.shr_px[s1],H.shr_py[s1],H.shr_pz[s1]);
        if(ds.Mag()>0) ds=ds.Unit();
        if(Eg>0&&Eg<10&&mL>0){
          TVector3 Ps=Plam+ds*Eg; const double E2=E+Eg;
          const double m2s=E2*E2-Ps.Mag2(); if(m2s>0) mS=std::sqrt(m2s);
        }
      }
    }
    // keep the events that actually fake a Lambda, when asked
    if(fakeOnly && !(mL>0 && std::fabs(mL-M_LAMBDA)<LAM_WIN)) continue;

    const std::string cat = modeName(H.interaction);
    if(got[cat]>=perMode) continue;
    got[cat]++; ndrawn++;

    // ── draw ─────────────────────────────────────────────────────────────────
    std::vector<Seg> segs;
    auto addTrk=[&](int k,int col){
      TVector3 a=start[k], b=a+dirOf(k)*std::max(2.0,(double)H.trk_len[k]);
      segs.push_back({a,b,col,3,1}); };
    addTrk(mu,kRed+1); addTrk(pr,kAzure+1); addTrk(pi,kGreen+2);
    TVector3 g0=start[s1], gd(H.shr_px[s1],H.shr_py[s1],H.shr_pz[s1]);
    if(gd.Mag()>0) gd=gd.Unit();
    segs.push_back({vtx,g0,kMagenta+1,2,2});          // vertex -> conversion point
    segs.push_back({g0,g0+gd*15.0,kMagenta+1,3,1});   // shower axis

    TCanvas* c=new TCanvas(Form("bev%d",ndrawn),"",1300,600); c->Divide(2,1);
    drawView(c->cd(1),segs,vtx,2,0,"Z [cm]","X [cm]");
    drawView(c->cd(2),segs,vtx,2,1,"Z [cm]","Y [cm]");

    c->cd(1);
    TLegend* lg=new TLegend(0.12,0.70,0.60,0.90);
    lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.032);
    lg->SetHeader("candidate  (true particle)");
    lg->AddEntry((TObject*)0,Form("#mu cand      #rightarrow  %s",pdgName(H.bkt_pdg[mu])),"");
    lg->AddEntry((TObject*)0,Form("p cand       #rightarrow  %s",pdgName(H.bkt_pdg[pr])),"");
    lg->AddEntry((TObject*)0,Form("#pi cand      #rightarrow  %s",pdgName(H.bkt_pdg[pi])),"");
    lg->AddEntry((TObject*)0,Form("#gamma shower  #rightarrow  %s",pdgName(H.bkt_pdg[s1])),"");
    lg->Draw();

    c->cd(2);
    TPaveText* pt=new TPaveText(0.16,0.68,0.72,0.90,"NDC");
    pt->SetBorderSize(0); pt->SetFillStyle(0); pt->SetTextAlign(12); pt->SetTextSize(0.032);
    pt->AddText(Form("BACKGROUND  %s   run %d  evt %d",cat.c_str(),H.run,H.evt));
    pt->AddText(Form("#nu pdg %d   CC/NC %d   topo %.2f",H.nu_pdg,H.ccnc,topo));
    pt->AddText(Form("n_{trk}=%d  n_{shr}=%d",ntr,nsh));
    pt->AddText(Form("m_{#Lambda} = %s GeV   m_{#Sigma^{0}} = %s GeV",
                     mL>0?Form("%.3f",mL):"--", mS>0?Form("%.3f",mS):"--"));
    pt->Draw();

    TString out=Form("figs/eventdisplays/hypbkg_%s_run%d_evt%d.png",cat.c_str(),H.run,H.evt);
    c->SaveAs(out);
    printf("  drew %-6s m_lam=%s  mu->%s p->%s pi->%s gam->%s   %s\n",
           cat.c_str(), mL>0?Form("%.3f",mL):"--",
           pdgName(H.bkt_pdg[mu]),pdgName(H.bkt_pdg[pr]),
           pdgName(H.bkt_pdg[pi]),pdgName(H.bkt_pdg[s1]), out.Data());
    fflush(stdout);
    delete c;
  }
  printf("\n  total displays: %d\n",ndrawn); fflush(stdout);
}
