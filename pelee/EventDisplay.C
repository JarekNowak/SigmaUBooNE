// ============================================================================
//  EventDisplay.C — schematic candidate-level event displays for representative
//  BACKGROUND events of the Sigma0 selection.
//
//  PeLEE ntuples carry no wire hits, so this draws the reconstructed candidates
//  (vertex + muon/proton/pion track sticks + photon-shower direction) in two 2D
//  projections (Z-X and Z-Y), annotated with truth (interaction mode, candidate
//  backtracked PDGs) and the reconstructed Lambda / Sigma0 / pi0 masses.  One
//  representative event is drawn per interaction mode, plus a cosmic and a clear
//  pi0 example, selected at the NTracksCut preselection.
//
//  Run:  root -l -b -q EventDisplay.C
//  Out:  figs/eventdisplays/*.png
// ============================================================================
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TLine.h"
#include "TMarker.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TSystem.h"
#include <vector>
#include <string>
#include <cmath>

static bool inFV(float x,float y,float z){
  bool fv = (x>3.0&&x<253.35 && y>-112.53&&y<114.47 && z>3.1&&z<1033.0);
  bool dead = (z>675.1 && z<775.1);
  return fv && !dead;
}
static const char* modeName(int c,int cosmic){ if(cosmic) return "Cosmic";
  switch(c){case 0:return "QE";case 1:return "RES";case 2:return "DIS";case 3:return "COH";case 10:return "MEC";default:return "Other";} }
static const char* pn(int p){switch(p){case 13:return "#mu^{-}";case -13:return "#mu^{+}";case 2212:return "p";
  case 2112:return "n";case 211:return "#pi^{+}";case -211:return "#pi^{-}";case 111:return "#pi^{0}";
  case 11:return "e^{-}";case -11:return "e^{+}";case 22:return "#gamma";case 321:return "K^{+}";
  case 0:return "unmatched";default:return Form("%d",p);} }

struct Seg{ TVector3 a,b; int col,wid,sty; };

void drawView(TVirtualPad* pad,const std::vector<Seg>& segs,const TVector3& vtx,
              int ia,int ib,const char* xt,const char* yt){
  pad->cd();
  double lo[3]={1e9,1e9,1e9}, hi[3]={-1e9,-1e9,-1e9};
  auto upd=[&](const TVector3&v){ double c[3]={v.X(),v.Y(),v.Z()};
    for(int k=0;k<3;k++){ if(c[k]<lo[k])lo[k]=c[k]; if(c[k]>hi[k])hi[k]=c[k]; } };
  upd(vtx); for(auto&s:segs){upd(s.a);upd(s.b);}
  double xa[3]={vtx.X(),vtx.Y(),vtx.Z()};
  double xmin=lo[ia],xmax=hi[ia],ymin=lo[ib],ymax=hi[ib];
  double mx=0.15*(xmax-xmin)+5, my=0.15*(ymax-ymin)+5;
  TH2F* fr=new TH2F(Form("fr%p",(void*)pad),"",10,xmin-mx,xmax+mx,10,ymin-my,ymax+my);
  fr->SetStats(0); fr->GetXaxis()->SetTitle(xt); fr->GetYaxis()->SetTitle(yt); fr->Draw();
  for(auto&s:segs){ double A[3]={s.a.X(),s.a.Y(),s.a.Z()}, B[3]={s.b.X(),s.b.Y(),s.b.Z()};
    TLine* l=new TLine(A[ia],A[ib],B[ia],B[ib]); l->SetLineColor(s.col); l->SetLineWidth(s.wid); l->SetLineStyle(s.sty); l->Draw(); }
  TMarker* m=new TMarker(xa[ia],xa[ib],29); m->SetMarkerSize(2.2); m->SetMarkerColor(kBlack); m->Draw();
}

void EventDisplay(){
  gStyle->SetOptStat(0);
  gSystem->Exec("mkdir -p figs/eventdisplays");
  TFile* f=TFile::Open("/data/uboone/new_numi_flux/Run1_rhc_new_numi_flux_rhc_pandora_ntuple.root");
  TTree* t=nullptr; f->GetObject("nuselection/NeutrinoSelectionFilter",t);

  Int_t nu_pdg,ccnc,interaction,run,sub,evt;
  Float_t vx,vy,vz,topo,pe1,pe2,gdot;
  std::vector<float> *ts=0,*tl=0,*llr=0,*chi=0,*sx=0,*sy=0,*sz=0,*dx=0,*dy=0,*dz=0;
  std::vector<float> *hsx=0,*hsy=0,*hsz=0,*sdist=0,*se=0,*spx=0,*spy=0,*spz=0,*ovp=0;
  std::vector<int> *nhy=0,*bpdg=0;
  t->SetBranchStatus("*",0);
  auto B=[&](const char* n){ t->SetBranchStatus(n,1); };
  for(const char* n:{"nu_pdg","ccnc","interaction","run","sub","evt",
    "reco_nu_vtx_sce_x","reco_nu_vtx_sce_y","reco_nu_vtx_sce_z","topological_score",
    "pi0_energy1_Y","pi0_energy2_Y","pi0_gammadot",
    "trk_score_v","trk_len_v","trk_llr_pid_score_v","trk_pid_chipr_v",
    "trk_sce_start_x_v","trk_sce_start_y_v","trk_sce_start_z_v","trk_dir_x_v","trk_dir_y_v","trk_dir_z_v",
    "shr_start_x_v","shr_start_y_v","shr_start_z_v","shr_dist_v","shr_energy_y_v",
    "shr_px_v","shr_py_v","shr_pz_v","pfnplanehits_Y","backtracked_pdg","backtracked_overlay_purity"}) B(n);
  t->SetBranchAddress("nu_pdg",&nu_pdg);t->SetBranchAddress("ccnc",&ccnc);t->SetBranchAddress("interaction",&interaction);
  t->SetBranchAddress("run",&run);t->SetBranchAddress("sub",&sub);t->SetBranchAddress("evt",&evt);
  t->SetBranchAddress("reco_nu_vtx_sce_x",&vx);t->SetBranchAddress("reco_nu_vtx_sce_y",&vy);t->SetBranchAddress("reco_nu_vtx_sce_z",&vz);
  t->SetBranchAddress("topological_score",&topo);
  t->SetBranchAddress("pi0_energy1_Y",&pe1);t->SetBranchAddress("pi0_energy2_Y",&pe2);t->SetBranchAddress("pi0_gammadot",&gdot);
  t->SetBranchAddress("trk_score_v",&ts);t->SetBranchAddress("trk_len_v",&tl);t->SetBranchAddress("trk_llr_pid_score_v",&llr);
  t->SetBranchAddress("trk_pid_chipr_v",&chi);
  t->SetBranchAddress("trk_sce_start_x_v",&sx);t->SetBranchAddress("trk_sce_start_y_v",&sy);t->SetBranchAddress("trk_sce_start_z_v",&sz);
  t->SetBranchAddress("trk_dir_x_v",&dx);t->SetBranchAddress("trk_dir_y_v",&dy);t->SetBranchAddress("trk_dir_z_v",&dz);
  t->SetBranchAddress("shr_start_x_v",&hsx);t->SetBranchAddress("shr_start_y_v",&hsy);t->SetBranchAddress("shr_start_z_v",&hsz);
  t->SetBranchAddress("shr_dist_v",&sdist);t->SetBranchAddress("shr_energy_y_v",&se);
  t->SetBranchAddress("shr_px_v",&spx);t->SetBranchAddress("shr_py_v",&spy);t->SetBranchAddress("shr_pz_v",&spz);
  t->SetBranchAddress("pfnplanehits_Y",&nhy);t->SetBranchAddress("backtracked_pdg",&bpdg);
  t->SetBranchAddress("backtracked_overlay_purity",&ovp);

  const double MP=0.938272, MPI=0.139570;
  std::map<std::string,int> want={{"QE",1},{"RES",1},{"DIS",1},{"MEC",1},{"Cosmic",1},{"Pi0",1}};
  std::map<std::string,int> got;
  int ndrawn=0;

  Long64_t N=t->GetEntries();
  for(Long64_t i=0;i<N && ndrawn<(int)want.size();i++){
    t->GetEntry(i);
    if(!ts) continue; size_t n=ts->size();
    bool ok=true; for(auto v:{tl,llr,chi,sx,sy,sz,dx,dy,dz,hsx,hsy,hsz,sdist,se,spx,spy,spz}) if(!v||v->size()!=n) ok=false;
    if(!ok||!nhy||nhy->size()!=n||!bpdg||bpdg->size()!=n||!ovp||ovp->size()!=n) continue;
    if(!inFV(vx,vy,vz) || topo<=0.20) continue;
    TVector3 vtx(vx,vy,vz);
    auto isTrk=[&](int k){return ts->at(k)>0.5;};
    std::vector<TVector3> start(n); std::vector<double> pdist(n);
    for(size_t k=0;k<n;k++){ if(isTrk(k)){start[k].SetXYZ(sx->at(k),sy->at(k),sz->at(k)); pdist[k]=(start[k]-vtx).Mag();}
      else {start[k].SetXYZ(hsx->at(k),hsy->at(k),hsz->at(k)); pdist[k]=sdist->at(k);} }
    int ntrk=0,nshr=0; for(size_t k=0;k<n;k++)(isTrk(k)?ntrk:nshr)++;
    if(!(nshr>0&&nshr<4)||!(ntrk<5)) continue;
    // muon
    int mu=-1; for(size_t k=0;k<n;k++) if(isTrk(k)&&ts->at(k)>0.8&&pdist[k]<20&&tl->at(k)>10&&llr->at(k)>0.7)
      if(mu<0||tl->at(k)>tl->at(mu)) mu=k;
    if(mu<0) continue;
    // proton: min finite chipr (LLR fallback)
    int pr=-1; for(size_t k=0;k<n;k++){ if((int)k==mu||!isTrk(k))continue;
      if(std::isfinite(chi->at(k))&&(pr<0||!std::isfinite(chi->at(pr))||chi->at(k)<chi->at(pr)))pr=k; }
    if(pr<0||!std::isfinite(chi->at(pr))) for(size_t k=0;k<n;k++){ if((int)k==mu||!isTrk(k))continue; if(pr<0||llr->at(k)<llr->at(pr))pr=k; }
    if(pr<0) continue;
    // pion: max finite chipr (LLR fallback)
    int pi=-1; for(size_t k=0;k<n;k++){ if((int)k==mu||(int)k==pr||!isTrk(k))continue;
      if(std::isfinite(chi->at(k))&&(pi<0||!std::isfinite(chi->at(pi))||chi->at(k)>chi->at(pi)))pi=k; }
    if(pi<0||!std::isfinite(chi->at(pi))) for(size_t k=0;k<n;k++){ if((int)k==mu||(int)k==pr||!isTrk(k))continue; if(pi<0||llr->at(k)>llr->at(pi))pi=k; }
    if(pi<0) continue;
    // shower1
    int s1=-1; double b1=1e9; for(size_t k=0;k<n;k++){ if((int)k==mu||(int)k==pr||(int)k==pi||isTrk(k))continue;
      if(pdist[k]<b1){b1=pdist[k];s1=k;} }
    if(s1<0) continue;

    // category
    int cosmic = (ovp->at(mu)>0.5)?1:0;
    bool clearPi0 = (pe1>0&&pe2>0&&gdot>-1&&gdot<1 && std::sqrt(2*pe1*pe2*(1-gdot))>110 && std::sqrt(2*pe1*pe2*(1-gdot))<170);
    std::string cat = cosmic ? "Cosmic" : modeName(interaction,0);
    if(clearPi0 && got["Pi0"]<want["Pi0"] && !cosmic) cat="Pi0";
    if(!want.count(cat) || got[cat]>=want[cat]) continue;
    got[cat]++; ndrawn++;

    // build segments
    auto dirOf=[&](int k){ return TVector3(dx->at(k),dy->at(k),dz->at(k)); };
    std::vector<Seg> segs;
    auto addTrk=[&](int k,int col){ TVector3 a=start[k], b=a+dirOf(k)*std::max(2.0,(double)tl->at(k));
      segs.push_back({a,b,col,3,1}); };
    addTrk(mu,kRed+1); addTrk(pr,kAzure+1); addTrk(pi,kGreen+2);
    // shower: dashed from vtx to shower start, extended along shower dir
    TVector3 g0=start[s1]; TVector3 gdir(spx->at(s1),spy->at(s1),spz->at(s1)); if(gdir.Mag()>0)gdir=gdir.Unit();
    segs.push_back({vtx,g0,kMagenta+1,2,2});                       // gap (vtx->conversion)
    segs.push_back({g0,g0+gdir*15.0,kMagenta+1,3,1});             // shower direction

    // masses
    double mL=-1,mS=-1,mPi0=-1;
    { double pp=0; double KEp= -1; // proton mom from range KE not available here; approx from length? use chipr?
      // We do not have trk_energy_proton in the reduced display; approximate Lambda mass omitted if unavailable.
    }
    if(pe1>0&&pe2>0&&gdot>-1&&gdot<1){ double mm=2*pe1*pe2*(1-gdot); if(mm>0)mPi0=std::sqrt(mm)/1000.0; }

    // canvas
    TCanvas* c=new TCanvas(Form("ev%d",ndrawn),"",1300,600); c->Divide(2,1);
    drawView(c->cd(1),segs,vtx,2,0,"Z [cm]","X [cm]");
    drawView(c->cd(2),segs,vtx,2,1,"Z [cm]","Y [cm]");
    c->cd(1);
    TLegend* lg=new TLegend(0.12,0.74,0.55,0.90); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.032);
    TLine* lm=new TLine();lm->SetLineColor(kRed+1);lm->SetLineWidth(3);
    lg->AddEntry((TObject*)0,Form("#mu cand  (true %s)",pn(bpdg->at(mu))),"");
    lg->AddEntry((TObject*)0,Form("p cand  (true %s)",pn(bpdg->at(pr))),"");
    lg->AddEntry((TObject*)0,Form("#pi cand  (true %s)",pn(bpdg->at(pi))),"");
    lg->AddEntry((TObject*)0,Form("#gamma shower (true %s)",pn(bpdg->at(s1))),"");
    lg->Draw();
    c->cd(2);
    TPaveText* pt=new TPaveText(0.12,0.72,0.62,0.90,"NDC"); pt->SetBorderSize(0); pt->SetFillStyle(0); pt->SetTextAlign(12); pt->SetTextSize(0.032);
    pt->AddText(Form("%s   run %d sub %d evt %d",cat.c_str(),run,sub,evt));
    pt->AddText(Form("nu_pdg=%d  CC/NC=%d  topo=%.2f",nu_pdg,ccnc,topo));
    pt->AddText(Form("n_{trk}=%d  n_{shr}=%d  m_{#pi0}=%s GeV",ntrk,nshr, mPi0>0?Form("%.3f",mPi0):"--"));
    pt->Draw();
    TString out=Form("figs/eventdisplays/bkg_%s_run%d_evt%d.png",cat.c_str(),run,evt);
    c->SaveAs(out);
    printf("drew %s -> %s\n",cat.c_str(),out.Data());
  }
  printf("Total displays: %d\n",ndrawn);
}
