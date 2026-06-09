// ============================================================
//  BDTThreshScan.C
//  The pre-selection is NOT the efficiency bottleneck (see SelOptScan):
//  every loosened variant gives the identical final yield because the
//  BDTg>0.96 cut rejects the recovered events.  The real knob is the
//  BDTg threshold itself.
//
//  This macro runs the NOMINAL pre-selection once, stores the weighted
//  BDTg-score distributions for signal and background, then scans the
//  threshold to show S, B, efficiency, purity and S/sqrt(B).
//
//  Usage:  root -b -l -q BDTThreshScan.C
// ============================================================

#define Sigma_cxx
#include "Sigma.h"
#include "TMVA/Reader.h"
#include "TChain.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include <iostream>
#include <iomanip>

using namespace std;
template<typename T> static inline float SafeGet(const vector<T>* v,int i,float s=-999.f)
{ return (!v||i<0||i>=(int)v->size())?s:(float)v->at(i); }
template<typename T> static inline int SafeGetI(const vector<T>* v,int i,int s=-999)
{ return (!v||i<0||i>=(int)v->size())?s:(int)v->at(i); }

Float_t a_MuonTrkLength,a_MuonTrkDistVtx,a_MuonTrkScore,a_MuonNHitsY,a_MuonPhotonPID;
Float_t a_PionPID,a_PionTrkLength,a_PionTrkDistVtx,a_PionTrkScore,a_PionNHitsY,a_PionDeDxP2,a_PionPhotonPID;
Float_t a_ProtonPID,a_ProtonTrkLength,a_ProtonTrkDistVtx,a_ProtonTrkScore,a_ProtonNHitsY,a_ProtonDeDxP2,a_ProtonPhotonPID;
Float_t a_Shr1NHitsY,a_Shr1VtxDist,a_Shr1OpenAngle,a_Shr2NHitsY,a_Shr2VtxDist;
Float_t a_AngleMuonPion,a_AngleProtonPion,a_AngleProtonMuon,a_AngleShr1Muon,a_AngleShr1Pion,a_AngleShr1Proton;
Float_t a_DistMuonPion,a_DistProtonPion,a_DistProtonMuon,a_Weight;

void BDTThreshScan()
{
  gStyle->SetOptStat(0);
  TChain *chain=new TChain("OutputTree","");
  chain->Add("/data/nowak/Sigma/analysisOutputRHC_cttest_Hyperon_ana_hyperon_sigmazero_mSlice_prodgenie_run3b.root/OutputTree");
  chain->Add("/data/nowak/Sigma/analysisOutputRHC_cttest_Background_prodgenie_numi_overlay_rhc_run3b_reco2.root/OutputTree");
  chain->Add("/data/nowak/Sigma/analysisOutputRHC_cttest_Dirt_prodgenie_numi_uboone_overlay_dirt_rhc_mcc9_run3b.root/OutputTree");
  Sigma S(chain);

  double PoT[5]={552.0,552.0,9.280,17.30,4.90};
  double scale[5]; for(int i=0;i<4;i++)scale[i]=PoT[4]/PoT[i]; scale[4]=1.0;

  TMVA::Reader *R=new TMVA::Reader("!Color:!Silent");
  R->AddVariable("MuonTrkLength",&a_MuonTrkLength);   R->AddVariable("MuonTrkDistVtx",&a_MuonTrkDistVtx);
  R->AddVariable("MuonTrkScore",&a_MuonTrkScore);     R->AddVariable("MuonNHitsY",&a_MuonNHitsY);
  R->AddVariable("MuonPhotonPID",&a_MuonPhotonPID);   R->AddVariable("PionPID",&a_PionPID);
  R->AddVariable("PionTrkLength",&a_PionTrkLength);   R->AddVariable("PionTrkDistVtx",&a_PionTrkDistVtx);
  R->AddVariable("PionTrkScore",&a_PionTrkScore);     R->AddVariable("PionNHitsY",&a_PionNHitsY);
  R->AddVariable("PionDeDxPlane2",&a_PionDeDxP2);     R->AddVariable("PionPhotonPID",&a_PionPhotonPID);
  R->AddVariable("ProtonPID",&a_ProtonPID);           R->AddVariable("ProtonTrkLength",&a_ProtonTrkLength);
  R->AddVariable("ProtonTrkDistVtx",&a_ProtonTrkDistVtx); R->AddVariable("ProtonTrkScore",&a_ProtonTrkScore);
  R->AddVariable("ProtonNHitsY",&a_ProtonNHitsY);     R->AddVariable("ProtonDeDxPlane2",&a_ProtonDeDxP2);
  R->AddVariable("ProtonPhotonPID",&a_ProtonPhotonPID); R->AddVariable("Shr1NHitsY",&a_Shr1NHitsY);
  R->AddVariable("Shr1VtxDist",&a_Shr1VtxDist);       R->AddVariable("Shr1OpenAngle",&a_Shr1OpenAngle);
  R->AddVariable("Shr2NHitsY",&a_Shr2NHitsY);         R->AddVariable("Shr2VtxDist",&a_Shr2VtxDist);
  R->AddVariable("AngleMuonPion",&a_AngleMuonPion);   R->AddVariable("AngleProtonPion",&a_AngleProtonPion);
  R->AddVariable("AngleProtonMuon",&a_AngleProtonMuon); R->AddVariable("AngleShr1Muon",&a_AngleShr1Muon);
  R->AddVariable("AngleShr1Pion",&a_AngleShr1Pion);   R->AddVariable("AngleShr1Proton",&a_AngleShr1Proton);
  R->AddVariable("DistMuonPion",&a_DistMuonPion);     R->AddVariable("DistProtonPion",&a_DistProtonPion);
  R->AddVariable("DistProtonMuon",&a_DistProtonMuon); R->AddSpectator("Weight",&a_Weight);
  R->BookMVA("BDTg","SigmaDataset/weights/TMVASigma_BDTg.weights.xml");

  // weighted BDTg-score distributions, fine bins on [-1,1]
  const int NB=400;
  TH1D *hS=new TH1D("hS","BDTg score;BDTg;weighted events",NB,-1,1);
  TH1D *hB=new TH1D("hB","",NB,-1,1);
  hS->Sumw2(); hB->Sumw2();
  double denomSig=0; // weighted signal entering pre-selection denominator (true signal in FV)

  Long64_t n=chain->GetEntries();
  cout<<"Entries: "<<n<<endl;
  for(Long64_t je=0; je<n; je++){
    S.LoadTree(je); S.fChain->GetEntry(je);
    if(je%100000==0) cout<<"\r "<<100*je/n<<"%"<<flush;

    double pmom=0,pimom=0; int mcdec=0;
    for(int i=0;i<(int)S.mc_decay_pdg->size();i++){
      if(S.mc_decay_pdg->at(i)==2212){mcdec++;pmom=SafeGet(S.mc_decay_mom,i);}
      if(S.mc_decay_pdg->at(i)==-211){mcdec++;pimom=SafeGet(S.mc_decay_mom,i);}
    }
    bool inFV=(S.mc_nu_pos_x>3.00&&S.mc_nu_pos_x<253.35&&S.mc_nu_pos_y>-112.53&&
               S.mc_nu_pos_y<114.47&&S.mc_nu_pos_z>3.10&&S.mc_nu_pos_z<1033.0);
    bool dead=(S.mc_nu_pos_z>675.1&&S.mc_nu_pos_z<775.1);
    bool IsSignal=(S.mc_nu_pdg==-14)&&(S.mc_hyperon_pdg==3212)&&(mcdec==2)&&pmom>0.3&&pimom>0.1;
    int tn=S.fChain->GetTreeNumber();
    int s=(tn==0)?(IsSignal?0:1):(tn==1?2:3);
    if((s==0||s==1)&&(!inFV||dead)) continue;
    double w=scale[s];
    if(s==0) denomSig+=w;

    // ---- nominal pre-selection ----
    TVector3 vtx(S.reco_primary_vtx_x,S.reco_primary_vtx_y,S.reco_primary_vtx_z);
    if(!S.reco_primary_vtx_inFV) continue;
    int ntr=0,nsh=0;
    for(int i=0;i<(int)S.pfp_trk_shr_score->size();i++)(S.pfp_trk_shr_score->at(i)<0.5)?nsh++:ntr++;
    if(!(nsh>0&&nsh<4)) continue;
    int iMu=-1,iPi=-1,iPr=-1,iS1=-1,iS2=-1; bool selMu=false,selPi=false,selPr=false;
    for(int i=0;i<(int)S.pfp_trk_shr_score->size();i++){
      TVector3 st(S.trk_start_x->at(i),S.trk_start_y->at(i),S.trk_start_z->at(i));
      if(S.pfp_trk_shr_score->at(i)>0.8&&(st-vtx).Mag()<20&&S.trk_length->at(i)>10&&S.trk_llrpid->at(i)>0.7){
        selMu=true; if(iMu==-1||S.trk_length->at(i)>S.trk_length->at(iMu))iMu=i;}
    }
    if(!selMu)continue;
    for(int i=0;i<(int)S.pfp_trk_shr_score->size();i++){if(i==iMu)continue;
      if(S.pfp_trk_shr_score->at(i)>0.5){selPi=true;if(iPi==-1||S.trk_llrpid->at(i)>S.trk_llrpid->at(iPi))iPi=i;}}
    if(!selPi)continue;
    for(int i=0;i<(int)S.pfp_trk_shr_score->size();i++){if(i==iMu||i==iPi)continue;
      if(S.pfp_trk_shr_score->at(i)>0.5){selPr=true;if(iPr==-1||S.trk_llrpid->at(i)<S.trk_llrpid->at(iPr))iPr=i;}}
    if(!selPr)continue;
    if(!(ntr<5))continue;
    double b1=1e9,b2=1e9;
    for(int i=0;i<(int)S.pfp_trk_shr_score->size();i++){if(i==iMu||i==iPi||i==iPr)continue;
      if(S.pfp_trk_shr_score->at(i)>=0.5)continue;
      TVector3 st(S.trk_start_x->at(i),S.trk_start_y->at(i),S.trk_start_z->at(i));
      double d=(st-vtx).Mag(); if(d<b1){b2=b1;iS2=iS1;b1=d;iS1=i;}else if(d<b2){b2=d;iS2=i;}}
    bool s1=false; if(iS1>-1){int h=S.pfpnplanehits_Y->at(iS1);s1=(h>10&&h<300);}
    if(!s1)continue;
    TVector3 muSt,piSt,prSt,s1St,s2St,muD,piD,prD,s1D;
    muSt.SetXYZ(S.trk_start_x->at(iMu),S.trk_start_y->at(iMu),S.trk_start_z->at(iMu));
    muD.SetXYZ(S.trk_dir_x->at(iMu),S.trk_dir_y->at(iMu),S.trk_dir_z->at(iMu));
    piSt.SetXYZ(S.trk_start_x->at(iPi),S.trk_start_y->at(iPi),S.trk_start_z->at(iPi));
    piD.SetXYZ(S.trk_dir_x->at(iPi),S.trk_dir_y->at(iPi),S.trk_dir_z->at(iPi));
    prSt.SetXYZ(S.trk_start_x->at(iPr),S.trk_start_y->at(iPr),S.trk_start_z->at(iPr));
    prD.SetXYZ(S.trk_dir_x->at(iPr),S.trk_dir_y->at(iPr),S.trk_dir_z->at(iPr));
    s1St.SetXYZ(S.trk_start_x->at(iS1),S.trk_start_y->at(iS1),S.trk_start_z->at(iS1));
    s1D.SetXYZ(S.trk_dir_x->at(iS1),S.trk_dir_y->at(iS1),S.trk_dir_z->at(iS1));
    if(iS2>-1)s2St.SetXYZ(S.trk_start_x->at(iS2),S.trk_start_y->at(iS2),S.trk_start_z->at(iS2));
    double dMuPi=(muSt-piSt).Mag(),dPrPi=(prSt-piSt).Mag(),dPrMu=(prSt-muSt).Mag();
    double dS1Mu=(s1St-muSt).Mag(),dS1Pi=(s1St-piSt).Mag(),dS1Pr=(s1St-prSt).Mag();
    if(!(dMuPi<30&&dPrMu<20&&dPrPi<20&&dS1Mu<100&&dS1Pi<100&&dS1Pr<100))continue;

    a_MuonTrkLength=SafeGet(S.trk_length,iMu);a_MuonTrkDistVtx=(Float_t)(muSt-vtx).Mag();
    a_MuonTrkScore=SafeGet(S.pfp_trk_shr_score,iMu);a_MuonNHitsY=(Float_t)SafeGetI(S.pfpnplanehits_Y,iMu);
    a_MuonPhotonPID=SafeGet(S.pfp_photon_bdt_score,iMu);
    a_PionPID=SafeGet(S.trk_llrpid,iPi);a_PionTrkLength=SafeGet(S.trk_length,iPi);
    a_PionTrkDistVtx=(Float_t)(piSt-vtx).Mag();a_PionTrkScore=SafeGet(S.pfp_trk_shr_score,iPi);
    a_PionNHitsY=(Float_t)SafeGetI(S.pfpnplanehits_Y,iPi);a_PionDeDxP2=SafeGet(S.trk_mean_dedx_plane2,iPi);
    a_PionPhotonPID=SafeGet(S.pfp_photon_bdt_score,iPi);
    a_ProtonPID=SafeGet(S.trk_llrpid,iPr);a_ProtonTrkLength=SafeGet(S.trk_length,iPr);
    a_ProtonTrkDistVtx=(Float_t)(prSt-vtx).Mag();a_ProtonTrkScore=SafeGet(S.pfp_trk_shr_score,iPr);
    a_ProtonNHitsY=(Float_t)SafeGetI(S.pfpnplanehits_Y,iPr);a_ProtonDeDxP2=SafeGet(S.trk_mean_dedx_plane2,iPr);
    a_ProtonPhotonPID=SafeGet(S.pfp_photon_bdt_score,iPr);
    a_Shr1NHitsY=(Float_t)SafeGetI(S.pfpnplanehits_Y,iS1);a_Shr1VtxDist=(Float_t)(s1St-vtx).Mag();
    a_Shr1OpenAngle=SafeGet(S.shr_open_angle,iS1);
    a_Shr2NHitsY=(iS2>-1)?(Float_t)SafeGetI(S.pfpnplanehits_Y,iS2):-999.f;
    a_Shr2VtxDist=(iS2>-1)?(Float_t)(s2St-vtx).Mag():-999.f;
    a_AngleMuonPion=(Float_t)muD.Angle(piD);a_AngleProtonPion=(Float_t)prD.Angle(piD);
    a_AngleProtonMuon=(Float_t)prD.Angle(muD);a_AngleShr1Muon=(Float_t)s1D.Angle(muD);
    a_AngleShr1Pion=(Float_t)s1D.Angle(piD);a_AngleShr1Proton=(Float_t)s1D.Angle(prD);
    a_DistMuonPion=(Float_t)dMuPi;a_DistProtonPion=(Float_t)dPrPi;a_DistProtonMuon=(Float_t)dPrMu;
    a_Weight=1.0f;
    float bdtg=R->EvaluateMVA("BDTg");
    if(s==0) hS->Fill(bdtg,w); else hB->Fill(bdtg,w);
  }
  cout<<"\rdone     "<<endl;

  // ---- scan thresholds: S, B above threshold ----
  cout<<fixed<<setprecision(3);
  cout<<"\n  denominator (true signal, in FV): "<<denomSig<<" weighted events\n";
  cout<<"  total pre-selected: S="<<hS->Integral()<<"  B="<<hB->Integral()<<"\n\n";
  cout<<setw(8)<<"thr"<<setw(10)<<"S"<<setw(10)<<"B"
      <<setw(10)<<"eff%"<<setw(11)<<"purity%"<<setw(11)<<"S/sqrtB"<<setw(12)<<"S/sqrt(S+B)"<<"\n";
  double thrs[]={0.80,0.85,0.88,0.90,0.92,0.94,0.95,0.96,0.97,0.98};
  for(double t: thrs){
    int b0=hS->FindBin(t+1e-6);
    double Ss=hS->Integral(b0,NB), Bb=hB->Integral(b0,NB);
    double eff=100*Ss/denomSig;
    double pur=(Ss+Bb>0)?100*Ss/(Ss+Bb):0;
    double ssb=(Bb>0)?Ss/sqrt(Bb):0;
    double sssb=(Ss+Bb>0)?Ss/sqrt(Ss+Bb):0;
    cout<<setw(8)<<t<<setw(10)<<Ss<<setw(10)<<Bb
        <<setw(10)<<eff<<setw(11)<<pur<<setw(11)<<ssb<<setw(12)<<sssb
        <<(fabs(t-0.96)<1e-6?"   <-- current":"")<<"\n";
  }

  // overlay plot (area-normalised shapes)
  TCanvas *c=new TCanvas("c","BDTg",800,550);
  c->SetLogy();
  hB->SetLineColor(kRed+1); hB->SetFillColorAlpha(kRed+1,0.25); hB->SetLineWidth(2);
  hS->SetLineColor(kBlue+1);hS->SetFillColorAlpha(kBlue+1,0.25);hS->SetLineWidth(2);
  hB->GetYaxis()->SetTitle("weighted events / bin");
  hB->Draw("HIST"); hS->Draw("HIST SAME");
  TLegend *l=new TLegend(0.15,0.75,0.45,0.88); l->SetBorderSize(0);
  l->AddEntry(hS,"Signal","f"); l->AddEntry(hB,"Background","f"); l->Draw();
  c->SaveAs("figs/EfficiencyPlots/BDTg_ScoreSB.png");
  cout<<"\nSaved figs/EfficiencyPlots/BDTg_ScoreSB.png\n";
}
