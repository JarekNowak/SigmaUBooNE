// ============================================================
//  SelOptScan.C
//  Diagnostic: where is signal lost in the pre-selection, and which
//  loosenings recover signal while the BDTg>0.96 cut keeps background
//  at the same level?
//
//  For each selection "variant" we re-run the FULL candidate selection
//  with one parameter changed, re-evaluate BDTg, and tabulate the
//  weighted signal (S) and background (B = other hyperons + bkg + dirt)
//  at two levels: pre-selection (Cuts[10]) and final (BDTg>0.96).
//
//  Also dumps the TRUE Lambda flight distance for signal events to test
//  whether the muon<->Lambda-daughter distance cuts are the bottleneck.
//
//  Usage:  cd /home/t2k/nowak/MicroBooNE/Sigma
//          root -b -l -q SelOptScan.C
// ============================================================

#define Sigma_cxx
#include "Sigma.h"
#include "TMVA/Reader.h"
#include "TChain.h"
#include "TH1F.h"
#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

template<typename T>
static inline float SafeGet(const vector<T>* v, int i, float s=-999.f)
{ return (!v||i<0||i>=(int)v->size()) ? s : (float)v->at(i); }
template<typename T>
static inline int SafeGetI(const vector<T>* v, int i, int s=-999)
{ return (!v||i<0||i>=(int)v->size()) ? s : (int)v->at(i); }

// ---- Selection parameters (nominal values from Sigma.C) ----
struct Params {
  string name;
  float muonScore   = 0.8f;   // muon trk_shr_score >
  float muonLLR     = 0.7f;   // muon trk_llrpid    >
  float muonDist    = 20.f;   // muon start within  (cm)
  float muonLen     = 10.f;   // muon trk_length    > (cm)
  float dMuPiMax    = 30.f;   // distMuonPion   <
  float dPrMuMax    = 20.f;   // distProtonMuon <
  float dPrPiMax    = 20.f;   // distProtonPion <
  float dShr1Max    = 100.f;  // shower1 to mu/pi/pr <
  int   shr1HitsMin = 10;     // 1st shower hits range
  int   shr1HitsMax = 300;
  int   nShrMax     = 4;      // nshowers in (0, nShrMax)
  int   nTrkMax     = 5;      // ntracks < nTrkMax
};

// globals for the TMVA reader
Float_t g_MuonTrkLength,g_MuonTrkDistVtx,g_MuonTrkScore,g_MuonNHitsY,g_MuonPhotonPID;
Float_t g_PionPID,g_PionTrkLength,g_PionTrkDistVtx,g_PionTrkScore,g_PionNHitsY,g_PionDeDxP2,g_PionPhotonPID;
Float_t g_ProtonPID,g_ProtonTrkLength,g_ProtonTrkDistVtx,g_ProtonTrkScore,g_ProtonNHitsY,g_ProtonDeDxP2,g_ProtonPhotonPID;
Float_t g_Shr1NHitsY,g_Shr1VtxDist,g_Shr1OpenAngle,g_Shr2NHitsY,g_Shr2VtxDist;
Float_t g_AngleMuonPion,g_AngleProtonPion,g_AngleProtonMuon,g_AngleShr1Muon,g_AngleShr1Pion,g_AngleShr1Proton;
Float_t g_DistMuonPion,g_DistProtonPion,g_DistProtonMuon,g_Weight;

// ============================================================
void SelOptScan()
// ============================================================
{
  // ---- chain: 0=signal/hyperon, 1=bkg, 2=dirt ----
  TChain *chain = new TChain("OutputTree","");
  chain->Add("/data/nowak/Sigma/analysisOutputRHC_cttest_Hyperon_ana_hyperon_sigmazero_mSlice_prodgenie_run3b.root/OutputTree");
  chain->Add("/data/nowak/Sigma/analysisOutputRHC_cttest_Background_prodgenie_numi_overlay_rhc_run3b_reco2.root/OutputTree");
  chain->Add("/data/nowak/Sigma/analysisOutputRHC_cttest_Dirt_prodgenie_numi_uboone_overlay_dirt_rhc_mcc9_run3b.root/OutputTree");
  Sigma S(chain);

  // ---- PoT scaling (from Sigma.C) ----
  double PoT[5] = {552.0, 552.0, 9.280, 17.30, 4.90};
  double scale[5]; for (int i=0;i<4;i++) scale[i]=PoT[4]/PoT[i]; scale[4]=1.0;

  // ---- TMVA reader (same 33 vars + spectator) ----
  TMVA::Reader *R = new TMVA::Reader("!Color:!Silent");
  R->AddVariable("MuonTrkLength",&g_MuonTrkLength);   R->AddVariable("MuonTrkDistVtx",&g_MuonTrkDistVtx);
  R->AddVariable("MuonTrkScore",&g_MuonTrkScore);     R->AddVariable("MuonNHitsY",&g_MuonNHitsY);
  R->AddVariable("MuonPhotonPID",&g_MuonPhotonPID);   R->AddVariable("PionPID",&g_PionPID);
  R->AddVariable("PionTrkLength",&g_PionTrkLength);   R->AddVariable("PionTrkDistVtx",&g_PionTrkDistVtx);
  R->AddVariable("PionTrkScore",&g_PionTrkScore);     R->AddVariable("PionNHitsY",&g_PionNHitsY);
  R->AddVariable("PionDeDxPlane2",&g_PionDeDxP2);     R->AddVariable("PionPhotonPID",&g_PionPhotonPID);
  R->AddVariable("ProtonPID",&g_ProtonPID);           R->AddVariable("ProtonTrkLength",&g_ProtonTrkLength);
  R->AddVariable("ProtonTrkDistVtx",&g_ProtonTrkDistVtx); R->AddVariable("ProtonTrkScore",&g_ProtonTrkScore);
  R->AddVariable("ProtonNHitsY",&g_ProtonNHitsY);     R->AddVariable("ProtonDeDxPlane2",&g_ProtonDeDxP2);
  R->AddVariable("ProtonPhotonPID",&g_ProtonPhotonPID); R->AddVariable("Shr1NHitsY",&g_Shr1NHitsY);
  R->AddVariable("Shr1VtxDist",&g_Shr1VtxDist);       R->AddVariable("Shr1OpenAngle",&g_Shr1OpenAngle);
  R->AddVariable("Shr2NHitsY",&g_Shr2NHitsY);         R->AddVariable("Shr2VtxDist",&g_Shr2VtxDist);
  R->AddVariable("AngleMuonPion",&g_AngleMuonPion);   R->AddVariable("AngleProtonPion",&g_AngleProtonPion);
  R->AddVariable("AngleProtonMuon",&g_AngleProtonMuon); R->AddVariable("AngleShr1Muon",&g_AngleShr1Muon);
  R->AddVariable("AngleShr1Pion",&g_AngleShr1Pion);   R->AddVariable("AngleShr1Proton",&g_AngleShr1Proton);
  R->AddVariable("DistMuonPion",&g_DistMuonPion);     R->AddVariable("DistProtonPion",&g_DistProtonPion);
  R->AddVariable("DistProtonMuon",&g_DistProtonMuon); R->AddSpectator("Weight",&g_Weight);
  R->BookMVA("BDTg","SigmaDataset/weights/TMVASigma_BDTg.weights.xml");

  // ---- variants ----
  vector<Params> V;
  { Params p; p.name="0 nominal";                         V.push_back(p); }
  { Params p; p.name="1 muonScore 0.8->0.5"; p.muonScore=0.5f; V.push_back(p); }
  { Params p; p.name="2 muonLLR 0.7->0.0";   p.muonLLR=0.0f;   V.push_back(p); }
  { Params p; p.name="3 muonDist 20->30";    p.muonDist=30.f;  V.push_back(p); }
  { Params p; p.name="4 dPrMu/dPrPi 20->30"; p.dPrMuMax=30.f; p.dPrPiMax=30.f; V.push_back(p); }
  { Params p; p.name="5 dMuPi30/dPrMu/dPrPi ->50/40/40"; p.dMuPiMax=50.f; p.dPrMuMax=40.f; p.dPrPiMax=40.f; V.push_back(p); }
  { Params p; p.name="6 shr1HitsMax 300->600"; p.shr1HitsMax=600; V.push_back(p); }
  { Params p; p.name="7 nTrkMax 5->6";       p.nTrkMax=6;      V.push_back(p); }
  const int NV = V.size();

  // accumulators: [variant] -> weighted S/B at pre & final
  vector<double> Spre(NV,0), Bpre(NV,0), Sfin(NV,0), Bfin(NV,0);
  // raw signal counts (unweighted) for stat sanity
  vector<long>   NSpre(NV,0), NSfin(NV,0);

  // true Lambda flight distance for signal: lost-by-distance vs kept
  TH1F *hFlightAll  = new TH1F("hFlightAll","true #Lambda flight;flight [cm];signal evts",60,0,60);
  TH1F *hFlightLost = new TH1F("hFlightLost","lost by dist cut",60,0,60);

  // ---- selection lambda: returns final pass, sets pre via ref ----
  auto runSel = [&](const Params& P, bool& passPre)->float {
    passPre=false;
    TVector3 vtx(S.reco_primary_vtx_x,S.reco_primary_vtx_y,S.reco_primary_vtx_z);
    if (!S.reco_primary_vtx_inFV) return -999.f;

    int ntracks=0,nshowers=0;
    for (int i=0;i<(int)S.pfp_trk_shr_score->size();i++)
      (S.pfp_trk_shr_score->at(i)<0.5)?nshowers++:ntracks++;
    if (!(nshowers>0 && nshowers<P.nShrMax)) return -999.f;

    int iMu=-1,iPi=-1,iPr=-1,iS1=-1,iS2=-1;
    bool selMu=false,selPi=false,selPr=false;
    for (int i=0;i<(int)S.pfp_trk_shr_score->size();i++){
      TVector3 st(S.trk_start_x->at(i),S.trk_start_y->at(i),S.trk_start_z->at(i));
      if (S.pfp_trk_shr_score->at(i)>P.muonScore && (st-vtx).Mag()<P.muonDist
          && S.trk_length->at(i)>P.muonLen && S.trk_llrpid->at(i)>P.muonLLR){
        selMu=true;
        if (iMu==-1||S.trk_length->at(i)>S.trk_length->at(iMu)) iMu=i;
      }
    }
    if (!selMu) return -999.f;
    for (int i=0;i<(int)S.pfp_trk_shr_score->size();i++){
      if (i==iMu) continue;
      if (S.pfp_trk_shr_score->at(i)>0.5){ selPi=true;
        if (iPi==-1||S.trk_llrpid->at(i)>S.trk_llrpid->at(iPi)) iPi=i; }
    }
    if (!selPi) return -999.f;
    for (int i=0;i<(int)S.pfp_trk_shr_score->size();i++){
      if (i==iMu||i==iPi) continue;
      if (S.pfp_trk_shr_score->at(i)>0.5){ selPr=true;
        if (iPr==-1||S.trk_llrpid->at(i)<S.trk_llrpid->at(iPr)) iPr=i; }
    }
    if (!selPr) return -999.f;
    if (!(ntracks<P.nTrkMax)) return -999.f;

    double b1=1e9,b2=1e9;
    for (int i=0;i<(int)S.pfp_trk_shr_score->size();i++){
      if (i==iMu||i==iPi||i==iPr) continue;
      if (S.pfp_trk_shr_score->at(i)>=0.5) continue;
      TVector3 st(S.trk_start_x->at(i),S.trk_start_y->at(i),S.trk_start_z->at(i));
      double d=(st-vtx).Mag();
      if (d<b1){b2=b1;iS2=iS1;b1=d;iS1=i;} else if(d<b2){b2=d;iS2=i;}
    }
    bool sel1stShr=false;
    if (iS1>-1){int h=S.pfpnplanehits_Y->at(iS1); sel1stShr=(h>P.shr1HitsMin&&h<P.shr1HitsMax);}
    if (!sel1stShr) return -999.f;

    TVector3 muSt,piSt,prSt,s1St,s2St,muD,piD,prD,s1D;
    if(iMu>-1){muSt.SetXYZ(S.trk_start_x->at(iMu),S.trk_start_y->at(iMu),S.trk_start_z->at(iMu));
               muD.SetXYZ(S.trk_dir_x->at(iMu),S.trk_dir_y->at(iMu),S.trk_dir_z->at(iMu));}
    if(iPi>-1){piSt.SetXYZ(S.trk_start_x->at(iPi),S.trk_start_y->at(iPi),S.trk_start_z->at(iPi));
               piD.SetXYZ(S.trk_dir_x->at(iPi),S.trk_dir_y->at(iPi),S.trk_dir_z->at(iPi));}
    if(iPr>-1){prSt.SetXYZ(S.trk_start_x->at(iPr),S.trk_start_y->at(iPr),S.trk_start_z->at(iPr));
               prD.SetXYZ(S.trk_dir_x->at(iPr),S.trk_dir_y->at(iPr),S.trk_dir_z->at(iPr));}
    if(iS1>-1){s1St.SetXYZ(S.trk_start_x->at(iS1),S.trk_start_y->at(iS1),S.trk_start_z->at(iS1));
               s1D.SetXYZ(S.trk_dir_x->at(iS1),S.trk_dir_y->at(iS1),S.trk_dir_z->at(iS1));}
    if(iS2>-1){s2St.SetXYZ(S.trk_start_x->at(iS2),S.trk_start_y->at(iS2),S.trk_start_z->at(iS2));}

    double dMuPi=(muSt-piSt).Mag(), dPrPi=(prSt-piSt).Mag(), dPrMu=(prSt-muSt).Mag();
    double dS1Mu=(iS1>-1)?(s1St-muSt).Mag():9e4, dS1Pi=(iS1>-1)?(s1St-piSt).Mag():9e4,
           dS1Pr=(iS1>-1)?(s1St-prSt).Mag():9e4;
    bool selDist=(dMuPi<P.dMuPiMax && dPrMu<P.dPrMuMax && dPrPi<P.dPrPiMax &&
                  dS1Mu<P.dShr1Max && dS1Pi<P.dShr1Max && dS1Pr<P.dShr1Max);
    if (!selDist) return -999.f;

    passPre=true;
    // fill reader
    g_MuonTrkLength=SafeGet(S.trk_length,iMu); g_MuonTrkDistVtx=(Float_t)(muSt-vtx).Mag();
    g_MuonTrkScore=SafeGet(S.pfp_trk_shr_score,iMu); g_MuonNHitsY=(Float_t)SafeGetI(S.pfpnplanehits_Y,iMu);
    g_MuonPhotonPID=SafeGet(S.pfp_photon_bdt_score,iMu);
    g_PionPID=SafeGet(S.trk_llrpid,iPi); g_PionTrkLength=SafeGet(S.trk_length,iPi);
    g_PionTrkDistVtx=(Float_t)(piSt-vtx).Mag(); g_PionTrkScore=SafeGet(S.pfp_trk_shr_score,iPi);
    g_PionNHitsY=(Float_t)SafeGetI(S.pfpnplanehits_Y,iPi); g_PionDeDxP2=SafeGet(S.trk_mean_dedx_plane2,iPi);
    g_PionPhotonPID=SafeGet(S.pfp_photon_bdt_score,iPi);
    g_ProtonPID=SafeGet(S.trk_llrpid,iPr); g_ProtonTrkLength=SafeGet(S.trk_length,iPr);
    g_ProtonTrkDistVtx=(Float_t)(prSt-vtx).Mag(); g_ProtonTrkScore=SafeGet(S.pfp_trk_shr_score,iPr);
    g_ProtonNHitsY=(Float_t)SafeGetI(S.pfpnplanehits_Y,iPr); g_ProtonDeDxP2=SafeGet(S.trk_mean_dedx_plane2,iPr);
    g_ProtonPhotonPID=SafeGet(S.pfp_photon_bdt_score,iPr);
    g_Shr1NHitsY=(Float_t)SafeGetI(S.pfpnplanehits_Y,iS1); g_Shr1VtxDist=(Float_t)(s1St-vtx).Mag();
    g_Shr1OpenAngle=SafeGet(S.shr_open_angle,iS1);
    g_Shr2NHitsY=(iS2>-1)?(Float_t)SafeGetI(S.pfpnplanehits_Y,iS2):-999.f;
    g_Shr2VtxDist=(iS2>-1)?(Float_t)(s2St-vtx).Mag():-999.f;
    g_AngleMuonPion=(Float_t)muD.Angle(piD); g_AngleProtonPion=(Float_t)prD.Angle(piD);
    g_AngleProtonMuon=(Float_t)prD.Angle(muD);
    g_AngleShr1Muon=(iS1>-1)?(Float_t)s1D.Angle(muD):-999.f;
    g_AngleShr1Pion=(iS1>-1)?(Float_t)s1D.Angle(piD):-999.f;
    g_AngleShr1Proton=(iS1>-1)?(Float_t)s1D.Angle(prD):-999.f;
    g_DistMuonPion=(Float_t)dMuPi; g_DistProtonPion=(Float_t)dPrPi; g_DistProtonMuon=(Float_t)dPrMu;
    g_Weight=1.0f;
    return R->EvaluateMVA("BDTg");
  };

  Long64_t n = chain->GetEntries();
  cout<<"Entries: "<<n<<endl;
  for (Long64_t je=0; je<n; je++){
    S.LoadTree(je); S.fChain->GetEntry(je);
    if (je%100000==0) cout<<"\r "<<100*je/n<<"%"<<flush;

    // truth + sample
    double pmom=0,pimom=0; int mcdec=0;
    for (int i=0;i<(int)S.mc_decay_pdg->size();i++){
      if (S.mc_decay_pdg->at(i)==2212){mcdec++;pmom=SafeGet(S.mc_decay_mom,i);}
      if (S.mc_decay_pdg->at(i)==-211){mcdec++;pimom=SafeGet(S.mc_decay_mom,i);}
    }
    bool inFV=(S.mc_nu_pos_x>3.00&&S.mc_nu_pos_x<253.35&&S.mc_nu_pos_y>-112.53&&
               S.mc_nu_pos_y<114.47&&S.mc_nu_pos_z>3.10&&S.mc_nu_pos_z<1033.0);
    bool dead=(S.mc_nu_pos_z>675.1&&S.mc_nu_pos_z<775.1);
    bool IsSignal=(S.mc_nu_pdg==-14)&&(S.mc_hyperon_pdg==3212)&&(mcdec==2)&&pmom>0.3&&pimom>0.1;
    int tn=S.fChain->GetTreeNumber();
    int s = (tn==0)?(IsSignal?0:1):(tn==1?2:3);
    if ((s==0||s==1)&&(!inFV||dead)) continue;

    double w=scale[s];
    bool isSig=(s==0);
    bool isBkg=(s==1||s==2||s==3);

    // true Lambda flight (signal only): decay start - primary vtx
    double flight=-1;
    if (isSig && S.mc_decay_start_x->size()>0){
      double dx=S.mc_decay_start_x->at(0)-S.mc_nu_pos_x;
      double dy=S.mc_decay_start_y->at(0)-S.mc_nu_pos_y;
      double dz=S.mc_decay_start_z->at(0)-S.mc_nu_pos_z;
      flight=sqrt(dx*dx+dy*dy+dz*dz);
    }

    for (int v=0; v<NV; v++){
      bool passPre=false;
      float bdtg=runSel(V[v],passPre);
      bool passFin = passPre && (bdtg>0.96f);
      if (passPre){ if(isSig){Spre[v]+=w;NSpre[v]++;} else if(isBkg)Bpre[v]+=w; }
      if (passFin){ if(isSig){Sfin[v]+=w;NSfin[v]++;} else if(isBkg)Bfin[v]+=w; }

      // flight diagnostic only under nominal
      if (v==0 && isSig && flight>=0){
        hFlightAll->Fill(flight);
        if (!passPre){
          // was it the distance cut specifically? re-run with relaxed dist
          Params relaxed=V[0]; relaxed.dMuPiMax=999; relaxed.dPrMuMax=999; relaxed.dPrPiMax=999; relaxed.dShr1Max=9999;
          bool pp=false; runSel(relaxed,pp);
          if (pp) hFlightLost->Fill(flight); // passes everything except distance
        }
      }
    }
  }
  cout<<"\rdone     "<<endl;

  // ---- report ----
  cout<<fixed<<setprecision(2);
  cout<<"\n================ FINAL (BDTg>0.96) ================\n";
  cout<<setw(34)<<left<<"variant"<<right
      <<setw(10)<<"S"<<setw(10)<<"B"<<setw(10)<<"dS%"<<setw(10)<<"dB%"
      <<setw(11)<<"S/sqrtB"<<setw(9)<<"NS_raw"<<"\n";
  double S0=Sfin[0],B0=Bfin[0];
  for (int v=0;v<NV;v++){
    double ssb=(Bfin[v]>0)?Sfin[v]/sqrt(Bfin[v]):0;
    cout<<setw(34)<<left<<V[v].name<<right
        <<setw(10)<<Sfin[v]<<setw(10)<<Bfin[v]
        <<setw(10)<<100*(Sfin[v]-S0)/S0<<setw(10)<<(B0>0?100*(Bfin[v]-B0)/B0:0)
        <<setw(11)<<ssb<<setw(9)<<NSfin[v]<<"\n";
  }
  cout<<"\n================ PRE-SELECTION (Cuts[10]) ================\n";
  cout<<setw(34)<<left<<"variant"<<right
      <<setw(10)<<"S"<<setw(10)<<"B"<<setw(10)<<"dS%"<<setw(10)<<"dB%"<<"\n";
  double Sp0=Spre[0],Bp0=Bpre[0];
  for (int v=0;v<NV;v++)
    cout<<setw(34)<<left<<V[v].name<<right
        <<setw(10)<<Spre[v]<<setw(10)<<Bpre[v]
        <<setw(10)<<100*(Spre[v]-Sp0)/Sp0<<setw(10)<<(Bp0>0?100*(Bpre[v]-Bp0)/Bp0:0)<<"\n";

  cout<<"\n========= TRUE Lambda flight distance (signal) =========\n";
  cout<<"  mean(all)="<<hFlightAll->GetMean()<<" cm   "
      <<"integral(all)="<<hFlightAll->Integral()<<"\n";
  cout<<"  lost ONLY by distance cut: integral="<<hFlightLost->Integral()
      <<"  mean="<<hFlightLost->GetMean()<<" cm\n";
  cout<<"  bin  flight[cm]   all   lost_by_dist\n";
  for (int b=1;b<=12;b++)
    cout<<"  "<<setw(3)<<b<<setw(10)<<hFlightAll->GetBinCenter(b)
        <<setw(8)<<hFlightAll->GetBinContent(b)
        <<setw(12)<<hFlightLost->GetBinContent(b)<<"\n";
}
