// ============================================================
//  EvalCompare.C
//  Apples-to-apples comparison of the PRODUCTION BDTg (33 vars) vs the
//  retrained v2 BDTg (33 + LambdaMass + DistShr1{Muon,Pion,Proton}),
//  evaluated on the SAME pre-selected events (TreeS/TreeB in
//  TMVA_input_v3.root, whose Weight branch is the PoT scale factor).
//
//  Reports weighted signal S, background B, efficiency (S / true-signal
//  denominator) and S/sqrt(B) as a function of the BDTg threshold, for
//  both classifiers.
//
//  Usage:  root -b -l -q EvalCompare.C
// ============================================================
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMVA/Reader.h"
#include <iostream>
#include <iomanip>
using namespace std;

// true-signal denominator (weighted, in FV) measured by BDTThreshScan
static const double DENOM = 17.869;

// all variables needed by either reader
float MuonTrkLength,MuonTrkDistVtx,MuonTrkScore,MuonNHitsY,MuonPhotonPID;
float PionPID,PionTrkLength,PionTrkDistVtx,PionTrkScore,PionNHitsY,PionDeDxPlane2,PionPhotonPID;
float ProtonPID,ProtonTrkLength,ProtonTrkDistVtx,ProtonTrkScore,ProtonNHitsY,ProtonDeDxPlane2,ProtonPhotonPID;
float Shr1NHitsY,Shr1VtxDist,Shr1OpenAngle,Shr2NHitsY,Shr2VtxDist;
float AngleMuonPion,AngleProtonPion,AngleProtonMuon,AngleShr1Muon,AngleShr1Pion,AngleShr1Proton;
float DistMuonPion,DistProtonPion,DistProtonMuon;
float LambdaMass,DistShr1Muon,DistShr1Pion,DistShr1Proton;
float Weight;

static void add33(TMVA::Reader* R){
  R->AddVariable("MuonTrkLength",&MuonTrkLength);   R->AddVariable("MuonTrkDistVtx",&MuonTrkDistVtx);
  R->AddVariable("MuonTrkScore",&MuonTrkScore);     R->AddVariable("MuonNHitsY",&MuonNHitsY);
  R->AddVariable("MuonPhotonPID",&MuonPhotonPID);   R->AddVariable("PionPID",&PionPID);
  R->AddVariable("PionTrkLength",&PionTrkLength);   R->AddVariable("PionTrkDistVtx",&PionTrkDistVtx);
  R->AddVariable("PionTrkScore",&PionTrkScore);     R->AddVariable("PionNHitsY",&PionNHitsY);
  R->AddVariable("PionDeDxPlane2",&PionDeDxPlane2); R->AddVariable("PionPhotonPID",&PionPhotonPID);
  R->AddVariable("ProtonPID",&ProtonPID);           R->AddVariable("ProtonTrkLength",&ProtonTrkLength);
  R->AddVariable("ProtonTrkDistVtx",&ProtonTrkDistVtx); R->AddVariable("ProtonTrkScore",&ProtonTrkScore);
  R->AddVariable("ProtonNHitsY",&ProtonNHitsY);     R->AddVariable("ProtonDeDxPlane2",&ProtonDeDxPlane2);
  R->AddVariable("ProtonPhotonPID",&ProtonPhotonPID); R->AddVariable("Shr1NHitsY",&Shr1NHitsY);
  R->AddVariable("Shr1VtxDist",&Shr1VtxDist);       R->AddVariable("Shr1OpenAngle",&Shr1OpenAngle);
  R->AddVariable("Shr2NHitsY",&Shr2NHitsY);         R->AddVariable("Shr2VtxDist",&Shr2VtxDist);
  R->AddVariable("AngleMuonPion",&AngleMuonPion);   R->AddVariable("AngleProtonPion",&AngleProtonPion);
  R->AddVariable("AngleProtonMuon",&AngleProtonMuon); R->AddVariable("AngleShr1Muon",&AngleShr1Muon);
  R->AddVariable("AngleShr1Pion",&AngleShr1Pion);   R->AddVariable("AngleShr1Proton",&AngleShr1Proton);
  R->AddVariable("DistMuonPion",&DistMuonPion);     R->AddVariable("DistProtonPion",&DistProtonPion);
  R->AddVariable("DistProtonMuon",&DistProtonMuon);
}
static void bind(TTree* t){
  t->SetBranchAddress("MuonTrkLength",&MuonTrkLength);t->SetBranchAddress("MuonTrkDistVtx",&MuonTrkDistVtx);
  t->SetBranchAddress("MuonTrkScore",&MuonTrkScore);t->SetBranchAddress("MuonNHitsY",&MuonNHitsY);
  t->SetBranchAddress("MuonPhotonPID",&MuonPhotonPID);t->SetBranchAddress("PionPID",&PionPID);
  t->SetBranchAddress("PionTrkLength",&PionTrkLength);t->SetBranchAddress("PionTrkDistVtx",&PionTrkDistVtx);
  t->SetBranchAddress("PionTrkScore",&PionTrkScore);t->SetBranchAddress("PionNHitsY",&PionNHitsY);
  t->SetBranchAddress("PionDeDxPlane2",&PionDeDxPlane2);t->SetBranchAddress("PionPhotonPID",&PionPhotonPID);
  t->SetBranchAddress("ProtonPID",&ProtonPID);t->SetBranchAddress("ProtonTrkLength",&ProtonTrkLength);
  t->SetBranchAddress("ProtonTrkDistVtx",&ProtonTrkDistVtx);t->SetBranchAddress("ProtonTrkScore",&ProtonTrkScore);
  t->SetBranchAddress("ProtonNHitsY",&ProtonNHitsY);t->SetBranchAddress("ProtonDeDxPlane2",&ProtonDeDxPlane2);
  t->SetBranchAddress("ProtonPhotonPID",&ProtonPhotonPID);t->SetBranchAddress("Shr1NHitsY",&Shr1NHitsY);
  t->SetBranchAddress("Shr1VtxDist",&Shr1VtxDist);t->SetBranchAddress("Shr1OpenAngle",&Shr1OpenAngle);
  t->SetBranchAddress("Shr2NHitsY",&Shr2NHitsY);t->SetBranchAddress("Shr2VtxDist",&Shr2VtxDist);
  t->SetBranchAddress("AngleMuonPion",&AngleMuonPion);t->SetBranchAddress("AngleProtonPion",&AngleProtonPion);
  t->SetBranchAddress("AngleProtonMuon",&AngleProtonMuon);t->SetBranchAddress("AngleShr1Muon",&AngleShr1Muon);
  t->SetBranchAddress("AngleShr1Pion",&AngleShr1Pion);t->SetBranchAddress("AngleShr1Proton",&AngleShr1Proton);
  t->SetBranchAddress("DistMuonPion",&DistMuonPion);t->SetBranchAddress("DistProtonPion",&DistProtonPion);
  t->SetBranchAddress("DistProtonMuon",&DistProtonMuon);
  t->SetBranchAddress("LambdaMass",&LambdaMass);t->SetBranchAddress("DistShr1Muon",&DistShr1Muon);
  t->SetBranchAddress("DistShr1Pion",&DistShr1Pion);t->SetBranchAddress("DistShr1Proton",&DistShr1Proton);
  t->SetBranchAddress("Weight",&Weight);
}

void EvalCompare(){
  TFile f("TMVA_input_v3.root");
  TTree* tS=(TTree*)f.Get("TreeS"); TTree* tB=(TTree*)f.Get("TreeB");

  // production reader: 33 vars
  TMVA::Reader* Rp=new TMVA::Reader("!Color:!Silent"); add33(Rp);
  Rp->AddSpectator("Weight",&Weight);
  Rp->BookMVA("BDTg","SigmaDataset/weights/TMVASigma_BDTg.weights.xml");

  // v2 reader: 33 + LambdaMass + DistShr1{Mu,Pi,Pr}  (order must match RunTMVA_v2.C)
  TMVA::Reader* Rv=new TMVA::Reader("!Color:!Silent"); add33(Rv);
  Rv->AddVariable("LambdaMass",&LambdaMass);
  Rv->AddVariable("DistShr1Muon",&DistShr1Muon);
  Rv->AddVariable("DistShr1Pion",&DistShr1Pion);
  Rv->AddVariable("DistShr1Proton",&DistShr1Proton);
  Rv->AddSpectator("Weight",&Weight);
  Rv->BookMVA("BDTg","SigmaDataset_v2/weights/TMVASigma_BDTg.weights.xml");

  const int NB=2000;
  TH1D *pS=new TH1D("pS","",NB,-1,1),*pB=new TH1D("pB","",NB,-1,1);
  TH1D *vS=new TH1D("vS","",NB,-1,1),*vB=new TH1D("vB","",NB,-1,1);

  for(int pass=0;pass<2;pass++){
    TTree* t=pass?tB:tS; bind(t);
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      float bp=Rp->EvaluateMVA("BDTg"), bv=Rv->EvaluateMVA("BDTg");
      if(pass){pB->Fill(bp,Weight);vB->Fill(bv,Weight);}
      else    {pS->Fill(bp,Weight);vS->Fill(bv,Weight);}
    }
  }

  auto scan=[&](const char* tag,TH1D*hS,TH1D*hB){
    cout<<"\n== "<<tag<<" ==\n";
    cout<<setw(8)<<"thr"<<setw(10)<<"S"<<setw(10)<<"B"<<setw(10)<<"eff%"
        <<setw(11)<<"S/sqrtB"<<"\n";
    double thr[]={0.80,0.85,0.90,0.92,0.94,0.95,0.96,0.97,0.98,0.99};
    cout<<fixed<<setprecision(3);
    for(double tt:thr){int b0=hS->FindBin(tt+1e-6);
      double Ss=hS->Integral(b0,NB),Bb=hB->Integral(b0,NB);
      cout<<setw(8)<<tt<<setw(10)<<Ss<<setw(10)<<Bb
          <<setw(10)<<100*Ss/DENOM<<setw(11)<<(Bb>0?Ss/sqrt(Bb):0)<<"\n";}
  };
  cout<<fixed<<setprecision(3);
  cout<<"\nPre-selected totals:  S(TreeS)="<<pS->Integral()<<"  B(TreeB)="<<pB->Integral()<<"\n";
  scan("PRODUCTION BDTg (33 vars)",pS,pB);
  scan("v2 BDTg (33 + LambdaMass + DistShr1*)",vS,vB);

  // headline: at the production working-point background (B at thr=0.96), what eff does each give?
  int b96=pS->FindBin(0.96+1e-6);
  double Bref=pB->Integral(b96,NB), Sref=pS->Integral(b96,NB);
  // find v2 threshold giving B <= Bref
  double thrV=1.0,Sv=0;
  for(int b=NB;b>=1;b--){ double Bb=vB->Integral(b,NB); if(Bb>Bref){break;} thrV=vS->GetBinLowEdge(b); Sv=vS->Integral(b,NB);}
  cout<<"\n=== HEADLINE (matched background) ===\n";
  cout<<"  production @0.96 : S="<<Sref<<"  B="<<Bref<<"  eff="<<100*Sref/DENOM<<"%\n";
  cout<<"  v2 @thr="<<thrV<<" : S="<<Sv<<"  B<="<<Bref
      <<"  eff="<<100*Sv/DENOM<<"%   (gain "<<100*(Sv-Sref)/Sref<<"% signal)\n";
}
