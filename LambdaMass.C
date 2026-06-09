// ============================================================
//  LambdaMass.C
//  Reconstruct the Lambda -> p pi- invariant mass from quantities
//  ALREADY stored in TMVA_input.root:
//     proton momentum  <- range (ProtonTrkLength)
//     pion momentum    <- range (PionTrkLength)
//     opening angle    <- AngleProtonPion
//  and test its signal/background separation.
//
//  (The Sigma0 mass needs the photon energy, but shr_energy_plane* is
//   100% sentinel in this ntuple, so it cannot be reconstructed here.)
//
//  Range->momentum: power-law fits to PSTAR/muon range tables in LAr.
//     proton: KE = (R/8.02e-4)^(1/1.82)        [R in cm, KE in MeV]
//     pion  : KE = (R/1.134e-2)^(1/1.586)      (muon-like range approx)
//
//  Usage:  root -b -l -q LambdaMass.C
// ============================================================
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include <cmath>
#include <iostream>

static double pFromRange(double R,double m,double a,double b){
  if(R<=0) return 0;
  double KE=pow(R/a,1.0/b);            // MeV
  return sqrt(KE*KE+2.0*KE*m);         // MeV/c
}
static double Emass(double p,double m){ return sqrt(p*p+m*m); }

void LambdaMass()
{
  gStyle->SetOptStat(0);
  const double mp=938.272, mpi=139.570, mLam=1115.683;

  TFile f("TMVA_input.root");
  TTree*tS=(TTree*)f.Get("TreeS");
  TTree*tB=(TTree*)f.Get("TreeB");

  float prL,piL,ang,w;
  auto setup=[&](TTree*t){
    t->SetBranchAddress("ProtonTrkLength",&prL);
    t->SetBranchAddress("PionTrkLength",&piL);
    t->SetBranchAddress("AngleProtonPion",&ang);
    t->SetBranchAddress("Weight",&w);
  };

  TH1D*hS=new TH1D("hS","#Lambda#rightarrowp#pi^{-} invariant mass;M(p#pi^{-}) [MeV];norm.",60,1000,1600);
  TH1D*hB=new TH1D("hB","",60,1000,1600);
  TH1D*hSraw=new TH1D("hSraw","",60,1000,1600); // unweighted for stats

  for(int pass=0;pass<2;pass++){
    TTree*t=pass?tB:tS; setup(t);
    for(Long64_t i=0;i<t->GetEntries();i++){
      t->GetEntry(i);
      if(prL<=0||piL<=0||ang<-998) continue;
      double pp=pFromRange(prL,mp,8.02e-4,1.82);
      double ppi=pFromRange(piL,mpi,1.134e-2,1.586);
      double E=Emass(pp,mp)+Emass(ppi,mpi);
      double m2=mp*mp+mpi*mpi+2.0*(Emass(pp,mp)*Emass(ppi,mpi)-pp*ppi*cos(ang));
      double M=(m2>0)?sqrt(m2):0;
      if(pass){hB->Fill(M,w);} else {hS->Fill(M,w); hSraw->Fill(M);}
    }
  }

  // TMVA-style separation <s^2> on area-normalised shapes
  TH1D*nS=(TH1D*)hS->Clone("nS"); TH1D*nB=(TH1D*)hB->Clone("nB");
  if(nS->Integral()>0)nS->Scale(1./nS->Integral());
  if(nB->Integral()>0)nB->Scale(1./nB->Integral());
  double sep=0;
  for(int b=1;b<=nS->GetNbinsX();b++){
    double s=nS->GetBinContent(b),bb=nB->GetBinContent(b);
    if(s+bb>0) sep+=0.5*(s-bb)*(s-bb)/(s+bb);
  }
  // signal window efficiency/purity around PDG mass +-50 MeV
  int b1=hS->FindBin(mLam-50), b2=hS->FindBin(mLam+50);
  double Sin=hS->Integral(b1,b2), Stot=hS->Integral();
  double Bin=hB->Integral(b1,b2), Btot=hB->Integral();

  std::cout.setf(std::ios::fixed); std::cout.precision(3);
  std::cout<<"\n=== Lambda invariant mass M(p pi-) ===\n";
  std::cout<<"  signal  : mean="<<hS->GetMean()<<"  RMS="<<hS->GetRMS()
           <<"  (PDG="<<mLam<<")  raw entries="<<hSraw->GetEntries()<<"\n";
  std::cout<<"  backgnd : mean="<<hB->GetMean()<<"  RMS="<<hB->GetRMS()<<"\n";
  std::cout<<"  TMVA separation <s^2> = "<<sep<<"   (1=perfect, 0=none)\n";
  std::cout<<"  window "<<mLam-50<<"-"<<mLam+50<<" MeV : "
           <<"signal keeps "<<100*Sin/Stot<<"% , bkg keeps "<<100*Bin/Btot<<"%\n";

  TCanvas*c=new TCanvas("c","LambdaMass",800,550);
  nB->SetLineColor(kRed+1); nB->SetFillColorAlpha(kRed+1,0.25); nB->SetLineWidth(2);
  nS->SetLineColor(kBlue+1);nS->SetFillColorAlpha(kBlue+1,0.25);nS->SetLineWidth(2);
  nS->SetTitle("#Lambda#rightarrowp#pi^{-} invariant mass;M(p#pi^{-}) [MeV];area-normalised");
  nS->SetMaximum(1.3*std::max(nS->GetMaximum(),nB->GetMaximum()));
  nS->Draw("HIST"); nB->Draw("HIST SAME");
  TLine*l=new TLine(mLam,0,mLam,nS->GetMaximum()); l->SetLineStyle(2);l->SetLineWidth(2);l->Draw();
  TLegend*lg=new TLegend(0.62,0.72,0.88,0.88);lg->SetBorderSize(0);
  lg->AddEntry(nS,"Signal","f");lg->AddEntry(nB,"Background","f");
  lg->AddEntry(l,"#Lambda PDG mass","l");lg->Draw();
  c->SaveAs("figs/EfficiencyPlots/LambdaMass.png");
  std::cout<<"\nSaved figs/EfficiencyPlots/LambdaMass.png\n";
}
