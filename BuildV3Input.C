// ============================================================
//  BuildV3Input.C
//  Clone TreeS/TreeB from TMVA_input.root into TMVA_input_v3.root,
//  adding one new branch: LambdaMass = M(p pi-) reconstructed from
//  range-based momenta + the stored proton-pion opening angle.
//  (Everything else is copied unchanged.)
//
//  Usage:  root -b -l -q BuildV3Input.C
// ============================================================
#include "TFile.h"
#include "TTree.h"
#include <cmath>
#include <iostream>

static double pFromRange(double R,double m,double a,double b){
  if(R<=0) return 0; double KE=pow(R/a,1.0/b); return sqrt(KE*KE+2.0*KE*m);
}
static double Em(double p,double m){ return sqrt(p*p+m*m); }

static void clone(TFile* in, TFile* out, const char* name){
  const double mp=938.272, mpi=139.570;
  TTree* t=(TTree*)in->Get(name);
  out->cd();
  TTree* nt=t->CloneTree(0);          // empty clone, same branches
  float prL,piL,ang, LambdaMass;
  t->SetBranchAddress("ProtonTrkLength",&prL);
  t->SetBranchAddress("PionTrkLength",&piL);
  t->SetBranchAddress("AngleProtonPion",&ang);
  TBranch* bM=nt->Branch("LambdaMass",&LambdaMass,"LambdaMass/F");
  for(Long64_t i=0;i<t->GetEntries();i++){
    t->GetEntry(i);
    if(prL>0&&piL>0&&ang>-998){
      double pp=pFromRange(prL,mp,8.02e-4,1.82);
      double ppi=pFromRange(piL,mpi,1.134e-2,1.586);
      double m2=mp*mp+mpi*mpi+2.0*(Em(pp,mp)*Em(ppi,mpi)-pp*ppi*cos(ang));
      LambdaMass=(m2>0)?(float)sqrt(m2):-999.f;
    } else LambdaMass=-999.f;
    nt->Fill();
  }
  nt->Write();
  std::cout<<name<<" -> "<<nt->GetEntries()<<" entries, +LambdaMass\n";
}

void BuildV3Input(){
  TFile* in =TFile::Open("TMVA_input.root","READ");
  TFile* out=TFile::Open("TMVA_input_v3.root","RECREATE");
  clone(in,out,"TreeS");
  clone(in,out,"TreeB");
  out->Close(); in->Close();
  std::cout<<"Wrote TMVA_input_v3.root\n";
}
