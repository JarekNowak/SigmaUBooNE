// apply_split.C -- evaluate the per-class BDT weights on BDTInput.root and write the
// three scores per event to SplitScores.root, for eval_split.C to compare against the
// single BDT.
//
// The two-fold cross-validation is preserved exactly as in the selection: fold A is
// trained on even events and applied to odd, fold B the reverse, so no event is ever
// scored by a BDT that saw it in training. Getting this wrong inflates the ROC and
// the apparent gain, which is the whole reason the split is being tested.
//
//   usage: root -l -b -q 'apply_split.C'
#include <vector>
#include <string>
#include "TMVA/Reader.h"
R__LOAD_LIBRARY(libTMVA)

void apply_split(){
  const char* CLS[3] = {"pi0","nunopi0","cosmic"};

  TFile* fin=TFile::Open("BDTInput.root","READ");
  if(!fin||fin->IsZombie()){ printf("  cannot open BDTInput.root\n"); return; }

  // variable list, read from the training tree so it cannot drift from the trainer
  TTree* probe=(TTree*)fin->Get("TreeS");
  std::vector<std::string> vars;
  TObjArray* br=probe->GetListOfBranches();
  for(int i=0;i<br->GetEntries();++i){
    std::string n=br->At(i)->GetName();
    if(n=="Weight"||n=="IsOdd"||n=="BdtScore") continue;
    // truth/diagnostic branches are not BDT inputs
    static const std::vector<std::string> skip={"Sample","NuPdg","Ccnc","Interaction",
      "NuE","HasPi0","CosmicMu","MuBktPdg","ProtonBktPdg","PionBktPdg","Shr1BktPdg"};
    if(std::find(skip.begin(),skip.end(),n)!=skip.end()) continue;
    vars.push_back(n);
  }
  printf("  %zu BDT input variables\n", vars.size());

  std::vector<float> val(vars.size(), 0.f);
  float wt=0.f, isodd=0.f;
  TMVA::Reader* rdr[3][2] = {{nullptr,nullptr},{nullptr,nullptr},{nullptr,nullptr}};
  for(int c=0;c<3;++c) for(int fo=0;fo<2;++fo){
    TString w=TString::Format("SigmaBDT_%s_%c/weights/TMVASigma_BDTG.weights.xml",
                              CLS[c], fo?'B':'A');
    if(gSystem->AccessPathName(w)!=kFALSE){ printf("  missing %s\n", w.Data()); continue; }
    rdr[c][fo]=new TMVA::Reader("!Color:!Silent");
    for(size_t k=0;k<vars.size();++k) rdr[c][fo]->AddVariable(vars[k].c_str(), &val[k]);
    rdr[c][fo]->AddSpectator("Weight",&wt);
    rdr[c][fo]->AddSpectator("IsOdd",&isodd);
    rdr[c][fo]->BookMVA("BDTG", w);
  }

  TFile* fout=new TFile("SplitScores.root","RECREATE");
  for(const char* tn : {"TreeS","TreeB"}){
    TTree* t=(TTree*)fin->Get(tn);
    if(!t) continue;
    for(size_t k=0;k<vars.size();++k) t->SetBranchAddress(vars[k].c_str(), &val[k]);
    t->SetBranchAddress("Weight",&wt);
    t->SetBranchAddress("IsOdd",&isodd);
    fout->cd();
    TTree* ot=new TTree(strcmp(tn,"TreeS")?"ScoreB":"ScoreS","per-class BDT scores");
    float sc[3]; float w_out=0;
    for(int c=0;c<3;++c) ot->Branch(Form("Score_%s",CLS[c]), &sc[c]);
    ot->Branch("Weight",&w_out);
    Long64_t N=t->GetEntries();
    for(Long64_t i=0;i<N;++i){
      t->GetEntry(i);
      // odd events -> fold A (trained on even); even -> fold B. Never self-scored.
      int fo = (isodd>0.5f) ? 0 : 1;
      for(int c=0;c<3;++c)
        sc[c] = rdr[c][fo] ? rdr[c][fo]->EvaluateMVA("BDTG") : -99.f;
      w_out=wt;
      ot->Fill();
    }
    ot->Write();
    printf("  %s -> %lld events scored\n", tn, N);
    t->ResetBranchAddresses();
  }
  fout->Close(); fin->Close();
  printf("  wrote SplitScores.root\n");
}
