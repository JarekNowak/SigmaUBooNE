// ClampBDTInputs.C -- copy TreeS/TreeB with each listed variable clamped to [L, H]. L/H are the
// 0.1%/99.9% quantiles of its values inside the bulk window [q1 - 10 w, q99 + 10 w],
// w = q99 - q1, of the values with |x| < 1e30 (signal and background together), so a
// sentinel or garbage tail holding more than 0.1% of events does not set the range.
// Values below L (including the -FLT_MAX "missing" sentinel) go to L - 1% of (H - L) and
// values above H to H + 1%, so both stay distinct values. Monotonic, so only TMVA's cut
// grid (equidistant between a node's min and max) can see the difference.
// capsCSV ("var=H,...") fixes H by hand where more than 1% of the values are garbage, so
// that even q99 sits in the tail (calorimetric energies, hypana shower energy).
//
// Test of 2026-10-19 (README, "Not adopted: input ranges"), PeLEE tree:
//   V=$(sed -n '/kVars\[\] = {/,/^};/p' TrainBDT.C | grep -v '//' | grep -o '"[A-Za-z0-9_]*"' | tr -d '"' | paste -sd,)
//   root -l -b -q "ClampBDTInputs.C+(\"BDTInput.root\",\"BDTInput_clamp.root\",\"$V\",\"PionCaloEnergyY=2000,ProtonCaloEnergyY=2000\")"
//   root -l -b -q 'TrainBDT.C("BDTInput_clamp.root", 0.05, "SigmaBDT_clamp")'
//   root -l -b -q 'RescoreBDT.C+("BDTInput_clamp.root","scored_clamp.root","SigmaBDT_clamp_A","SigmaBDT_clamp_B")'
//   root -l -b -q 'CompareBDT.C+("BDTInput.root","scored_clamp.root")'
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjArray.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
void ClampBDTInputs(const char* in, const char* out, const char* varsCSV, const char* capsCSV=""){
  std::map<TString,double> caps; { TObjArray* a = TString(capsCSV).Tokenize(",");
    for (auto* o : *a) { TString kv=o->GetName(); const int e=kv.Index("=");
      caps[kv(0,e)]=TString(kv(e+1,kv.Length())).Atof(); } delete a; }
  std::vector<TString> vars; { TObjArray* a = TString(varsCSV).Tokenize(",");
    for (auto* o : *a) vars.push_back(o->GetName()); delete a; }
  TFile fi(in); TTree* ts=(TTree*)fi.Get("TreeS"); TTree* tb=(TTree*)fi.Get("TreeB");
  std::map<TString,std::pair<double,double>> rng;
  for (auto& v : vars) { std::vector<double> x;
    for (TTree* t : {ts,tb}) { t->Draw(v,"","goff"); for (Long64_t i=0;i<t->GetSelectedRows();++i)
      if (std::fabs(t->GetV1()[i]) < 1e30) x.push_back(t->GetV1()[i]); }
    std::sort(x.begin(),x.end());
    const double q1=x[(size_t)(0.01*(x.size()-1))], q99=x[(size_t)(0.99*(x.size()-1))], w=q99-q1;
    std::vector<double> y; for (double z : x) if (z>=q1-10*w && z<=q99+10*w) y.push_back(z);
    const double L=y[(size_t)(0.001*(y.size()-1))];
    const double H=caps.count(v)?caps[v]:y[(size_t)(0.999*(y.size()-1))];
    rng[v]={L,H}; }
  TFile fo(out,"RECREATE");
  for (TTree* t : {ts,tb}) {
    std::vector<float> buf(vars.size());
    for (size_t i=0;i<vars.size();++i) t->SetBranchAddress(vars[i],&buf[i]);
    fo.cd(); TTree* o=t->CloneTree(0);
    std::vector<Long64_t> nlo(vars.size(),0), nhi(vars.size(),0);
    for (Long64_t e=0;e<t->GetEntries();++e){ t->GetEntry(e);
      for (size_t i=0;i<vars.size();++i){ const double L=rng[vars[i]].first, H=rng[vars[i]].second;
        if (!(buf[i]>=L)) { buf[i]=(float)(L-0.01*(H-L)); nlo[i]++; }   // also catches NaN
        else if (buf[i]>H) { buf[i]=(float)(H+0.01*(H-L)); nhi[i]++; } }
      o->Fill(); }
    o->Write();
    printf("  %s: %lld events\n", t->GetName(), t->GetEntries());
    for (size_t i=0;i<vars.size();++i) if (nlo[i]+nhi[i] > 0.003*t->GetEntries())
      printf("    %-20s [%-10.4g, %-10.4g]  low %lld  high %lld\n", vars[i].Data(), rng[vars[i]].first, rng[vars[i]].second, nlo[i], nhi[i]);
    t->ResetBranchAddresses();
  }
  fo.Close();
}
