// eval_split.C -- does the per-class BDT beat the single one?
//
// The comparison must be at MATCHED SIGNAL EFFICIENCY, not at a fixed score, because
// the two discriminants have different score distributions and a common numerical
// threshold means nothing. For each single-BDT operating point, find the per-class
// thresholds giving the same signal efficiency and compare the background left.
//
// The per-class scores are combined by independent cuts, as gLEE does, rather than by
// a second-stage discriminant: with ~1000 signal events a stacked classifier would
// overtrain, and independent cuts can each be validated in their own sideband.
//
//   usage: root -l -b -q 'eval_split.C'
#include <vector>
#include <string>
#include <algorithm>

void eval_split(){
  TFile* f=TFile::Open("BDTInput.root");
  TTree* S=(TTree*)f->Get("TreeS"); TTree* B=(TTree*)f->Get("TreeB");
  if(!S||!B){ printf("  trees missing\n"); return; }
  if(!S->GetBranch("BdtScore")){ printf("  BdtScore missing -- run pass 2 first\n"); return; }
  const double ns=S->GetEntries(), nb=B->GetEntries();
  printf("\n  preselection  S=%.0f  B=%.0f\n", ns, nb);

  // --- reference: the single BDT ---
  printf("\n  single BDT, for reference\n");
  printf("  %10s %8s %10s %8s %9s\n","cut","S","B","eff[%]","FoM");
  std::vector<double> refEff, refFoM, refB;
  for(double c : {-0.95,-0.90,-0.80,-0.60,-0.40,-0.20,0.0,0.2}){
    TString q=TString::Format("BdtScore>%g",c);
    double s=S->GetEntries(q), b=B->GetEntries(q);
    if(s<=0) continue;
    double fom=s/std::sqrt(s+b);
    printf("  %10.2f %8.0f %10.0f %8.1f %9.4f\n", c,s,b,100*s/ns,fom);
    refEff.push_back(s/ns); refFoM.push_back(fom); refB.push_back(b);
  }

  // --- per-class scores, if the friend trees are present ---
  // The evaluation tree is written by apply_split.C; if absent, say so plainly rather
  // than silently reporting nothing.
  TFile* fs=TFile::Open("SplitScores.root","READ");
  if(!fs||fs->IsZombie()){
    printf("\n  SplitScores.root not found.\n"
           "  Run apply_split.C after training to evaluate the per-class weights on\n"
           "  BDTInput.root, then re-run this macro.\n");
    return;
  }
  TTree* TS=(TTree*)fs->Get("ScoreS"); TTree* TB=(TTree*)fs->Get("ScoreB");
  if(!TS||!TB){ printf("  score trees missing\n"); return; }

  // Scan a common quantile threshold on each class score, matched to each reference
  // signal efficiency, and report the background surviving all three cuts.
  const char* cls[3]={"pi0","nunopi0","cosmic"};
  printf("\n  per-class BDTs, thresholds matched to the single-BDT signal efficiency\n");
  printf("  %8s | %10s %10s | %s\n","eff[%]","B single","B split","gain in FoM");
  for(size_t k=0;k<refEff.size();++k){
    double target=refEff[k];
    // find, for each class, the threshold keeping the target efficiency; then apply
    // all three together. Because the cuts are correlated the joint efficiency is
    // lower than the target, so tighten iteratively on a common quantile.
    double lo=0., hi=1., qq=0.5;
    for(int it=0; it<40; ++it){
      qq=0.5*(lo+hi);
      TString cut;
      for(int c=0;c<3;++c){
        // quantile qq of the signal score for this class
        TS->Draw(Form("Score_%s>>hq(2000,-1.2,1.2)",cls[c]),"","goff");
        TH1D* h=(TH1D*)gDirectory->Get("hq");
        double pr[1]={1.0-qq}, out[1]; h->GetQuantiles(1,out,pr);
        cut += TString::Format("%sScore_%s>%g", c?" && ":"", cls[c], out[0]);
      }
      double s=TS->GetEntries(cut);
      if(s/ns > target) hi=qq; else lo=qq;
    }
    TString cut;
    for(int c=0;c<3;++c){
      TS->Draw(Form("Score_%s>>hq2(2000,-1.2,1.2)",cls[c]),"","goff");
      TH1D* h=(TH1D*)gDirectory->Get("hq2");
      double pr[1]={1.0-qq}, out[1]; h->GetQuantiles(1,out,pr);
      cut += TString::Format("%sScore_%s>%g", c?" && ":"", cls[c], out[0]);
    }
    double s=TS->GetEntries(cut), b=TB->GetEntries(cut);
    if(s<=0) continue;
    double fom=s/std::sqrt(s+b);
    printf("  %8.1f | %10.0f %10.0f | x%.3f\n", 100*s/ns, refB[k], b, fom/refFoM[k]);
  }
}
