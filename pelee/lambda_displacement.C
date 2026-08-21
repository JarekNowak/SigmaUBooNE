// lambda_displacement.C -- would a Lambda displaced-vertex requirement work as a
// PRESELECTION cut, before any BDT?
//
// Lambda has ctau = 7.89 cm, so a genuine Lambda -> p pi- decays measurably away from
// the primary vertex. Combinatorial background built from two primary tracks should
// not. LambdaDecayDist (the p/pi common origin, taken as the midpoint of the two
// track starts, measured from the reco vertex) and LambdaPointAngle (the angle
// between the reconstructed Lambda momentum and its flight direction) are already
// computed but only ever used as BDT inputs.
void lambda_displacement(){
  TFile* f=TFile::Open("BDTInput.root");
  TTree* S=(TTree*)f->Get("TreeS"); TTree* B=(TTree*)f->Get("TreeB");
  double ns=S->GetEntries(), nb=B->GetEntries();
  double base=ns/std::sqrt(ns+nb);
  printf("\n  preselection: S=%.0f  B=%.0f  FoM=%.4f\n", ns, nb, base);

  printf("\n  --- LambdaDecayDist scan ---\n");
  printf("  %8s %8s %9s %8s %8s %8s\n","cut[cm]","S","B","eff[%]","rej[%]","FoM x");
  for(double c : {1.,2.,3.,4.,5.,7.,10.,15.,20.,30.}){
    TString cut=TString::Format("LambdaDecayDist>%g",c);
    double s=S->GetEntries(cut), b=B->GetEntries(cut);
    if(s+b<=0) continue;
    printf("  %8.0f %8.0f %9.0f %8.1f %8.1f %8.3f\n",
           c,s,b,100*s/ns,100*(1-b/nb),(s/std::sqrt(s+b))/base);
  }
  printf("\n  --- LambdaPointAngle scan (collinearity of p_Lambda with its flight) ---\n");
  printf("  %8s %8s %9s %8s %8s %8s\n","cut[rad]","S","B","eff[%]","rej[%]","FoM x");
  for(double c : {0.2,0.3,0.4,0.5,0.7,1.0,1.5}){
    TString cut=TString::Format("LambdaPointAngle>=0 && LambdaPointAngle<%g",c);
    double s=S->GetEntries(cut), b=B->GetEntries(cut);
    if(s+b<=0) continue;
    printf("  %8.2f %8.0f %9.0f %8.1f %8.1f %8.3f\n",
           c,s,b,100*s/ns,100*(1-b/nb),(s/std::sqrt(s+b))/base);
  }
  printf("\n  --- combined, and against the pi0 class specifically ---\n");
  const char* combos[] = {
    "LambdaDecayDist>3",
    "LambdaDecayDist>5",
    "LambdaDecayDist>3 && LambdaPointAngle<0.5",
    "LambdaDecayDist>5 && LambdaPointAngle<0.5",
    "LambdaDecayDist>2 && LambdaPointAngle<0.7" };
  printf("  %-42s %7s %9s %7s %7s %7s\n","cut","S","B","eff[%]","rej[%]","FoM x");
  for(auto c:combos){
    double s=S->GetEntries(c), b=B->GetEntries(c);
    if(s+b<=0) continue;
    printf("  %-42s %7.0f %9.0f %7.1f %7.1f %7.3f\n",
           c,s,b,100*s/ns,100*(1-b/nb),(s/std::sqrt(s+b))/base);
  }
  printf("\n  mean LambdaDecayDist:  signal %.2f cm   background %.2f cm   (ctau = 7.89 cm)\n",
    (S->Draw("LambdaDecayDist>>hs(200,0,100)","LambdaDecayDist>=0","goff"), ((TH1D*)gDirectory->Get("hs"))->GetMean()),
    (B->Draw("LambdaDecayDist>>hb(200,0,100)","LambdaDecayDist>=0","goff"), ((TH1D*)gDirectory->Get("hb"))->GetMean()));
}
