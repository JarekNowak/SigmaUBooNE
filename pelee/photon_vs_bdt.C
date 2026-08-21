// Shower1Energy is already one of the 49 BDT inputs, so the question is not whether
// a soft-photon cut separates, but whether it adds anything the BDT has not already
// taken. Compare the FoM at the BDT operating point with and without the cut.
void photon_vs_bdt(){
  TFile f("BDTInput.root");
  TTree* S=(TTree*)f.Get("TreeS"); TTree* B=(TTree*)f.Get("TreeB");
  const double BDTCUT=-0.90;
  auto fom=[&](const char* cut)->std::vector<double>{
    double s=S->GetEntries(cut), b=B->GetEntries(cut);
    return {s,b,(s+b>0)?s/std::sqrt(s+b):0};
  };
  auto row=[&](const char* cut,const char* nm){
    auto r=fom(cut);
    printf("  %-44s  S %6.0f  B %7.0f  FoM %.4f\n", nm, r[0], r[1], r[2]);
    return r[2];
  };
  printf("\n  at preselection (no BDT)\n");
  double p0=row("1","no cut");
  row("Shower1Energy<0.12","E_gamma < 120 MeV");
  printf("\n  at the BDT operating point (BdtScore > -0.90)\n");
  double b0=row("BdtScore>-0.90","BDT only");
  double b1=row("BdtScore>-0.90 && Shower1Energy<0.12","BDT + E_gamma < 120 MeV");
  double b2=row("BdtScore>-0.90 && Shower1Energy<0.10","BDT + E_gamma < 100 MeV");
  double b3=row("BdtScore>-0.90 && Shower1Energy>0.005 && Shower1Energy<0.10","BDT + E in [5,100] MeV");
  printf("\n  gain of the cut ON TOP OF the BDT: x%.3f (120 MeV), x%.3f (100 MeV), x%.3f (window)\n",
    b1/b0, b2/b0, b3/b0);
  printf("  (preselection-only gain was x%.3f, so the BDT has already taken %s of it)\n",
    1.092, (b1/b0 < 1.02) ? "essentially all" : "only part");
}
