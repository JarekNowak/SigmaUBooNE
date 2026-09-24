// ScanPeleeBDT.C -- locate the FoM peak of the PeLEE BDTG from BDTInput.root.
// BdtScore is already cross-validated (even->foldB, odd->foldA), so the scan is
// out-of-sample. Sample 0 is the Sigma0 signal; everything else is background.
void ScanPeleeBDT(const char* fn="BDTInput.root", const char* tag="new"){
  TFile f(fn); if(f.IsZombie()){printf("  cannot open %s\n",fn);return;}
  TTree* TS=(TTree*)f.Get("TreeS"); TTree* TB=(TTree*)f.Get("TreeB");
  std::vector<double> sS,bS,sW,bW;
  for(int isSig=0;isSig<2;++isSig){
    TTree* t=isSig?TS:TB; if(!t) continue;
    float w=0,sc=0,smp=0;   // Sample is stored as Float_t
    t->SetBranchAddress("Weight",&w); t->SetBranchAddress("BdtScore",&sc);
    t->SetBranchAddress("Sample",&smp);
    for(Long64_t e=0;e<t->GetEntries();++e){ t->GetEntry(e);
      if(isSig){ sS.push_back(sc); sW.push_back(w); }
      else     { bS.push_back(sc); bW.push_back(w); } }
    t->ResetBranchAddresses();
  }
  double S0=0,B0=0; for(double x:sW)S0+=x; for(double x:bW)B0+=x;
  auto at=[&](double c,double&s,double&b){ s=0;b=0;
    for(size_t i=0;i<sS.size();++i) if(sS[i]>c) s+=sW[i];
    for(size_t i=0;i<bS.size();++i) if(bS[i]>c) b+=bW[i]; };
  double best=-1,bc=-2;
  for(double c=-1.0;c<=1.0;c+=0.002){ double s,b; at(c,s,b);
    if(s<=0||b<=0) continue; double fm=s/std::sqrt(s+b); if(fm>best){best=fm;bc=c;} }
  double s90,b90; at(-0.90,s90,b90);
  printf("\n  === %s assignment (%s) ===\n",tag,fn);
  printf("  preselection      S=%.3f B=%.1f pur=%.4f%% FoM=%.4f\n",
         S0,B0,100*S0/(S0+B0),S0/std::sqrt(S0+B0));
  printf("  at BDTG > -0.90   S=%.3f B=%.1f pur=%.4f%% FoM=%.4f\n",
         s90,b90,100*s90/(s90+b90),s90/std::sqrt(s90+b90));
  double sp,bp; at(bc,sp,bp);
  printf("  FoM peak  > %+.3f S=%.3f B=%.1f pur=%.4f%% FoM=%.4f  (eff %.1f%% of presel)\n",
         bc,sp,bp,100*sp/(sp+bp),best,100*sp/S0);
}
