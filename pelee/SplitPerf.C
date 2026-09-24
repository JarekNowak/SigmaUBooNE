// How well does the BDT do on the paired vs pairless populations separately?
void SplitPerf(){
  TFile f("BDTInput.root"); TTree*TS=(TTree*)f.Get("TreeS"); TTree*TB=(TTree*)f.Get("TreeB");
  for(int mode=0;mode<3;++mode){
    const char* nm = mode==0?"all":(mode==1?"paired (HasPair=1)":"pairless (HasPair=0)");
    std::vector<double> sS,bS,sW,bW;
    for(int isSig=0;isSig<2;++isSig){
      TTree*t=isSig?TS:TB; float w=0,sc=0,hp=0;
      t->SetBranchAddress("Weight",&w); t->SetBranchAddress("BdtScore",&sc);
      t->SetBranchAddress("HasPair",&hp);
      for(Long64_t e=0;e<t->GetEntries();++e){ t->GetEntry(e);
        if(mode==1&&hp<0.5) continue; if(mode==2&&hp>0.5) continue;
        if(isSig){sS.push_back(sc);sW.push_back(w);} else {bS.push_back(sc);bW.push_back(w);} }
      t->ResetBranchAddresses(); }
    if(sS.size()<10||bS.size()<10){ printf("  %-22s too few\n",nm); continue; }
    // AUC
    std::vector<double> bs=bS; std::sort(bs.begin(),bs.end()); double tot=0;
    for(double x:sS) tot += (std::lower_bound(bs.begin(),bs.end(),x)-bs.begin())
        + 0.5*(std::upper_bound(bs.begin(),bs.end(),x)-std::lower_bound(bs.begin(),bs.end(),x));
    double auc=tot/(double(sS.size())*bs.size()); if(auc<0.5) auc=1-auc;
    double S0=0,B0=0; for(double w:sW)S0+=w; for(double w:bW)B0+=w;
    auto at=[&](double c,double&s,double&b){ s=0;b=0;
      for(size_t i=0;i<sS.size();++i) if(sS[i]>c) s+=sW[i];
      for(size_t i=0;i<bS.size();++i) if(bS[i]>c) b+=bW[i]; };
    double best=-1,bc=0; for(double c=-1;c<=1;c+=0.005){ double s,b; at(c,s,b);
      if(s<=0||b<=0) continue; double fm=s/std::sqrt(s+b); if(fm>best){best=fm;bc=c;} }
    double sp,bp; at(bc,sp,bp);
    printf("  %-22s AUC=%.3f | presel S=%6.3f B=%8.1f FoM=%.4f | peak %+.2f S=%5.3f B=%6.1f FoM=%.4f\n",
           nm,auc,S0,B0,S0/std::sqrt(S0+B0),bc,sp,bp,best);
  }
}
