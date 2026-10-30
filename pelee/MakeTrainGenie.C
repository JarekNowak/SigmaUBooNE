// MakeTrainGenie.C -- training/evaluation inputs from BDTInput.root (TreeS NuWro, TreeB,
// TreeG GENIE Sigma0 QE signal, written by SigmaSelection.C since 2026-10-30).
// mode 1: TreeS = NuWro + GENIE, own weights.  mode 2: same, each generator scaled to equal total weight.
// mode 9: evaluation file with TreeS = GENIE (TreeG), TreeB = TreeB.
void MakeTrainGenie(const char* in, const char* out, int mode){
  TFile fi(in); TTree* ts=(TTree*)fi.Get("TreeS"); TTree* tb=(TTree*)fi.Get("TreeB"); TTree* tg=(TTree*)fi.Get("TreeG");
  auto sumW=[](TTree* t){ t->Draw("Weight","","goff"); double s=0; for(Long64_t i=0;i<t->GetSelectedRows();++i) s+=t->GetV1()[i]; return s; };
  const double ws=sumW(ts), wg=sumW(tg);
  TFile fo(out,"RECREATE");
  if (mode==9) { TTree* a=tg->CloneTree(-1,"fast"); a->Write("TreeS"); tb->CloneTree(-1,"fast")->Write("TreeB"); fo.Close(); return; }
  const double kg = (mode==2) ? ws/wg : 1.0;
  printf("NuWro sum %.3f (%lld)  GENIE sum %.3f (%lld)  GENIE scaled by %.4f\n",ws,ts->GetEntries(),wg,tg->GetEntries(),kg);
  TChain ch("TreeS"); ch.Add(Form("%s/TreeS",in)); ch.Add(Form("%s/TreeG",in));
  float w=0; ch.SetBranchAddress("Weight",&w);
  fo.cd(); TTree* o=ch.CloneTree(0); o->SetName("TreeS");
  for(Long64_t i=0;i<ch.GetEntries();++i){ ch.GetEntry(i); if (ch.GetTreeNumber()==1) w*=kg; o->Fill(); }
  o->Write("TreeS"); fo.cd(); tb->CloneTree(-1,"fast")->Write("TreeB");
  printf("wrote %s: TreeS %lld\n",out,o->GetEntries());
  fo.Close();
}
// Test of 2026-10-30 (README, "BDT training with GENIE signal"):
//   root -l -b -q 'MakeTrainGenie.C("BDTInput.root","in_m1.root",1)'
//   root -l -b -q 'MakeTrainGenie.C("BDTInput.root","in_m9.root",9)'
//   root -l -b -q 'TrainBDT.C("in_m1.root",0.05,"B_m1_s100",100)'
//   root -l -b -q 'RescoreBDT.C+("BDTInput.root","nuwro_m1.root","B_m1_s100_A","B_m1_s100_B")'  # NuWro
//   root -l -b -q 'RescoreBDT.C+("in_m9.root","genie_m1.root","B_m1_s100_A","B_m1_s100_B")'    # GENIE
//   root -l -b -q 'CompareBDT.C+("nuwro_m1.root","genie_m1.root")'
//   root -l -b -q 'Sensitivity.C+("nuwro_m1.root","SigmaHists.root",0.759,0.20)'   # NuWro limits
