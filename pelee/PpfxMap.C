// PpfxMap.C -- mean ppfx_cv against true E_nu for numubar CC in the true FV over all RHC
// overlay runs (ppfx_numubar_map.root), for SigmaSelection.C to weight samples that have
// no flux-parent record (the GENIE Sigma0 QE sample). Also prints the map-predicted and
// the measured <ppfx_cv> on the overlay's own Sigma0 signal.
// Run:  root -l -b -q 'PpfxMap.C("ppfx_numubar_map.root")'
void PpfxMap(const char* out){
  TChain ch("nuselection/NeutrinoSelectionFilter");
  const char* N="/data/uboone/new_numi_flux/";
  for (const char* f : {"Run1_rhc_new_numi_flux_rhc_pandora_ntuple.root","Run2_rhc_new_numi_flux_rhc_pandora_ntuple.root",
       "Run4a_rhc_new_numi_flux_rhc_pandora_ntuple.root","Run4b_rhc_new_numi_flux_rhc_pandora_ntuple.root","Run4c_rhc_new_numi_flux_rhc_pandora_ntuple.root",
       "Run3_rhc_new_numi_flux_rhc_pandora_ntuple_aa.root","Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ab.root","Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ac.root",
       "Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ad.root","Run3_rhc_new_numi_flux_rhc_pandora_ntuple_ae.root"}) ch.Add(Form("%s%s",N,f));
  const char* fv="true_nu_vtx_x>3&&true_nu_vtx_x<253.35&&true_nu_vtx_y>-112.53&&true_nu_vtx_y<114.47&&true_nu_vtx_z>3.1&&true_nu_vtx_z<1033&&!(true_nu_vtx_z>675.1&&true_nu_vtx_z<775.1)";
  const char* ok="ppfx_cv>0&&ppfx_cv<30";
  TFile fo(out,"RECREATE");
  TProfile* p=new TProfile("ppfx_numubar_cc",";true E_{#nu} [GeV];<ppfx_cv>",40,0,10);
  ch.Draw("ppfx_cv:nu_e>>ppfx_numubar_cc",Form("nu_pdg==-14&&ccnc==0&&%s&&%s",fv,ok),"goff prof");
  ch.Draw("ppfx_cv",Form("nu_pdg==-14&&ccnc==0&&Sum$(mc_pdg==3212)>0&&%s&&%s",fv,ok),"goff");
  double s=0; for(Long64_t i=0;i<ch.GetSelectedRows();i++) s+=ch.GetV1()[i];
  printf("signal: n %lld <ppfx_cv> %.4f\n",ch.GetSelectedRows(),s/ch.GetSelectedRows());
  ch.Draw("ppfx_cv:nu_e",Form("nu_pdg==-14&&ccnc==0&&Sum$(mc_pdg==3212)>0&&%s&&%s",fv,ok),"goff");
  double pr=0; for(Long64_t i=0;i<ch.GetSelectedRows();i++) pr+=p->GetBinContent(p->FindBin(ch.GetV2()[i]));
  printf("signal: map-predicted <ppfx_cv> %.4f\n",pr/ch.GetSelectedRows());
  for(int b=1;b<=20;b++) printf("  E %.2f-%.2f  <w> %.4f +- %.4f (%g)\n",p->GetXaxis()->GetBinLowEdge(b),p->GetXaxis()->GetBinUpEdge(b),p->GetBinContent(b),p->GetBinError(b),p->GetBinEntries(b));
  p->Write(); fo.Close();
}
