// fv_scan.C -- what does the fiducial volume actually buy, and would a different box
// do better?
//
// The current box (SigmaFV.h) is x [3.00, 253.35], y [-112.53, 114.47],
// z [3.1, 1033.0] excluding the z dead region (675.1, 775.1), inherited verbatim from
// the archived Sigma.C. It has never been optimised for this selection. The active
// TPC is x [0,256], y [-116,116], z [0,1036], so the current box is barely fiducial
// at all -- it trims only ~3 cm in x and y.
//
// Scanned on the reco vertex, at the AllEvents stage, so the result is the FV cut in
// isolation rather than folded through later cuts.
#include <vector>
void fv_scan(){
  const char* SIG="/data/uboone/prod_nuwro_sigmazero_overlay_numi_rhc_run3b_extunbiased_mcc9/ana/Hyperons_Nuwro_run3b_WithCR_pelee_2500.root";
  // one RHC overlay run stands in for the background shape; absolute normalisation
  // cancels in the figure of merit ratio used below
  const char* BKG="/data/uboone/new_numi_flux/Run1_rhc_new_numi_flux_rhc_pandora_ntuple.root";

  TFile* fs=TFile::Open(SIG); TTree* S=(TTree*)fs->Get("nuselection/NeutrinoSelectionFilter");
  TFile* fb=TFile::Open(BKG); TTree* B=(TTree*)fb->Get("nuselection/NeutrinoSelectionFilter");

  // signal must also be true-signal in the FV sense the analysis uses
  const char* SIGSEL = "nu_pdg==-14 && ccnc==0 && Sum$(mc_pdg==3212)>0";

  auto box=[&](double x0,double x1,double y0,double y1,double z0,double z1,bool dead){
    TString c=TString::Format(
      "reco_nu_vtx_sce_x>%g && reco_nu_vtx_sce_x<%g && "
      "reco_nu_vtx_sce_y>%g && reco_nu_vtx_sce_y<%g && "
      "reco_nu_vtx_sce_z>%g && reco_nu_vtx_sce_z<%g", x0,x1,y0,y1,z0,z1);
    if(dead) c += " && !(reco_nu_vtx_sce_z>675.1 && reco_nu_vtx_sce_z<775.1)";
    return c;
  };
  double ns0=S->GetEntries(SIGSEL), nb0=B->GetEntries();
  printf("\n  before any FV cut:  signal %.0f   background %.0f\n", ns0, nb0);

  struct C { const char* nm; double x0,x1,y0,y1,z0,z1; bool dead; };
  std::vector<C> cands = {
    {"current (SigmaFV.h)",      3.00,253.35,-112.53,114.47,  3.1,1033.0, true },
    {"current, no dead-region",  3.00,253.35,-112.53,114.47,  3.1,1033.0, false},
    {"5 cm inset",               5.00,251.00,-111.00,111.00,  5.0,1031.0, true },
    {"10 cm inset",             10.00,246.00,-106.00,106.00, 10.0,1026.0, true },
    {"20 cm inset",             20.00,236.00, -96.00, 96.00, 20.0,1016.0, true },
    {"30 cm inset",             30.00,226.00, -86.00, 86.00, 30.0,1006.0, true },
    {"CCpi analysis box",       10.00,246.00,-101.00,101.00, 10.0, 986.0, true },
    {"tight z only",             3.00,253.35,-112.53,114.47, 30.0,1000.0, true },
    {"tight x only",            20.00,236.00,-112.53,114.47,  3.1,1033.0, true },
    {"tight y only",             3.00,253.35, -96.00, 96.00,  3.1,1033.0, true },
  };
  printf("\n  %-24s %8s %10s %8s %8s %9s\n","fiducial volume","S","B","effS[%]","rejB[%]","FoM rel");
  printf("  %s\n", std::string(72,'-').c_str());
  double ref=-1;
  for(auto& c : cands){
    TString cut=box(c.x0,c.x1,c.y0,c.y1,c.z0,c.z1,c.dead);
    double s=S->GetEntries(TString::Format("(%s) && (%s)",SIGSEL,cut.Data()));
    double b=B->GetEntries(cut);
    if(s+b<=0) continue;
    double fom=s/std::sqrt(s+b);
    if(ref<0) ref=fom;
    printf("  %-24s %8.0f %10.0f %8.1f %8.1f %9.3f\n",
           c.nm, s, b, 100*s/ns0, 100*(1-b/nb0), fom/ref);
  }
  printf("\n  (FoM relative to the current box; >1 is better)\n");
  fs->Close(); fb->Close();
}
