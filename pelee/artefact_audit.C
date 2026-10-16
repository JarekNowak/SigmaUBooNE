// artefact_audit.C -- check every BDT input for production-chain artefacts.
//
// A variable that separates the SAMPLES (hyperon signal production vs standard NuMI
// overlay) rather than the PHYSICS will be seized on by the BDT and will invert on
// real data. slnunhits did exactly this: 0.774 separation between samples, purely
// because the branch is not populated in the hyperon production.
//
// The test compares each variable between the two samples with NO selection applied,
// so any difference is provenance, not physics: both are neutrino interactions in the
// same detector. Overlap below ~50% or separation above ~0.2 is a flag.
//
// Caveat, stated because it limits the conclusion: the two samples do differ
// physically -- one is hyperon-enriched, the other is inclusive NuMI RHC -- so some
// genuine difference is expected in variables that track the final state (track
// multiplicity, PID). The flag is a prompt to look, not a verdict. What it catches
// cleanly is the pathological case: near-disjoint distributions.
//   usage: root -l -b -q 'artefact_audit.C'
#include <vector>
#include <string>
#include <algorithm>

void artefact_audit(){
  const char* HY="/data/uboone/prod_nuwro_sigmazero_overlay_numi_rhc_run3b_extunbiased_mcc9/ana/Hyperons_Nuwro_run3b_WithCR_pelee_2500.root";
  const char* OV="/data/uboone/new_numi_flux/Run1_rhc_new_numi_flux_rhc_pandora_ntuple.root";

  // Every ntuple branch the selection reads that feeds a BDT input, with the range to
  // histogram over. Per-PFP vectors are summarised by their mean over the event.
  struct V { const char* expr; const char* nm; int nb; double lo, hi; };
  std::vector<V> vs = {
    {"topological_score","topological_score",100,0,1},
    {"slnhits","slnhits",100,0,6000},
    {"slnunhits","slnunhits (known bad)",100,0,3000},
    {"slclustfrac","slclustfrac",100,0,1.05},
    {"secondshower_Y_nhit","secondshower_Y_nhit",100,0,400},
    {"secondshower_Y_dot","secondshower_Y_dot",100,-1.1,1.1},
    {"shrsubclusters0+shrsubclusters1+shrsubclusters2","shrsubclusters (sum)",100,0,60},
    {"Sum$(trk_score_v)/Length$(trk_score_v)","trk_score_v (mean)",100,0,1},
    {"Sum$(trk_len_v)/Length$(trk_len_v)","trk_len_v (mean)",100,0,200},
    {"Sum$(trk_llr_pid_score_v)/Length$(trk_llr_pid_score_v)","trk_llr_pid (mean)",100,-1,1},
    {"Sum$(trk_pid_chipr_v)/Length$(trk_pid_chipr_v)","trk_pid_chipr (mean)",100,0,400},
    {"Sum$(trk_pid_chipi_v)/Length$(trk_pid_chipi_v)","trk_pid_chipi (mean)",100,0,400},
    {"Sum$(trk_energy_proton_v)/Length$(trk_energy_proton_v)","trk_energy_proton (mean)",100,0,1},
    {"Sum$(trk_range_muon_mom_v)/Length$(trk_range_muon_mom_v)","trk_range_muon_mom (mean)",100,0,2},
    {"Sum$(shr_energy_y_v)/Length$(shr_energy_y_v)","shr_energy_y_v (mean)",100,0,1000},
    {"Sum$(shr_moliere_avg_v)/Length$(shr_moliere_avg_v)","shr_moliere_avg (mean)",100,0,180},
    {"Sum$(shr_tkfit_dedx_y_v)/Length$(shr_tkfit_dedx_y_v)","shr_tkfit_dedx_y (mean)",100,0,15},
    {"Sum$(pfnplanehits_Y)/Length$(pfnplanehits_Y)","pfnplanehits_Y (mean)",100,0,600},
    {"Sum$(pfnhits)/Length$(pfnhits)","pfnhits (mean)",100,0,800},
    {"pi0_energy1_Y","pi0_energy1_Y",100,0,600},
    {"pi0_energy2_Y","pi0_energy2_Y",100,0,600},
    {"pi0_gammadot","pi0_gammadot",100,-1.1,1.1},
    {"n_tracks_contained","n_tracks_contained",12,0,12},
    {"n_showers_contained","n_showers_contained",10,0,10},
    {"reco_nu_vtx_sce_x","reco_nu_vtx_sce_x",100,-50,320},
    {"reco_nu_vtx_sce_y","reco_nu_vtx_sce_y",100,-150,150},
    {"reco_nu_vtx_sce_z","reco_nu_vtx_sce_z",100,-50,1100},
  };

  TFile* fh=TFile::Open(HY); TTree* th=(TTree*)fh->Get("nuselection/NeutrinoSelectionFilter");
  TFile* fo=TFile::Open(OV); TTree* to=(TTree*)fo->Get("nuselection/NeutrinoSelectionFilter");
  printf("\n  hyperon %lld entries, overlay %lld entries\n", th->GetEntries(), to->GetEntries());
  printf("\n  %-30s %9s %11s   %s\n","variable","overlap","separation","");
  printf("  %s\n", std::string(66,'-').c_str());

  struct R { double ov, sep; std::string nm; };
  std::vector<R> res;
  for (auto& v : vs) {
    if (!th->GetBranch(TString(v.expr).Contains("$")?"trk_score_v":v.expr) && !TString(v.expr).Contains("$")) {
      printf("  %-30s  (branch missing)\n", v.nm); continue;
    }
    TH1D a("a","",v.nb,v.lo,v.hi), b("b","",v.nb,v.lo,v.hi);
    th->Draw(Form("%s>>a",v.expr),"","goff");
    to->Draw(Form("%s>>b",v.expr),"","goff");
    if (a.Integral()<=0 || b.Integral()<=0) { printf("  %-30s  (empty)\n", v.nm); continue; }
    a.Scale(1./a.Integral()); b.Scale(1./b.Integral());
    double ov=0, sep=0;
    for (int i=1;i<=v.nb;++i){ double x=a.GetBinContent(i), y=b.GetBinContent(i);
      ov += std::min(x,y); if (x+y>0) sep += 0.5*(x-y)*(x-y)/(x+y); }
    res.push_back({ov,sep,v.nm});
  }
  std::sort(res.begin(),res.end(),[](const R&x,const R&y){return x.sep>y.sep;});
  for (auto& r : res) {
    const char* flag = (r.sep>0.5) ? "  <<< PATHOLOGICAL"
                     : (r.sep>0.2) ? "  <<  check"
                     : "";
    printf("  %-30s %8.1f%% %11.4f%s\n", r.nm.c_str(), 100*r.ov, r.sep, flag);
  }
  fh->Close(); fo->Close();
}
