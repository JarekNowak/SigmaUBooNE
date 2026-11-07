// ============================================================================
//  CvWeights.C -- per-event weights of an island cache from an old-GENIE, old-flux sample
//  (the hyperon-enriched detvar CV, cthorpe_make_hyperon_events_numi_rhc_run3b_detvar_reco2_CV):
//    w_flux  old -> current NuMI flux, from the reco2_flux_weights.C table
//            ("run sub evt nu_pdg E angle w"; E_nu is also taken from it)
//    w_xsec  GENIE Pais -> full trace for numubar CC QE Lambda / Sigma0 / Sigma- (HypXsecWeight.h,
//            struck nucleon at rest, E_nu from the flux table, Q2 = mc_nu_q2), 1 otherwise
//  -- the same weights the overlay background gets (bkg_hypqe.tsv, built from PeLEE truth).
//  Output "run sub evt w_flux w_xsec found" for every isl entry (found = 0: no flux entry, w_flux 1).
//  Run:  root -l -b -q 'CvWeights.C+("isl_dvcv_bkg.root","/data/uboone/detvars/hypana/dv_hyp_rhc_CV_fluxweights.tsv","w_dvcv_bkg.tsv")'
// ============================================================================
#include "HypXsecWeight.h"
#include "FluxWeights.h"
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <map>
#include <utility>
#include <cstdio>

void CvWeights(const char* isl, const char* fluxTable, const char* out,
               const char* tables = "hyp_xsec_weights.root")
{
  HypXsecWeight::Table xw; if (!xw.load(tables)) return;
  std::map<unsigned long long, std::pair<double, double>> fw;   // key -> (E_nu, w)
  { std::ifstream in(fluxTable); std::string l;
    while (std::getline(in, l)) {
      if (l.empty() || l[0] == '#') continue;
      unsigned r, s, e; int pdg; double E, ang, w;
      if (sscanf(l.c_str(), "%u %u %u %d %lf %lf %lf", &r, &s, &e, &pdg, &E, &ang, &w) == 7)
        fw[rseKey(r, s, e)] = {E, w};
    } }
  printf("  flux table: %zu events from %s\n", fw.size(), fluxTable);
  TFile f(isl); TTree* t = (TTree*)f.Get("isl");
  if (!t || !t->GetBranch("hypPdg")) { printf("  %s: no isl tree with hypPdg\n", isl); return; }
  int run, sub, evt, mode, hyp, nu, ccnc; double q2;
  t->SetBranchAddress("run", &run); t->SetBranchAddress("sub", &sub); t->SetBranchAddress("evt", &evt);
  t->SetBranchAddress("mode", &mode); t->SetBranchAddress("hypPdg", &hyp); t->SetBranchAddress("nuPdg", &nu);
  t->SetBranchAddress("ccnc", &ccnc); t->SetBranchAddress("q2", &q2);
  std::ofstream fo(out);
  fo << "# run sub evt w_flux w_xsec found  -- CvWeights.C (" << isl << ", " << fluxTable << ")\n";
  long n = 0, nMiss = 0, nQE = 0; double sf = 0, sx = 0;
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    auto it = fw.find(rseKey(run, sub, evt));
    const bool found = it != fw.end();
    const double wf = found ? it->second.second : 1.0, E = found ? it->second.first : -1;
    double wx = 1.0;
    if (nu < 0 && ccnc == 0 && mode == 0 && (hyp == 3122 || hyp == 3212 || hyp == 3112) && E > 0) {
      wx = xw.weight(hyp, E, q2); sx += wx; ++nQE; }
    fo << run << " " << sub << " " << evt << " " << wf << " " << wx << " " << found << "\n";
    ++n; if (!found) ++nMiss; else sf += wf;
  }
  printf("  %ld entries: flux <w> %.4f (%ld without a flux entry), %ld numubar CC QE hyperons <w_xsec> %.4f -> %s\n",
         n, (n - nMiss) ? sf / (n - nMiss) : 0., nMiss, nQE, nQE ? sx / nQE : 0., out);
}
