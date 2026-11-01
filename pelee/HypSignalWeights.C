// ============================================================================
//  HypSignalWeights.C -- per-event hyperon weights for a GENIE signal production,
//  keyed on (run, subrun, event), for chains that cannot form them themselves
//  (the hypana BDT trees, PeakBDT.C).
//
//  w = w_xsec * w_pol:  w_xsec from HypXsecWeight.h (GENIE Pais -> full trace, struck
//  nucleon at rest, E_nu and Q2 from the PeLEE truth), w_pol from PolWeights.C's table.
//  Output "run sub evt w w_xsec w_pol" (FluxWeights.h reads the first weight).
//
//  Run:  root -l -b -q 'HypSignalWeights.C+("<pelee.root>","pol_genie_qe.tsv","hypw_genie_qe.tsv")'
// ============================================================================
#include "HypXsecWeight.h"
#include "FluxWeights.h"
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <vector>

void HypSignalWeights(const char* pelee, const char* polTable, const char* out,
                      const char* tables = "hyp_xsec_weights.root")
{
  HypXsecWeight::Table xw; if (!xw.load(tables)) return;
  const auto pol = loadFluxWeights(polTable);
  TFile f(pelee); TTree* t = (TTree*)f.Get("nuselection/NeutrinoSelectionFilter");
  if (!t) { printf("  cannot read %s\n", pelee); return; }
  Int_t run, sub, evt, nupdg, ccnc; Float_t e, px, py, pz;
  std::vector<int>* pd = nullptr; std::vector<float> *mx = nullptr, *my = nullptr, *mz = nullptr;
  t->SetBranchStatus("*", 0);
  for (auto b : {"run", "sub", "evt", "nu_pdg", "ccnc", "nu_e", "true_nu_px", "true_nu_py", "true_nu_pz",
                 "mc_pdg", "mc_px", "mc_py", "mc_pz"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("run", &run); t->SetBranchAddress("sub", &sub); t->SetBranchAddress("evt", &evt);
  t->SetBranchAddress("nu_pdg", &nupdg); t->SetBranchAddress("ccnc", &ccnc); t->SetBranchAddress("nu_e", &e);
  t->SetBranchAddress("true_nu_px", &px); t->SetBranchAddress("true_nu_py", &py); t->SetBranchAddress("true_nu_pz", &pz);
  t->SetBranchAddress("mc_pdg", &pd); t->SetBranchAddress("mc_px", &mx); t->SetBranchAddress("mc_py", &my); t->SetBranchAddress("mc_pz", &mz);
  std::ofstream fo(out);
  fo << "# run sub evt w w_xsec w_pol  -- HypSignalWeights.C (" << pelee << ", " << polTable << ")\n";
  long n = 0; double sx = 0, sp = 0;
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    if (nupdg != -14 || ccnc != 0) continue;
    int km = -1, ky = -1;
    for (size_t k = 0; k < pd->size(); ++k) {
      if (pd->at(k) == -13 && km < 0) km = (int)k;
      if ((pd->at(k) == 3122 || pd->at(k) == 3212 || pd->at(k) == 3112) && ky < 0) ky = (int)k;
    }
    if (km < 0 || ky < 0) continue;
    TLorentzVector pn, pm, pnuc(0, 0, 0, HypXsecWeight::nucleonMass(pd->at(ky)));
    pn.SetPxPyPzE(px, py, pz, e);
    pm.SetXYZM(mx->at(km), my->at(km), mz->at(km), HypXsecWeight::MMU);
    const double wx = xw.weight(pd->at(ky), pn, pnuc, pm);
    auto ip = pol.find(rseKey(run, sub, evt));
    const double wp = ip != pol.end() ? ip->second : 1.0;
    fo << run << " " << sub << " " << evt << " " << wx * wp << " " << wx << " " << wp << "\n";
    n++; sx += wx; sp += wp;
  }
  printf("  %ld events, <w_xsec> %.4f, <w_pol> %.4f -> %s\n", n, n ? sx / n : 0., n ? sp / n : 0., out);
}
